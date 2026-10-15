use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, AsteroidId, CometCloudId, CometId, DebrisDiskId, MoonId, PlanetId,
    ProtoplanetaryDiskId, StarId, StarSystemId,
};
use crate::orbital::{OrbitalElements, propagate_kepler};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration};
use crate::units::{
    Angle, Length, Luminosity, METERS_PER_ASTRONOMICAL_UNIT, Mass, Position, Temperature, Velocity,
};

use super::{
    Asteroid, AsteroidBelt, Comet, CometCloud, DebrisDisk, HabitableZoneHistory,
    HabitableZoneWindow, Moon, MoonFormationHistory, MoonOrigin, ObjectError, Planet,
    PlanetAtmosphere, PlanetStatus, PlanetaryComposition, ProtoplanetaryDisk, Star,
};

const HABITABLE_INNER_AU_AT_SOLAR_LUMINOSITY: f64 = 0.95;
const HABITABLE_OUTER_AU_AT_SOLAR_LUMINOSITY: f64 = 1.37;

pub(super) fn modeled_habitable_zone_history(
    planet: &Planet,
    star: &Star,
    through: SimDuration,
) -> Result<HabitableZoneHistory, ObjectError> {
    if planet.star_id != star.id {
        return Err(ObjectError("planet and star do not match".to_owned()));
    }
    if matches!(
        planet.status,
        PlanetStatus::Unbound | PlanetStatus::Destroyed
    ) {
        return Err(ObjectError(
            "zone residence is unavailable for unbound or destroyed planets".to_owned(),
        ));
    }
    let start = planet.formed_at.seconds().max(star.formed_at.seconds());
    let end = through.seconds();
    let mut history = HabitableZoneHistory {
        windows: Vec::new(),
        total_duration: SimDuration::ZERO,
    };
    if end < start {
        return Ok(history);
    }

    let mut boundaries = vec![start, end];
    for (age, _) in super::stellar::stellar_schedule(star.initial_mass, star.metallicity)
        .map_err(|error| ObjectError(error.to_string()))?
    {
        if let Some(at) = star
            .formed_at
            .seconds()
            .checked_add(age.seconds().saturating_sub(star.epoch_age.seconds()))
            && (start..end).contains(&at)
        {
            boundaries.push(at);
        }
    }
    for limit in [0.03, 10_000.0] {
        if planet.migration_rate_au_per_year == 0.0 {
            break;
        }
        let years = (limit - planet.initial_orbital_radius.astronomical_units())
            / planet.migration_rate_au_per_year;
        let seconds = years * SECONDS_PER_YEAR as f64;
        if seconds.is_finite() && seconds >= 0.0 && seconds < u128::MAX as f64 {
            let offset = seconds.round() as u128;
            if let Some(at) = planet.formed_at.seconds().checked_add(offset)
                && (start..end).contains(&at)
            {
                boundaries.push(at);
            }
        }
    }
    boundaries.sort_unstable();
    boundaries.dedup();

    for pair in boundaries.windows(2) {
        let (segment_start, segment_end) = (pair[0], pair[1]);
        let midpoint = segment_start + (segment_end - segment_start) / 2;
        let stage = modeled_stellar_stage(star, midpoint)?;
        let luminosity_start = modeled_luminosity(star, stage, segment_start)?;
        let luminosity_end = modeled_luminosity(star, stage, segment_end)?;
        let orbit_start = modeled_orbit_au(planet, segment_start);
        let orbit_end = modeled_orbit_au(planet, segment_end);
        let mut cuts = vec![segment_start, segment_end];
        add_habitable_crossings(
            &mut cuts,
            segment_start,
            segment_end,
            orbit_start,
            orbit_end,
            luminosity_start,
            luminosity_end,
        );
        cuts.sort_unstable();
        cuts.dedup();
        for interval in cuts.windows(2) {
            let at = interval[0] + (interval[1] - interval[0]) / 2;
            let orbit = modeled_orbit_au(planet, at);
            let luminosity = modeled_luminosity(star, stage, at)?;
            if within_habitable_zone(luminosity, orbit) {
                if let Some(last) = history.windows.last_mut()
                    && last.exited_at == Some(SimDuration::from_seconds(interval[0]))
                {
                    last.exited_at = Some(SimDuration::from_seconds(interval[1]));
                } else {
                    history.windows.push(HabitableZoneWindow {
                        entered_at: SimDuration::from_seconds(interval[0]),
                        exited_at: Some(SimDuration::from_seconds(interval[1])),
                    });
                }
                history.total_duration = history
                    .total_duration
                    .checked_add(SimDuration::from_seconds(interval[1] - interval[0]))
                    .ok_or_else(|| ObjectError("habitable-zone duration overflow".to_owned()))?;
            }
        }
    }
    let at_end = modeled_stellar_stage(star, end)?;
    let end_luminosity = modeled_luminosity(star, at_end, end)?;
    let end_orbit = modeled_orbit_au(planet, end);
    let currently_inside = within_habitable_zone(end_luminosity, end_orbit);
    if currently_inside {
        if let Some(last) = history.windows.last_mut()
            && last.exited_at == Some(through)
        {
            last.exited_at = None;
        } else if start == end {
            history.windows.push(HabitableZoneWindow {
                entered_at: through,
                exited_at: None,
            });
        }
    }
    Ok(history)
}

fn modeled_stellar_stage(star: &Star, at: u128) -> Result<super::StellarStage, ObjectError> {
    let age = star
        .epoch_age
        .checked_add(SimDuration::from_seconds(
            at.saturating_sub(star.formed_at.seconds()),
        ))
        .ok_or_else(|| ObjectError("stellar age overflow".to_owned()))?;
    super::stellar::stellar_stage_at(star.initial_mass, star.metallicity, age)
}

fn modeled_luminosity(
    star: &Star,
    stage: super::StellarStage,
    at: u128,
) -> Result<f64, ObjectError> {
    let age = star
        .epoch_age
        .checked_add(SimDuration::from_seconds(
            at.saturating_sub(star.formed_at.seconds()),
        ))
        .ok_or_else(|| ObjectError("stellar age overflow".to_owned()))?;
    Ok(
        super::stellar::stellar_properties(star.initial_mass, star.metallicity, stage, age)?
            .1
            .solar_luminosities(),
    )
}

fn modeled_orbit_au(planet: &Planet, at: u128) -> f64 {
    let years = at.saturating_sub(planet.formed_at.seconds()) as f64 / SECONDS_PER_YEAR as f64;
    (planet.initial_orbital_radius.astronomical_units() + planet.migration_rate_au_per_year * years)
        .clamp(0.03, 10_000.0)
}

fn add_habitable_crossings(
    cuts: &mut Vec<u128>,
    start: u128,
    end: u128,
    orbit_start: f64,
    orbit_end: f64,
    luminosity_start: f64,
    luminosity_end: f64,
) {
    let delta_orbit = orbit_end - orbit_start;
    let delta_luminosity = luminosity_end - luminosity_start;
    for edge in [
        HABITABLE_INNER_AU_AT_SOLAR_LUMINOSITY,
        HABITABLE_OUTER_AU_AT_SOLAR_LUMINOSITY,
    ] {
        let a = delta_orbit * delta_orbit;
        let b = 2.0 * orbit_start * delta_orbit - edge * edge * delta_luminosity;
        let c = orbit_start * orbit_start - edge * edge * luminosity_start;
        let roots = if a == 0.0 {
            if b == 0.0 { vec![] } else { vec![-c / b] }
        } else {
            let discriminant = b.mul_add(b, -4.0 * a * c);
            if discriminant < 0.0 {
                vec![]
            } else {
                let square_root = discriminant.sqrt();
                let q = -0.5 * (b + square_root.copysign(b));
                if q == 0.0 {
                    vec![-b / (2.0 * a)]
                } else {
                    vec![q / a, c / q]
                }
            }
        };
        for root in roots {
            if root > 0.0 && root < 1.0 {
                let offset = ((end - start) as f64 * root).round() as u128;
                if offset < end - start {
                    cuts.push(start + offset);
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn generate_planet(
    seed: u64,
    id: PlanetId,
    system_id: StarSystemId,
    star_id: StarId,
    planet_position: Position,
    planet_velocity: Velocity,
    star_mass: Mass,
    star_luminosity: Luminosity,
    orbit_au: f64,
    formation_material: GasComposition,
    formed_at: SimDuration,
    config: GenerationConfig,
) -> Result<Planet, ObjectError> {
    let mut rng = StableRng::new(derive_object_seed(seed, id.value(), "planet-properties-v1"));
    let metallicity = formation_material.metallicity();
    let metal_mass_bias = (metallicity / 0.014).clamp(0.25, 2.0);
    let mass_min = config.planet_mass_earth_min.ln();
    let mass_max = config.planet_mass_earth_max.ln();
    let mass_earth = (mass_min + rng.next_f64().powf(1.5) * (mass_max - mass_min))
        .exp()
        .mul_add(metal_mass_bias.sqrt(), 0.0)
        .clamp(config.planet_mass_earth_min, config.planet_mass_earth_max);
    let orbital_radius = Length::from_astronomical_units(orbit_au)?;
    let atmosphere = generate_planet_atmosphere(&mut rng, mass_earth, metallicity);
    let radius = Length::from_meters(planet_radius_meters(
        Mass::from_earth_masses(mass_earth)?,
        atmosphere,
    ))?;
    let equilibrium_temperature =
        Temperature::from_kelvin(equilibrium_temperature_kelvin(star_luminosity, orbit_au))?;
    let bulk_composition = planetary_composition(formation_material, orbit_au, mass_earth);
    let orbital_elements = OrbitalElements::new(
        orbital_radius,
        0.25 * rng.next_f64().powi(2),
        Angle::from_degrees(super::random_signed(&mut rng) * 5.0)?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        SimDuration::ZERO,
    )?;
    let (position, velocity) = propagate_kepler(
        orbital_elements,
        star_mass,
        SimDuration::ZERO,
        planet_position,
        planet_velocity,
    )?;
    Ok(Planet {
        id,
        system_id,
        star_id,
        mass: Mass::from_earth_masses(mass_earth)?,
        radius,
        atmosphere,
        atmosphere_loss_fraction: 0.0,
        equilibrium_temperature,
        in_habitable_zone: is_in_habitable_zone(star_luminosity, orbit_au),
        bulk_composition,
        formed_at,
        initial_orbital_radius: orbital_radius,
        migration_rate_au_per_year: super::random_signed(&mut rng) * 2.0e-10,
        scattering_count: 0,
        orbital_radius,
        position,
        velocity,
        orbital_elements,
        status: PlanetStatus::Stable,
    })
}

pub(super) fn generate_moon_count(rng: &mut StableRng, planet_mass: Mass) -> u8 {
    let mass_earth = planet_mass.earth_masses();
    if mass_earth < 0.2 {
        return 0;
    }
    let occurrence = (0.10 + 0.25 * mass_earth.ln().max(0.0)).clamp(0.10, 0.85);
    if rng.next_f64() >= occurrence {
        return 0;
    }
    let maximum: u64 = if mass_earth >= 20.0 {
        4
    } else if mass_earth >= 2.0 {
        2
    } else {
        1
    };
    1 + u8::try_from(rng.next_u64_below(maximum).unwrap_or(0)).unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn generate_moon(
    id: MoonId,
    planet_id: PlanetId,
    system_id: StarSystemId,
    planet_position: Position,
    planet_velocity: Velocity,
    planet_mass: Mass,
    star_mass: Mass,
    star_luminosity: Luminosity,
    planet_orbit_au: f64,
    formed_at: SimDuration,
    rng: &mut StableRng,
) -> Result<Option<Moon>, ObjectError> {
    let planet_radius = planet_radius_meters(planet_mass, PlanetAtmosphere::None);
    let hill_radius = planet_orbit_au
        * METERS_PER_ASTRONOMICAL_UNIT
        * (planet_mass.kilograms() / (3.0 * star_mass.kilograms())).cbrt();
    let minimum_axis = planet_radius * 6.0;
    let maximum_axis = hill_radius * 0.25;
    if !maximum_axis.is_finite() || maximum_axis <= minimum_axis {
        return Ok(None);
    }
    let axis = super::random_range(rng, minimum_axis, maximum_axis);
    let mass_earth = (super::random_range(rng, 0.00001, 0.02) * planet_mass.earth_masses().sqrt())
        .clamp(0.000001, 0.08);
    let orbital_elements = OrbitalElements::new(
        Length::from_meters(axis)?,
        super::random_range(rng, 0.0, 0.15),
        Angle::from_degrees(super::random_signed(rng) * 8.0)?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        SimDuration::ZERO,
    )?;
    let moon_mass = Mass::from_earth_masses(mass_earth)?;
    let moon_atmosphere = if mass_earth >= 0.01 {
        if rng.next_f64() < 0.5 {
            PlanetAtmosphere::Thin
        } else {
            PlanetAtmosphere::None
        }
    } else {
        PlanetAtmosphere::None
    };
    let origin = if rng.next_f64() < 0.15 {
        MoonOrigin::Captured
    } else if rng.next_f64() < 0.30 {
        MoonOrigin::Impact
    } else {
        MoonOrigin::InSitu
    };
    let tidal_rate = (planet_mass.earth_masses().sqrt() / mass_earth.max(1e-6)).clamp(0.001, 0.1)
        * if origin == MoonOrigin::Captured {
            -1.0
        } else {
            1.0
        };
    let (position, velocity) = propagate_kepler(
        orbital_elements,
        planet_mass,
        SimDuration::ZERO,
        planet_position,
        planet_velocity,
    )?;
    Ok(Some(Moon {
        id,
        planet_id,
        system_id,
        mass: moon_mass,
        radius: Length::from_meters(planet_radius_meters(moon_mass, PlanetAtmosphere::None))?,
        atmosphere: moon_atmosphere,
        equilibrium_temperature: Temperature::from_kelvin(equilibrium_temperature_kelvin(
            star_luminosity,
            planet_orbit_au,
        ))?,
        formation_history: MoonFormationHistory::new(
            origin,
            formed_at,
            (origin == MoonOrigin::Captured).then_some(formed_at),
        ),
        tidal_migration_rate_m_per_year: tidal_rate,
        initial_orbital_radius: Length::from_meters(axis)?,
        orbital_radius: Length::from_meters(axis)?,
        position,
        velocity,
        orbital_elements,
    }))
}

pub(super) fn generate_asteroids(
    seed: u64,
    belt_id: AsteroidBeltId,
    system_id: StarSystemId,
    star_mass: Mass,
    origin: Position,
    origin_velocity: Velocity,
    belt: &AsteroidBelt,
) -> Result<Vec<Asteroid>, ObjectError> {
    (0..3)
        .map(|index| {
            let id = AsteroidId::new(derive_object_seed(
                seed,
                belt_id.value() ^ index,
                "asteroid-member-id-v1",
            ));
            let mut rng = StableRng::new(derive_object_seed(seed, id.value(), "asteroid-v1"));
            let orbit_au = super::random_range(
                &mut rng,
                belt.inner_radius.astronomical_units(),
                belt.outer_radius.astronomical_units(),
            );
            let elements = OrbitalElements::new(
                Length::from_astronomical_units(orbit_au)?,
                super::random_range(&mut rng, 0.0, 0.25),
                Angle::from_degrees(super::random_signed(&mut rng) * 10.0)?,
                Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
                SimDuration::ZERO,
            )?;
            let (position, velocity) = propagate_kepler(
                elements,
                star_mass,
                SimDuration::ZERO,
                origin,
                origin_velocity,
            )?;
            Ok(Asteroid {
                id,
                belt_id,
                system_id,
                mass: Mass::from_kilograms(super::random_range(&mut rng, 1.0e12, 1.0e18))?,
                radius: Length::from_meters(super::random_range(&mut rng, 500.0, 100_000.0))?,
                position,
                velocity,
                orbital_elements: elements,
            })
        })
        .collect()
}

pub(super) fn generate_comets(
    seed: u64,
    cloud_id: CometCloudId,
    system_id: StarSystemId,
    star_mass: Mass,
    origin: Position,
    origin_velocity: Velocity,
    cloud: &CometCloud,
) -> Result<Vec<Comet>, ObjectError> {
    (0..3)
        .map(|index| {
            let id = CometId::new(derive_object_seed(
                seed,
                cloud_id.value() ^ index,
                "comet-member-id-v1",
            ));
            let mut rng = StableRng::new(derive_object_seed(seed, id.value(), "comet-v1"));
            let orbit_au = super::random_range(
                &mut rng,
                cloud.inner_radius.astronomical_units(),
                cloud.outer_radius.astronomical_units(),
            );
            let elements = OrbitalElements::new(
                Length::from_astronomical_units(orbit_au)?,
                super::random_range(&mut rng, 0.1, 0.85),
                Angle::from_degrees(super::random_signed(&mut rng) * 30.0)?,
                Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
                SimDuration::ZERO,
            )?;
            let (position, velocity) = propagate_kepler(
                elements,
                star_mass,
                SimDuration::ZERO,
                origin,
                origin_velocity,
            )?;
            Ok(Comet {
                id,
                cloud_id,
                system_id,
                mass: Mass::from_kilograms(super::random_range(&mut rng, 1.0e11, 1.0e16))?,
                radius: Length::from_meters(super::random_range(&mut rng, 100.0, 50_000.0))?,
                volatile_fraction: super::random_range(&mut rng, 0.2, 0.9),
                position,
                velocity,
                orbital_elements: elements,
            })
        })
        .collect()
}

pub(super) fn generate_protoplanetary_disk(
    id: ProtoplanetaryDiskId,
    system_id: StarSystemId,
    star_mass: Mass,
    star_age: SimDuration,
    rng: &mut StableRng,
) -> Result<ProtoplanetaryDisk, ObjectError> {
    let inner_radius = Length::from_astronomical_units(super::random_range(rng, 0.05, 0.2))?;
    let outer_radius = Length::from_astronomical_units(super::random_range(rng, 20.0, 80.0))?;
    let initial_mass = Mass::from_solar_masses(
        (0.01 * star_mass.solar_masses() * super::random_range(rng, 0.5, 2.0)).clamp(0.001, 0.2),
    )?;
    let dissipates_at =
        SimDuration::from_years(10_000_000).map_err(|error| ObjectError(error.to_string()))?;
    Ok(ProtoplanetaryDisk {
        id,
        system_id,
        inner_radius,
        outer_radius,
        initial_mass,
        remaining_mass: if star_age < dissipates_at {
            initial_mass
        } else {
            Mass::from_kilograms(0.0)?
        },
        formed_at: SimDuration::ZERO,
        dissipates_at,
    })
}

pub(super) fn generate_asteroid_belt(
    id: AsteroidBeltId,
    system_id: StarSystemId,
    rng: &mut StableRng,
    max_planet_orbit_au: f64,
) -> Result<AsteroidBelt, ObjectError> {
    let inner_au = if max_planet_orbit_au > 0.0 {
        max_planet_orbit_au * super::random_range(rng, 1.25, 2.0)
    } else {
        super::random_range(rng, 1.0, 3.0)
    };
    let outer_au = inner_au * super::random_range(rng, 1.4, 3.0);
    Ok(AsteroidBelt {
        id,
        system_id,
        inner_radius: Length::from_astronomical_units(inner_au)?,
        outer_radius: Length::from_astronomical_units(outer_au)?,
        estimated_body_count: 1_000_000 + rng.next_u64_below(9_000_000_000).unwrap_or(0),
        total_mass: Mass::from_earth_masses(super::random_range(rng, 0.00001, 0.25))?,
    })
}

pub(super) fn generate_comet_cloud(
    id: CometCloudId,
    system_id: StarSystemId,
    rng: &mut StableRng,
) -> Result<CometCloud, ObjectError> {
    let inner_au = super::random_range(rng, 30.0, 80.0);
    let outer_au = inner_au * super::random_range(rng, 4.0, 20.0);
    Ok(CometCloud {
        id,
        system_id,
        inner_radius: Length::from_astronomical_units(inner_au)?,
        outer_radius: Length::from_astronomical_units(outer_au)?,
        estimated_body_count: 100_000 + rng.next_u64_below(10_000_000_000).unwrap_or(0),
        total_mass: Mass::from_earth_masses(super::random_range(rng, 0.001, 2.0))?,
    })
}

pub(super) fn generate_debris_disk(
    id: DebrisDiskId,
    system_id: StarSystemId,
    rng: &mut StableRng,
) -> Result<DebrisDisk, ObjectError> {
    let inner_au = super::random_range(rng, 0.05, 0.5);
    let outer_au = inner_au * super::random_range(rng, 5.0, 30.0);
    Ok(DebrisDisk {
        id,
        system_id,
        inner_radius: Length::from_astronomical_units(inner_au)?,
        outer_radius: Length::from_astronomical_units(outer_au)?,
        estimated_body_count: 10_000 + rng.next_u64_below(100_000_000).unwrap_or(0),
        total_mass: Mass::from_earth_masses(super::random_range(rng, 0.000001, 0.05))?,
    })
}

fn generate_planet_atmosphere(
    rng: &mut StableRng,
    mass_earth: f64,
    metallicity: f64,
) -> PlanetAtmosphere {
    if mass_earth < 0.01 {
        PlanetAtmosphere::None
    } else if mass_earth >= 15.0 {
        PlanetAtmosphere::HydrogenHelium
    } else if mass_earth < 0.1 || rng.next_f64() > (0.45 + metallicity * 4.0).clamp(0.45, 0.90) {
        PlanetAtmosphere::Thin
    } else if rng.next_f64() < 0.70 {
        PlanetAtmosphere::NitrogenRich
    } else {
        PlanetAtmosphere::CarbonDioxideRich
    }
}

pub(super) fn planet_radius_meters(mass: Mass, atmosphere: PlanetAtmosphere) -> f64 {
    let mass_earth = mass.earth_masses().max(0.000001);
    if atmosphere == PlanetAtmosphere::HydrogenHelium {
        11.2e6 * (mass_earth / 318.0).powf(0.01).clamp(0.75, 1.15)
    } else {
        6.371e6 * mass_earth.powf(0.27)
    }
}

pub(super) fn equilibrium_temperature_kelvin(luminosity: Luminosity, orbit_au: f64) -> f64 {
    let luminosity_solar = luminosity.solar_luminosities().max(0.0);
    (278.5 * luminosity_solar.powf(0.25) / orbit_au.max(0.001).sqrt()).clamp(2.7, 10_000.0)
}

pub(super) fn is_in_habitable_zone(luminosity: Luminosity, orbit_au: f64) -> bool {
    within_habitable_zone(luminosity.solar_luminosities(), orbit_au)
}

fn within_habitable_zone(luminosity_solar: f64, orbit_au: f64) -> bool {
    let root = luminosity_solar.max(0.0).sqrt();
    (HABITABLE_INNER_AU_AT_SOLAR_LUMINOSITY * root..=HABITABLE_OUTER_AU_AT_SOLAR_LUMINOSITY * root)
        .contains(&orbit_au)
}

fn planetary_composition(
    material: GasComposition,
    orbit_au: f64,
    mass_earth: f64,
) -> PlanetaryComposition {
    let metals = material.metallicity().max(1e-9);
    let iron = (0.05 + 0.20 * material.iron_fraction() / metals).clamp(0.03, 0.30);
    let icy = if orbit_au >= 3.0 {
        (0.08 + 2.0 * metals).clamp(0.08, 0.45)
    } else {
        0.02
    };
    let hydrogen_helium = if mass_earth >= 15.0 {
        (0.65 + 0.1 * (1.0 - metals / 0.04)).clamp(0.55, 0.80)
    } else {
        0.0
    };
    let rocky = (1.0 - icy - iron - hydrogen_helium).max(0.0);
    let total = rocky + icy + iron + hydrogen_helium;
    PlanetaryComposition::new(
        rocky / total,
        icy / total,
        iron / total,
        hydrogen_helium / total,
    )
}

pub(super) fn validate_small_body_population(
    inner_radius: Length,
    outer_radius: Length,
    estimated_body_count: u64,
    kind: &str,
) -> Result<(), ObjectError> {
    if inner_radius >= outer_radius || estimated_body_count == 0 {
        return Err(ObjectError(format!(
            "{kind} must have an increasing radial range and a nonzero population"
        )));
    }
    Ok(())
}

pub(super) fn generate_planet_count(
    rng: &mut StableRng,
    config: GenerationConfig,
    star_mass: Mass,
    metallicity: f64,
) -> u8 {
    let occurrence = (0.25
        + 0.35 * (metallicity / 0.014).clamp(0.0, 2.0)
        + 0.15 * (1.0 / star_mass.solar_masses()).clamp(0.0, 2.0))
    .clamp(0.05, 0.95);
    if config.planets_per_system_min == 0 && rng.next_f64() > occurrence {
        return 0;
    }
    let span = u64::from(config.planets_per_system_max - config.planets_per_system_min + 1);
    config.planets_per_system_min
        + u8::try_from(
            rng.next_u64_below(span)
                .expect("planet count span is nonzero"),
        )
        .expect("planet count span fits in u8")
}

pub(super) fn stable_orbital_radii(
    rng: &mut StableRng,
    count: u8,
    config: GenerationConfig,
) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }
    if count == 1 {
        return vec![
            (config.orbit_au_min.ln()
                + rng.next_f64() * (config.orbit_au_max.ln() - config.orbit_au_min.ln()))
            .exp(),
        ];
    }
    let ratio = (config.orbit_au_max / config.orbit_au_min).powf(1.0 / f64::from(count - 1));
    (0..count)
        .map(|index| config.orbit_au_min * ratio.powi(i32::from(index)))
        .collect()
}
