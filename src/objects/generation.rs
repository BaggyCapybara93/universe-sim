use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, CometCloudId, DebrisDiskId, MoonId, PlanetId, ProtoplanetaryDiskId, StarId,
    StarSystemId,
};
use crate::orbital::{OrbitalElements, propagate_kepler};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration};
use crate::units::{Angle, Length, METERS_PER_PARSEC, Mass, Position, Velocity};

use super::{
    Asteroid, AsteroidBelt, Comet, CometCloud, DebrisDisk, GalacticComponent, GalaxyCatalogEntry,
    GalaxyMorphology, Moon, ObjectError, Planet, ProtoplanetaryDisk, Star, StarSystem,
    StellarStage,
};

pub(super) struct GeneratedSystem {
    pub(super) system: StarSystem,
    pub(super) stars: Vec<Star>,
    pub(super) planets: Vec<Planet>,
    pub(super) moons: Vec<Moon>,
    pub(super) asteroids: Vec<Asteroid>,
    pub(super) comets: Vec<Comet>,
    pub(super) asteroid_belts: Vec<AsteroidBelt>,
    pub(super) comet_clouds: Vec<CometCloud>,
    pub(super) debris_disks: Vec<DebrisDisk>,
    pub(super) protoplanetary_disks: Vec<ProtoplanetaryDisk>,
}

pub(super) fn generate_system(
    seed: u64,
    config: GenerationConfig,
    galaxy: GalaxyCatalogEntry,
    system_index: u32,
) -> Result<GeneratedSystem, ObjectError> {
    let galaxy_id = galaxy.id;
    let system_key = u64::from(system_index);
    let system_id = StarSystemId::new(derive_object_seed(seed, system_key, "star-system-id-v1"));
    let star_id = StarId::new(derive_object_seed(seed, system_key, "star-id-v1"));
    let mut system_rng = StableRng::new(derive_object_seed(
        seed,
        system_id.value(),
        "star-system-properties-v1",
    ));
    let population = generate_system_population(&mut system_rng, galaxy)?;
    let [gx, gy, gz] = galaxy.position.meters();
    let [rx, ry, rz] = population.position.meters();
    let position = Position::from_meters(gx + rx, gy + ry, gz + rz)?;
    let velocity = population.velocity;
    let mut primary = generate_star(
        seed,
        star_id,
        system_id,
        population.epoch_age,
        population.metallicity,
        config,
    )?;
    primary.position = position;
    primary.velocity = velocity;
    let primary_mass = primary.initial_mass;
    let current_primary_mass = primary.current_mass;
    let current_primary_luminosity = primary.luminosity;
    let primary_stage = primary.stage;
    let mut stars = vec![primary];

    let binary_probability =
        (0.20 + 0.15 * primary_mass.solar_masses().ln().max(0.0)).clamp(0.15, 0.70);
    let mut companion_star_ids = Vec::new();
    let mut inner_companion_axis_au = None;
    if system_rng.next_f64() < binary_probability {
        let triple_probability =
            (0.05 + 0.04 * primary_mass.solar_masses().sqrt()).clamp(0.05, 0.25);
        let companion_count = if system_rng.next_f64() < triple_probability {
            2
        } else {
            1
        };
        for companion_index in 0..companion_count {
            let namespace = format!("companion-star-id-v1-{companion_index}");
            let companion_id = StarId::new(derive_object_seed(seed, system_key, &namespace));
            let mass_ratio = super::random_range(&mut system_rng, 0.10, 1.0);
            let companion_mass_solar = (primary_mass.solar_masses() * mass_ratio)
                .clamp(config.star_mass_solar_min, config.star_mass_solar_max);
            let mut companion = generate_star_with_mass(
                companion_id,
                system_id,
                population.epoch_age,
                population.metallicity,
                companion_mass_solar,
            )?;
            let semi_major_axis_au = if let Some(inner_axis) = inner_companion_axis_au {
                inner_axis * super::random_range(&mut system_rng, 10.0, 30.0)
            } else {
                (super::random_range(&mut system_rng, 0.05_f64.ln(), 50.0_f64.ln())).exp()
            };
            inner_companion_axis_au.get_or_insert(semi_major_axis_au);
            let orbit = OrbitalElements::new(
                Length::from_astronomical_units(semi_major_axis_au)?,
                0.6 * system_rng.next_f64().powi(2),
                Angle::from_degrees(super::random_signed(&mut system_rng) * 30.0)?,
                Angle::from_radians(std::f64::consts::TAU * system_rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * system_rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * system_rng.next_f64())?,
                SimDuration::ZERO,
            )?;
            let interior_mass: f64 = stars.iter().map(|star| star.initial_mass.kilograms()).sum();
            let orbital_mass =
                Mass::from_kilograms(interior_mass + companion.initial_mass.kilograms())?;
            let (companion_position, companion_velocity) =
                propagate_kepler(orbit, orbital_mass, SimDuration::ZERO, position, velocity)?;
            companion.position = companion_position;
            companion.velocity = companion_velocity;
            companion.orbital_elements = Some(orbit);
            stars.push(companion);
            companion_star_ids.push(companion_id);
        }
    }

    let planet_count = if primary_stage == StellarStage::Disrupted {
        0
    } else {
        super::planetary::generate_planet_count(
            &mut system_rng,
            config,
            primary_mass,
            population.metallicity,
        )
    };
    let orbital_radii =
        super::planetary::stable_orbital_radii(&mut system_rng, planet_count, config);
    let max_planet_orbit_au = orbital_radii.iter().copied().fold(0.0_f64, f64::max);
    let mut planet_ids = Vec::with_capacity(usize::from(planet_count));
    let mut planets = Vec::with_capacity(usize::from(planet_count));
    let mut moon_ids = Vec::new();
    let mut moons = Vec::new();
    let formation_material = GasComposition::from_metallicity(population.metallicity)
        .map_err(|error| ObjectError(error.to_string()))?;
    for (planet_index, orbit_au) in (0..planet_count).zip(orbital_radii) {
        let planet_key = (system_key << 8) | u64::from(planet_index);
        let planet_id = PlanetId::new(derive_object_seed(seed, planet_key, "planet-id-v1"));
        let planet = super::planetary::generate_planet(
            seed,
            planet_id,
            system_id,
            star_id,
            position,
            velocity,
            current_primary_mass,
            current_primary_luminosity,
            orbit_au,
            formation_material,
            SimDuration::ZERO,
            config,
        )?;
        planet_ids.push(planet_id);
        let moon_count = super::planetary::generate_moon_count(&mut system_rng, planet.mass);
        for moon_index in 0..moon_count {
            let moon_key =
                (system_key << 16) | (u64::from(planet_index) << 8) | u64::from(moon_index);
            let moon_id = MoonId::new(derive_object_seed(seed, moon_key, "moon-id-v1"));
            if let Some(moon) = super::planetary::generate_moon(
                moon_id,
                planet.id,
                system_id,
                planet.position,
                planet.velocity,
                planet.mass,
                current_primary_mass,
                current_primary_luminosity,
                orbit_au,
                SimDuration::ZERO,
                &mut system_rng,
            )? {
                moons.push(moon);
                moon_ids.push(moon_id);
            }
        }
        planets.push(planet);
    }

    let mut asteroid_ids = Vec::new();
    let mut asteroids = Vec::new();

    let mut asteroid_belt_ids = Vec::new();
    let mut asteroid_belts = Vec::new();
    if primary_stage != StellarStage::Disrupted && system_rng.next_f64() < 0.70 {
        let id = AsteroidBeltId::new(derive_object_seed(seed, system_key, "asteroid-belt-id-v1"));
        let belt = super::planetary::generate_asteroid_belt(
            id,
            system_id,
            &mut system_rng,
            max_planet_orbit_au,
        )?;
        let members = super::planetary::generate_asteroids(
            seed,
            id,
            system_id,
            current_primary_mass,
            position,
            velocity,
            &belt,
        )?;
        asteroid_ids.extend(members.iter().map(|asteroid| asteroid.id));
        asteroids.extend(members);
        asteroid_belts.push(belt);
        asteroid_belt_ids.push(id);
    }
    let mut comet_cloud_ids = Vec::new();
    let mut comet_clouds = Vec::new();
    let mut comet_ids = Vec::new();
    let mut comets = Vec::new();
    if primary_stage != StellarStage::Disrupted && system_rng.next_f64() < 0.55 {
        let id = CometCloudId::new(derive_object_seed(seed, system_key, "comet-cloud-id-v1"));
        let cloud = super::planetary::generate_comet_cloud(id, system_id, &mut system_rng)?;
        let members = super::planetary::generate_comets(
            seed,
            id,
            system_id,
            current_primary_mass,
            position,
            velocity,
            &cloud,
        )?;
        comet_ids.extend(members.iter().map(|comet| comet.id));
        comets.extend(members);
        comet_clouds.push(cloud);
        comet_cloud_ids.push(id);
    }
    let mut debris_disk_ids = Vec::new();
    let mut debris_disks = Vec::new();
    if primary_stage != StellarStage::Disrupted
        && population.epoch_age.seconds() <= 200_000_000 * SECONDS_PER_YEAR
        && system_rng.next_f64() < 0.80
    {
        let id = DebrisDiskId::new(derive_object_seed(seed, system_key, "debris-disk-id-v1"));
        debris_disks.push(super::planetary::generate_debris_disk(
            id,
            system_id,
            &mut system_rng,
        )?);
        debris_disk_ids.push(id);
    }

    let mut protoplanetary_disk_ids = Vec::new();
    let mut protoplanetary_disks = Vec::new();
    if primary_stage != StellarStage::Disrupted {
        let id = ProtoplanetaryDiskId::new(derive_object_seed(
            seed,
            system_key,
            "protoplanetary-disk-id-v1",
        ));
        protoplanetary_disks.push(super::planetary::generate_protoplanetary_disk(
            id,
            system_id,
            current_primary_mass,
            population.epoch_age,
            &mut system_rng,
        )?);
        protoplanetary_disk_ids.push(id);
    }

    Ok(GeneratedSystem {
        system: StarSystem {
            id: system_id,
            galaxy_id,
            cluster_id: None,
            position,
            velocity,
            star_id,
            companion_star_ids,
            planet_ids,
            moon_ids,
            asteroid_ids,
            comet_ids,
            asteroid_belt_ids,
            comet_cloud_ids,
            debris_disk_ids,
            protoplanetary_disk_ids,
            component: population.component,
        },
        stars,
        planets,
        moons,
        asteroids,
        comets,
        asteroid_belts,
        comet_clouds,
        debris_disks,
        protoplanetary_disks,
    })
}

struct SystemPopulation {
    component: GalacticComponent,
    position: Position,
    velocity: Velocity,
    epoch_age: SimDuration,
    metallicity: f64,
}

fn generate_system_population(
    rng: &mut StableRng,
    galaxy: GalaxyCatalogEntry,
) -> Result<SystemPopulation, ObjectError> {
    let component = match (galaxy.morphology, rng.next_f64()) {
        (GalaxyMorphology::Elliptical, value) if value < 0.85 => GalacticComponent::Bulge,
        (GalaxyMorphology::Dwarf | GalaxyMorphology::Irregular, value) if value < 0.70 => {
            GalacticComponent::Disk
        }
        (_, value) if value < 0.80 => GalacticComponent::Disk,
        (_, value) if value < 0.95 => GalacticComponent::Bulge,
        _ => GalacticComponent::Halo,
    };
    let galaxy_radius = galaxy.radius_parsecs * METERS_PER_PARSEC;
    let (position, velocity, age_years, base_metallicity) = match component {
        GalacticComponent::Disk => {
            let radius = galaxy_radius * rng.next_f64().sqrt();
            let angle = std::f64::consts::TAU * rng.next_f64();
            let (sin_angle, cos_angle) = angle.sin_cos();
            let scale_height = galaxy_radius * 0.015;
            let vertical = super::random_signed(rng)
                * scale_height
                * (-(1.0 - rng.next_f64()).max(1e-12).ln());
            let circular_speed = 220_000.0;
            (
                Position::from_meters(radius * cos_angle, radius * sin_angle, vertical)?,
                Velocity::from_meters_per_second(
                    -circular_speed * sin_angle + super::random_signed(rng) * 20_000.0,
                    circular_speed * cos_angle + super::random_signed(rng) * 20_000.0,
                    super::random_signed(rng) * 10_000.0,
                )?,
                super::random_range(rng, 0.0, 10_000_000_000.0),
                0.024 * (-1.5 * radius / galaxy_radius).exp(),
            )
        }
        GalacticComponent::Bulge => {
            let radius = galaxy_radius * 0.15 * rng.next_f64().cbrt();
            let direction = super::random_unit_vector(rng);
            (
                Position::from_meters(
                    radius * direction[0],
                    radius * direction[1],
                    radius * direction[2],
                )?,
                Velocity::from_meters_per_second(
                    super::random_signed(rng) * 140_000.0,
                    super::random_signed(rng) * 140_000.0,
                    super::random_signed(rng) * 140_000.0,
                )?,
                super::random_range(rng, 8_000_000_000.0, 12_000_000_000.0),
                0.022,
            )
        }
        GalacticComponent::Halo => {
            let radius = galaxy_radius * rng.next_f64().cbrt();
            let direction = super::random_unit_vector(rng);
            (
                Position::from_meters(
                    radius * direction[0],
                    radius * direction[1],
                    radius * direction[2],
                )?,
                Velocity::from_meters_per_second(
                    super::random_signed(rng) * super::MAX_SYSTEM_SPEED_METERS_PER_SECOND,
                    super::random_signed(rng) * super::MAX_SYSTEM_SPEED_METERS_PER_SECOND,
                    super::random_signed(rng) * super::MAX_SYSTEM_SPEED_METERS_PER_SECOND,
                )?,
                super::random_range(rng, 10_000_000_000.0, 13_000_000_000.0),
                0.0015,
            )
        }
    };
    let (age_years, base_metallicity) = match galaxy.morphology {
        GalaxyMorphology::Spiral => (age_years, base_metallicity),
        GalaxyMorphology::Elliptical => (super::random_range(rng, 8.0e9, 13.0e9), 0.022),
        GalaxyMorphology::Irregular => (super::random_range(rng, 0.0, 8.0e9), 0.008),
        GalaxyMorphology::Dwarf => (super::random_range(rng, 0.0, 12.0e9), 0.003),
    };
    let age_metallicity_factor = (1.0 - 0.35 * age_years / 13_000_000_000.0).max(0.2);
    let metallicity = (base_metallicity * age_metallicity_factor).clamp(0.0001, 0.04);
    let epoch_age = SimDuration::from_years(age_years.round() as u128)
        .map_err(|error| ObjectError(error.to_string()))?;
    Ok(SystemPopulation {
        component,
        position,
        velocity,
        epoch_age,
        metallicity,
    })
}

fn generate_star(
    seed: u64,
    id: StarId,
    system_id: StarSystemId,
    epoch_age: SimDuration,
    metallicity: f64,
    config: GenerationConfig,
) -> Result<Star, ObjectError> {
    let mut rng = StableRng::new(derive_object_seed(seed, id.value(), "star-properties-v1"));
    let mass_solar = super::sample_stellar_mass(
        &mut rng,
        config.star_mass_solar_min,
        config.star_mass_solar_max,
    );
    generate_star_with_mass(id, system_id, epoch_age, metallicity, mass_solar)
}

fn generate_star_with_mass(
    id: StarId,
    system_id: StarSystemId,
    epoch_age: SimDuration,
    metallicity: f64,
    mass_solar: f64,
) -> Result<Star, ObjectError> {
    let initial_mass = Mass::from_solar_masses(mass_solar)?;
    let stage = super::stellar::stellar_stage_at(initial_mass, metallicity, epoch_age)?;
    let (current_mass, luminosity, temperature) =
        super::stellar::stellar_properties(initial_mass, metallicity, stage)?;
    Ok(Star {
        id,
        system_id,
        initial_mass,
        current_mass,
        epoch_age,
        formed_at: SimDuration::ZERO,
        age: epoch_age,
        metallicity,
        birth_composition: GasComposition::from_metallicity(metallicity)
            .map_err(|error| ObjectError(error.to_string()))?,
        luminosity,
        temperature,
        stage,
        position: Position::default(),
        velocity: Velocity::default(),
        orbital_elements: None,
    })
}

pub(super) fn generate_star_with_birth(
    id: StarId,
    system_id: StarSystemId,
    formed_at: SimDuration,
    birth_composition: GasComposition,
    mass_solar: f64,
) -> Result<Star, ObjectError> {
    let metallicity = birth_composition.metallicity();
    let initial_mass = Mass::from_solar_masses(mass_solar)?;
    let stage = StellarStage::Protostar;
    let (current_mass, luminosity, temperature) =
        super::stellar::stellar_properties(initial_mass, metallicity, stage)?;
    Ok(Star {
        id,
        system_id,
        initial_mass,
        current_mass,
        epoch_age: SimDuration::ZERO,
        formed_at,
        age: SimDuration::ZERO,
        metallicity,
        birth_composition,
        luminosity,
        temperature,
        stage,
        position: Position::default(),
        velocity: Velocity::default(),
        orbital_elements: None,
    })
}
