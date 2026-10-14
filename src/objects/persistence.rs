use std::fmt::Write as _;
use std::str::FromStr;

use crate::chemistry::GasComposition;
use crate::id::{
    AsteroidBeltId, AsteroidId, CometCloudId, CometId, DebrisDiskId, GalaxyId, MoonId, PlanetId,
    ProtoplanetaryDiskId, StarId, StarSystemId,
};
use crate::orbital::OrbitalElements;
use crate::time::SimDuration;
use crate::units::{Angle, Length, Luminosity, Mass, Position, Temperature, Velocity};

use super::{
    Asteroid, AsteroidBelt, Comet, CometCloud, DebrisDisk, GalacticComponent, Galaxy, Moon,
    MoonFormationHistory, MoonOrigin, ObjectError, ObjectRegistry, Planet, PlanetAtmosphere,
    PlanetStatus, PlanetaryComposition, ProtoplanetaryDisk, Star, StarSystem, StellarStage,
};

impl ObjectRegistry {
    pub(crate) fn encode(&self, output: &mut String) {
        output.push_str("[objects]\n");
        writeln!(
            output,
            "M|{}|{}|{}|{}|{}|{}|{}|{}",
            float_bits(self.gas_reservoir_solar_masses),
            self.formation_sequence,
            float_bits(self.gas_composition.hydrogen_fraction()),
            float_bits(self.gas_composition.helium_fraction()),
            float_bits(self.gas_composition.carbon_fraction()),
            float_bits(self.gas_composition.oxygen_fraction()),
            float_bits(self.gas_composition.iron_fraction()),
            float_bits(self.gas_composition.other_metals_fraction()),
        )
        .expect("writing to a String cannot fail");
        for galaxy in self.galaxies.values() {
            let [x, y, z] = galaxy.position.meters();
            writeln!(
                output,
                "G|{}|{}|{}|{}|{}",
                galaxy.id.value(),
                float_bits(x),
                float_bits(y),
                float_bits(z),
                encode_ids(&galaxy.star_system_ids)
            )
            .expect("writing to a String cannot fail");
        }
        for system in self.star_systems.values() {
            let [px, py, pz] = system.position.meters();
            let [vx, vy, vz] = system.velocity.meters_per_second();
            writeln!(
                output,
                "Y|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                system.id.value(),
                system.galaxy_id.value(),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                system.star_id.value(),
                encode_ids(&system.planet_ids),
                system.component.code(),
                encode_ids(&system.companion_star_ids),
                encode_ids(&system.moon_ids),
                encode_ids(&system.asteroid_ids),
                encode_ids(&system.comet_ids),
                encode_ids(&system.asteroid_belt_ids),
                encode_ids(&system.comet_cloud_ids),
                encode_ids(&system.debris_disk_ids),
                encode_ids(&system.protoplanetary_disk_ids),
            )
            .expect("writing to a String cannot fail");
        }
        for star in self.stars.values() {
            let [px, py, pz] = star.position.meters();
            let [vx, vy, vz] = star.velocity.meters_per_second();
            writeln!(
                output,
                "S|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                star.id.value(),
                star.system_id.value(),
                float_bits(star.initial_mass.kilograms()),
                float_bits(star.current_mass.kilograms()),
                star.epoch_age.seconds(),
                star.age.seconds(),
                float_bits(star.metallicity),
                float_bits(star.luminosity.watts()),
                float_bits(star.temperature.kelvin()),
                star.stage.code(),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                encode_optional_orbit(star.orbital_elements),
                star.formed_at.seconds(),
                float_bits(star.birth_composition.hydrogen_fraction()),
                float_bits(star.birth_composition.helium_fraction()),
                float_bits(star.birth_composition.carbon_fraction()),
                float_bits(star.birth_composition.oxygen_fraction()),
                float_bits(star.birth_composition.iron_fraction()),
                float_bits(star.birth_composition.other_metals_fraction()),
            )
            .expect("writing to a String cannot fail");
        }
        for planet in self.planets.values() {
            let [px, py, pz] = planet.position.meters();
            let [vx, vy, vz] = planet.velocity.meters_per_second();
            writeln!(
                output,
                "P|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                planet.id.value(),
                planet.system_id.value(),
                planet.star_id.value(),
                float_bits(planet.mass.kilograms()),
                float_bits(planet.orbital_radius.meters()),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                float_bits(planet.orbital_elements.semi_major_axis().meters()),
                float_bits(planet.orbital_elements.eccentricity()),
                float_bits(planet.orbital_elements.inclination().radians()),
                float_bits(
                    planet
                        .orbital_elements
                        .longitude_of_ascending_node()
                        .radians()
                ),
                float_bits(planet.orbital_elements.argument_of_periapsis().radians()),
                float_bits(planet.orbital_elements.mean_anomaly_at_epoch().radians()),
                planet.orbital_elements.epoch().seconds(),
                planet.status.code(),
                float_bits(planet.radius.meters()),
                float_bits(planet.equilibrium_temperature.kelvin()),
                planet.atmosphere.code(),
                float_bits(planet.atmosphere_loss_fraction),
                u8::from(planet.in_habitable_zone),
                float_bits(planet.bulk_composition.rocky_fraction()),
                float_bits(planet.bulk_composition.icy_fraction()),
                float_bits(planet.bulk_composition.iron_fraction()),
                float_bits(planet.bulk_composition.hydrogen_helium_fraction()),
                planet.formed_at.seconds(),
                float_bits(planet.initial_orbital_radius.meters()),
                float_bits(planet.migration_rate_au_per_year),
                planet.scattering_count,
            )
            .expect("writing to a String cannot fail");
        }
        for moon in self.moons.values() {
            let [px, py, pz] = moon.position.meters();
            let [vx, vy, vz] = moon.velocity.meters_per_second();
            writeln!(
                output,
                "N|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                moon.id.value(),
                moon.planet_id.value(),
                moon.system_id.value(),
                float_bits(moon.mass.kilograms()),
                float_bits(moon.radius.meters()),
                float_bits(moon.orbital_radius.meters()),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                float_bits(moon.orbital_elements.semi_major_axis().meters()),
                float_bits(moon.orbital_elements.eccentricity()),
                float_bits(moon.orbital_elements.inclination().radians()),
                float_bits(
                    moon.orbital_elements
                        .longitude_of_ascending_node()
                        .radians()
                ),
                float_bits(moon.orbital_elements.argument_of_periapsis().radians()),
                float_bits(moon.orbital_elements.mean_anomaly_at_epoch().radians()),
                moon.orbital_elements.epoch().seconds(),
                moon.atmosphere.code(),
                float_bits(moon.equilibrium_temperature.kelvin()),
                moon.formation_history.origin().code(),
                moon.formation_history.formed_at().seconds(),
                moon.formation_history
                    .captured_at()
                    .map_or_else(|| "-".to_owned(), |time| time.seconds().to_string()),
                float_bits(moon.tidal_migration_rate_m_per_year),
                float_bits(moon.initial_orbital_radius.meters()),
            )
            .expect("writing to a String cannot fail");
        }
        for belt in self.asteroid_belts.values() {
            writeln!(
                output,
                "A|{}|{}|{}|{}|{}|{}",
                belt.id.value(),
                belt.system_id.value(),
                float_bits(belt.inner_radius.meters()),
                float_bits(belt.outer_radius.meters()),
                belt.estimated_body_count,
                float_bits(belt.total_mass.kilograms()),
            )
            .expect("writing to a String cannot fail");
        }
        for cloud in self.comet_clouds.values() {
            writeln!(
                output,
                "C|{}|{}|{}|{}|{}|{}",
                cloud.id.value(),
                cloud.system_id.value(),
                float_bits(cloud.inner_radius.meters()),
                float_bits(cloud.outer_radius.meters()),
                cloud.estimated_body_count,
                float_bits(cloud.total_mass.kilograms()),
            )
            .expect("writing to a String cannot fail");
        }
        for asteroid in self.asteroids.values() {
            let [px, py, pz] = asteroid.position.meters();
            let [vx, vy, vz] = asteroid.velocity.meters_per_second();
            writeln!(
                output,
                "X|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                asteroid.id.value(),
                asteroid.belt_id.value(),
                asteroid.system_id.value(),
                float_bits(asteroid.mass.kilograms()),
                float_bits(asteroid.radius.meters()),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                float_bits(asteroid.orbital_elements.semi_major_axis().meters()),
                float_bits(asteroid.orbital_elements.eccentricity()),
                float_bits(asteroid.orbital_elements.inclination().radians()),
                float_bits(
                    asteroid
                        .orbital_elements
                        .longitude_of_ascending_node()
                        .radians()
                ),
                float_bits(asteroid.orbital_elements.argument_of_periapsis().radians()),
                float_bits(asteroid.orbital_elements.mean_anomaly_at_epoch().radians()),
                asteroid.orbital_elements.epoch().seconds(),
            )
            .expect("writing to a String cannot fail");
        }
        for comet in self.comets.values() {
            let [px, py, pz] = comet.position.meters();
            let [vx, vy, vz] = comet.velocity.meters_per_second();
            writeln!(
                output,
                "K|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                comet.id.value(),
                comet.cloud_id.value(),
                comet.system_id.value(),
                float_bits(comet.mass.kilograms()),
                float_bits(comet.radius.meters()),
                float_bits(comet.volatile_fraction),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                float_bits(comet.orbital_elements.semi_major_axis().meters()),
                float_bits(comet.orbital_elements.eccentricity()),
                float_bits(comet.orbital_elements.inclination().radians()),
                float_bits(
                    comet
                        .orbital_elements
                        .longitude_of_ascending_node()
                        .radians()
                ),
                float_bits(comet.orbital_elements.argument_of_periapsis().radians()),
                float_bits(comet.orbital_elements.mean_anomaly_at_epoch().radians()),
                comet.orbital_elements.epoch().seconds(),
            )
            .expect("writing to a String cannot fail");
        }
        for disk in self.debris_disks.values() {
            writeln!(
                output,
                "D|{}|{}|{}|{}|{}|{}",
                disk.id.value(),
                disk.system_id.value(),
                float_bits(disk.inner_radius.meters()),
                float_bits(disk.outer_radius.meters()),
                disk.estimated_body_count,
                float_bits(disk.total_mass.kilograms()),
            )
            .expect("writing to a String cannot fail");
        }
        for disk in self.protoplanetary_disks.values() {
            writeln!(
                output,
                "Q|{}|{}|{}|{}|{}|{}|{}|{}",
                disk.id.value(),
                disk.system_id.value(),
                float_bits(disk.inner_radius.meters()),
                float_bits(disk.outer_radius.meters()),
                float_bits(disk.initial_mass.kilograms()),
                float_bits(disk.remaining_mass.kilograms()),
                disk.formed_at.seconds(),
                disk.dissipates_at.seconds()
            )
            .expect("writing to a String cannot fail");
        }
    }

    pub(crate) fn decode<'a>(lines: impl Iterator<Item = &'a str>) -> Result<Self, ObjectError> {
        let mut registry = Self::default();
        for line in lines {
            if line.is_empty() {
                continue;
            }
            let fields: Vec<_> = line.split('|').collect();
            match fields.first().copied() {
                Some("M") => {
                    if !matches!(fields.len(), 3 | 9) {
                        return Err(ObjectError(format!(
                            "matter-reservoir record has {} fields; expected 3 or 9",
                            fields.len()
                        )));
                    }
                    registry.gas_reservoir_solar_masses =
                        parse_float_bits(fields[1], "gas reservoir mass")?;
                    registry.formation_sequence = parse_number(fields[2], "formation sequence")?;
                    if fields.len() == 9 {
                        registry.gas_composition = GasComposition::new(
                            parse_float_bits(fields[3], "hydrogen fraction")?,
                            parse_float_bits(fields[4], "helium fraction")?,
                            parse_float_bits(fields[5], "carbon fraction")?,
                            parse_float_bits(fields[6], "oxygen fraction")?,
                            parse_float_bits(fields[7], "iron fraction")?,
                            parse_float_bits(fields[8], "other-metals fraction")?,
                        )
                        .map_err(|error| ObjectError(error.to_string()))?;
                    }
                    if !registry.gas_reservoir_solar_masses.is_finite()
                        || registry.gas_reservoir_solar_masses < 0.0
                    {
                        return Err(ObjectError(
                            "gas reservoir mass must be finite and nonnegative".to_owned(),
                        ));
                    }
                }
                Some("G") => {
                    require_field_count(&fields, 6, "galaxy")?;
                    let galaxy = Galaxy {
                        id: GalaxyId::new(parse_number(fields[1], "galaxy ID")?),
                        position: parse_position(&fields[2..5])?,
                        star_system_ids: parse_ids(fields[5], StarSystemId::new)?,
                    };
                    super::insert_unique(&mut registry.galaxies, galaxy.id, galaxy, "galaxy")?;
                }
                Some("Y") => {
                    if !matches!(fields.len(), 11..=13 | 17 | 20) {
                        return Err(ObjectError(format!(
                            "star-system record has {} fields; expected 11, 12, 13, or 17",
                            fields.len()
                        )));
                    }
                    let system = StarSystem {
                        id: StarSystemId::new(parse_number(fields[1], "star-system ID")?),
                        galaxy_id: GalaxyId::new(parse_number(fields[2], "galaxy ID")?),
                        position: parse_position(&fields[3..6])?,
                        velocity: parse_velocity(&fields[6..9])?,
                        star_id: StarId::new(parse_number(fields[9], "star ID")?),
                        planet_ids: parse_ids(fields[10], PlanetId::new)?,
                        component: if fields.len() >= 12 {
                            GalacticComponent::parse_code(fields[11])?
                        } else {
                            GalacticComponent::Disk
                        },
                        companion_star_ids: if matches!(fields.len(), 13 | 17 | 20) {
                            parse_ids(fields[12], StarId::new)?
                        } else {
                            Vec::new()
                        },
                        moon_ids: if matches!(fields.len(), 17 | 20) {
                            parse_ids(fields[13], MoonId::new)?
                        } else {
                            Vec::new()
                        },
                        asteroid_ids: if fields.len() == 20 {
                            parse_ids(fields[14], AsteroidId::new)?
                        } else {
                            Vec::new()
                        },
                        comet_ids: if fields.len() == 20 {
                            parse_ids(fields[15], CometId::new)?
                        } else {
                            Vec::new()
                        },
                        asteroid_belt_ids: if fields.len() == 20 {
                            parse_ids(fields[16], AsteroidBeltId::new)?
                        } else if fields.len() == 17 {
                            parse_ids(fields[14], AsteroidBeltId::new)?
                        } else {
                            Vec::new()
                        },
                        comet_cloud_ids: if fields.len() == 20 {
                            parse_ids(fields[17], CometCloudId::new)?
                        } else if fields.len() == 17 {
                            parse_ids(fields[15], CometCloudId::new)?
                        } else {
                            Vec::new()
                        },
                        debris_disk_ids: if fields.len() == 20 {
                            parse_ids(fields[18], DebrisDiskId::new)?
                        } else if fields.len() == 17 {
                            parse_ids(fields[16], DebrisDiskId::new)?
                        } else {
                            Vec::new()
                        },
                        protoplanetary_disk_ids: if fields.len() == 20 {
                            parse_ids(fields[19], ProtoplanetaryDiskId::new)?
                        } else {
                            Vec::new()
                        },
                    };
                    super::insert_unique(
                        &mut registry.star_systems,
                        system.id,
                        system,
                        "star system",
                    )?;
                }
                Some("S") => {
                    if !matches!(fields.len(), 8 | 9 | 11 | 18 | 19 | 25) {
                        return Err(ObjectError(format!(
                            "star record has {} fields; expected 8, 9, 11, 18, 19, or 25",
                            fields.len()
                        )));
                    }
                    let is_current = matches!(fields.len(), 11 | 18 | 19 | 25);
                    let metallicity_index = if is_current { 7 } else { 5 };
                    let metallicity = parse_float_bits(fields[metallicity_index], "metallicity")?;
                    if !metallicity.is_finite() || metallicity < 0.0 {
                        return Err(ObjectError(
                            "star metallicity must be finite and nonnegative".to_owned(),
                        ));
                    }
                    let initial_mass =
                        Mass::from_kilograms(parse_float_bits(fields[3], "initial star mass")?)?;
                    let current_mass = if is_current {
                        Mass::from_kilograms(parse_float_bits(fields[4], "current star mass")?)?
                    } else {
                        initial_mass
                    };
                    let epoch_age = if is_current {
                        SimDuration::from_seconds(parse_number(fields[5], "epoch star age")?)
                    } else {
                        SimDuration::ZERO
                    };
                    let age_index = if is_current { 6 } else { 4 };
                    let age =
                        SimDuration::from_seconds(parse_number(fields[age_index], "star age")?);
                    let stage = if is_current {
                        StellarStage::parse_code(fields[10])?
                    } else if fields.len() == 9 {
                        StellarStage::parse_code(fields[8])?
                    } else {
                        super::stellar::stellar_stage_at(initial_mass, metallicity, age)?
                    };
                    let luminosity_index = if is_current { 8 } else { 6 };
                    let temperature_index = if is_current { 9 } else { 7 };
                    let star = Star {
                        id: StarId::new(parse_number(fields[1], "star ID")?),
                        system_id: StarSystemId::new(parse_number(fields[2], "system ID")?),
                        initial_mass,
                        current_mass,
                        epoch_age,
                        formed_at: if fields.len() >= 19 {
                            SimDuration::from_seconds(parse_number(fields[18], "formation time")?)
                        } else {
                            SimDuration::ZERO
                        },
                        age,
                        metallicity,
                        birth_composition: if fields.len() == 25 {
                            GasComposition::new(
                                parse_float_bits(fields[19], "birth hydrogen")?,
                                parse_float_bits(fields[20], "birth helium")?,
                                parse_float_bits(fields[21], "birth carbon")?,
                                parse_float_bits(fields[22], "birth oxygen")?,
                                parse_float_bits(fields[23], "birth iron")?,
                                parse_float_bits(fields[24], "birth other metals")?,
                            )
                        } else {
                            GasComposition::from_metallicity(metallicity)
                        }
                        .map_err(|error| ObjectError(error.to_string()))?,
                        luminosity: Luminosity::from_watts(parse_float_bits(
                            fields[luminosity_index],
                            "star luminosity",
                        )?)?,
                        temperature: Temperature::from_kelvin(parse_float_bits(
                            fields[temperature_index],
                            "star temperature",
                        )?)?,
                        stage,
                        position: if fields.len() >= 18 {
                            parse_position(&fields[11..14])?
                        } else {
                            Position::default()
                        },
                        velocity: if fields.len() >= 18 {
                            parse_velocity(&fields[14..17])?
                        } else {
                            Velocity::default()
                        },
                        orbital_elements: if fields.len() >= 18 {
                            parse_optional_orbit(fields[17])?
                        } else {
                            None
                        },
                    };
                    super::insert_unique(&mut registry.stars, star.id, star, "star")?;
                }
                Some("P") => {
                    if !matches!(fields.len(), 12 | 19 | 20 | 23 | 30 | 33) {
                        return Err(ObjectError(format!(
                            "planet record has {} fields; expected 12, 19, 20, or 23",
                            fields.len()
                        )));
                    }
                    let orbital_radius =
                        Length::from_meters(parse_float_bits(fields[5], "orbital radius")?)?;
                    let orbital_elements = if fields.len() >= 19 {
                        OrbitalElements::new(
                            Length::from_meters(parse_float_bits(fields[12], "semi-major axis")?)?,
                            parse_float_bits(fields[13], "eccentricity")?,
                            Angle::from_radians(parse_float_bits(fields[14], "inclination")?)?,
                            Angle::from_radians(parse_float_bits(fields[15], "ascending node")?)?,
                            Angle::from_radians(parse_float_bits(
                                fields[16],
                                "periapsis argument",
                            )?)?,
                            Angle::from_radians(parse_float_bits(fields[17], "mean anomaly")?)?,
                            SimDuration::from_seconds(parse_number(fields[18], "orbit epoch")?),
                        )?
                    } else {
                        OrbitalElements::new(
                            orbital_radius,
                            0.0,
                            Angle::default(),
                            Angle::default(),
                            Angle::default(),
                            Angle::default(),
                            SimDuration::ZERO,
                        )?
                    };
                    let mass = Mass::from_kilograms(parse_float_bits(fields[4], "planet mass")?)?;
                    let has_planetary_state = fields.len() == 30 || fields.len() == 33;
                    let has_migration_state = fields.len() == 33;
                    let atmosphere = if has_planetary_state || fields.len() == 23 {
                        PlanetAtmosphere::parse_code(fields[22])?
                    } else {
                        PlanetAtmosphere::Thin
                    };
                    let planet = Planet {
                        id: PlanetId::new(parse_number(fields[1], "planet ID")?),
                        system_id: StarSystemId::new(parse_number(fields[2], "system ID")?),
                        star_id: StarId::new(parse_number(fields[3], "star ID")?),
                        mass,
                        radius: if has_planetary_state || fields.len() == 23 {
                            Length::from_meters(parse_float_bits(fields[20], "planet radius")?)?
                        } else {
                            Length::from_meters(super::planetary::planet_radius_meters(
                                mass, atmosphere,
                            ))?
                        },
                        atmosphere,
                        atmosphere_loss_fraction: if has_planetary_state {
                            parse_float_bits(fields[23], "atmosphere loss fraction")?
                        } else {
                            0.0
                        },
                        equilibrium_temperature: if has_planetary_state || fields.len() == 23 {
                            Temperature::from_kelvin(parse_float_bits(
                                fields[21],
                                "planet equilibrium temperature",
                            )?)?
                        } else {
                            Temperature::from_kelvin(288.0)?
                        },
                        in_habitable_zone: if has_planetary_state {
                            parse_number::<u8>(fields[24], "habitable-zone flag")? != 0
                        } else {
                            false
                        },
                        bulk_composition: if has_planetary_state {
                            PlanetaryComposition::new(
                                parse_float_bits(fields[25], "rocky composition")?,
                                parse_float_bits(fields[26], "icy composition")?,
                                parse_float_bits(fields[27], "iron composition")?,
                                parse_float_bits(fields[28], "hydrogen-helium composition")?,
                            )
                        } else {
                            PlanetaryComposition::new(0.8, 0.1, 0.1, 0.0)
                        },
                        formed_at: if has_planetary_state {
                            SimDuration::from_seconds(parse_number(
                                fields[29],
                                "planet formation time",
                            )?)
                        } else {
                            SimDuration::ZERO
                        },
                        initial_orbital_radius: if has_migration_state {
                            Length::from_meters(parse_float_bits(
                                fields[30],
                                "initial planet orbital radius",
                            )?)?
                        } else {
                            orbital_radius
                        },
                        migration_rate_au_per_year: if has_migration_state {
                            parse_float_bits(fields[31], "planet migration rate")?
                        } else {
                            0.0
                        },
                        scattering_count: if has_migration_state {
                            parse_number(fields[32], "planet scattering count")?
                        } else {
                            0
                        },
                        orbital_radius,
                        position: parse_position(&fields[6..9])?,
                        velocity: parse_velocity(&fields[9..12])?,
                        orbital_elements,
                        status: if matches!(fields.len(), 20 | 23 | 30 | 33) {
                            PlanetStatus::parse_code(fields[19])?
                        } else {
                            PlanetStatus::Stable
                        },
                    };
                    super::insert_unique(&mut registry.planets, planet.id, planet, "planet")?;
                }
                Some("N") => {
                    if !matches!(fields.len(), 19 | 20 | 26 | 27) {
                        return Err(ObjectError(format!(
                            "moon record has {} fields; expected 19, 20, 26, or 27",
                            fields.len()
                        )));
                    }
                    let has_radius = matches!(fields.len(), 20 | 26 | 27);
                    let has_moon_state = fields.len() == 26 || fields.len() == 27;
                    let has_tidal_state = fields.len() == 27;
                    let mass = Mass::from_kilograms(parse_float_bits(fields[4], "moon mass")?)?;
                    let radius = if has_radius {
                        Length::from_meters(parse_float_bits(fields[5], "moon radius")?)?
                    } else {
                        Length::from_meters(super::planetary::planet_radius_meters(
                            mass,
                            PlanetAtmosphere::None,
                        ))?
                    };
                    let orbital_radius_index = if has_radius { 6 } else { 5 };
                    let position_start = if has_radius { 7 } else { 6 };
                    let orbit_start = if has_radius { 13 } else { 12 };
                    let orbital_radius = Length::from_meters(parse_float_bits(
                        fields[orbital_radius_index],
                        "moon orbital radius",
                    )?)?;
                    let orbital_elements = OrbitalElements::new(
                        Length::from_meters(parse_float_bits(
                            fields[orbit_start],
                            "moon semi-major axis",
                        )?)?,
                        parse_float_bits(fields[orbit_start + 1], "moon eccentricity")?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 2],
                            "moon inclination",
                        )?)?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 3],
                            "moon ascending node",
                        )?)?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 4],
                            "moon periapsis argument",
                        )?)?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 5],
                            "moon mean anomaly",
                        )?)?,
                        SimDuration::from_seconds(parse_number(
                            fields[orbit_start + 6],
                            "moon orbit epoch",
                        )?),
                    )?;
                    let moon = Moon {
                        id: MoonId::new(parse_number(fields[1], "moon ID")?),
                        planet_id: PlanetId::new(parse_number(fields[2], "moon planet ID")?),
                        system_id: StarSystemId::new(parse_number(fields[3], "moon system ID")?),
                        mass,
                        radius,
                        orbital_radius,
                        position: parse_position(&fields[position_start..position_start + 3])?,
                        velocity: parse_velocity(&fields[position_start + 3..position_start + 6])?,
                        orbital_elements,
                        atmosphere: if has_moon_state {
                            PlanetAtmosphere::parse_code(fields[20])?
                        } else {
                            PlanetAtmosphere::None
                        },
                        equilibrium_temperature: if has_moon_state {
                            Temperature::from_kelvin(parse_float_bits(
                                fields[21],
                                "moon equilibrium temperature",
                            )?)?
                        } else {
                            Temperature::from_kelvin(288.0)?
                        },
                        formation_history: if has_moon_state {
                            MoonFormationHistory::new(
                                MoonOrigin::parse_code(fields[22])?,
                                SimDuration::from_seconds(parse_number(
                                    fields[23],
                                    "moon formation time",
                                )?),
                                if fields[24] == "-" {
                                    None
                                } else {
                                    Some(SimDuration::from_seconds(parse_number(
                                        fields[24],
                                        "moon capture time",
                                    )?))
                                },
                            )
                        } else {
                            MoonFormationHistory::new(MoonOrigin::InSitu, SimDuration::ZERO, None)
                        },
                        tidal_migration_rate_m_per_year: if has_moon_state {
                            parse_float_bits(fields[25], "moon tidal migration rate")?
                        } else {
                            0.0
                        },
                        initial_orbital_radius: if has_tidal_state {
                            Length::from_meters(parse_float_bits(
                                fields[26],
                                "initial moon orbital radius",
                            )?)?
                        } else {
                            orbital_radius
                        },
                    };
                    super::insert_unique(&mut registry.moons, moon.id, moon, "moon")?;
                }
                Some("X") => {
                    require_field_count(&fields, 19, "asteroid")?;
                    let asteroid = Asteroid {
                        id: AsteroidId::new(parse_number(fields[1], "asteroid ID")?),
                        belt_id: AsteroidBeltId::new(parse_number(fields[2], "asteroid belt ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[3],
                            "asteroid system ID",
                        )?),
                        mass: Mass::from_kilograms(parse_float_bits(fields[4], "asteroid mass")?)?,
                        radius: Length::from_meters(parse_float_bits(
                            fields[5],
                            "asteroid radius",
                        )?)?,
                        position: parse_position(&fields[6..9])?,
                        velocity: parse_velocity(&fields[9..12])?,
                        orbital_elements: parse_orbit(&fields, 12, "asteroid")?,
                    };
                    super::insert_unique(
                        &mut registry.asteroids,
                        asteroid.id,
                        asteroid,
                        "asteroid",
                    )?;
                }
                Some("K") => {
                    require_field_count(&fields, 20, "comet")?;
                    let comet = Comet {
                        id: CometId::new(parse_number(fields[1], "comet ID")?),
                        cloud_id: CometCloudId::new(parse_number(fields[2], "comet cloud ID")?),
                        system_id: StarSystemId::new(parse_number(fields[3], "comet system ID")?),
                        mass: Mass::from_kilograms(parse_float_bits(fields[4], "comet mass")?)?,
                        radius: Length::from_meters(parse_float_bits(fields[5], "comet radius")?)?,
                        volatile_fraction: parse_float_bits(fields[6], "comet volatile fraction")?,
                        position: parse_position(&fields[7..10])?,
                        velocity: parse_velocity(&fields[10..13])?,
                        orbital_elements: parse_orbit(&fields, 13, "comet")?,
                    };
                    super::insert_unique(&mut registry.comets, comet.id, comet, "comet")?;
                }
                Some("A") => {
                    require_field_count(&fields, 7, "asteroid belt")?;
                    let belt = AsteroidBelt {
                        id: AsteroidBeltId::new(parse_number(fields[1], "asteroid-belt ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "asteroid-belt system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "asteroid-belt inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "asteroid-belt outer radius",
                        )?)?,
                        estimated_body_count: parse_number(fields[5], "asteroid-belt body count")?,
                        total_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "asteroid-belt mass",
                        )?)?,
                    };
                    super::insert_unique(
                        &mut registry.asteroid_belts,
                        belt.id,
                        belt,
                        "asteroid belt",
                    )?;
                }
                Some("C") => {
                    require_field_count(&fields, 7, "comet cloud")?;
                    let cloud = CometCloud {
                        id: CometCloudId::new(parse_number(fields[1], "comet-cloud ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "comet-cloud system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "comet-cloud inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "comet-cloud outer radius",
                        )?)?,
                        estimated_body_count: parse_number(fields[5], "comet-cloud body count")?,
                        total_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "comet-cloud mass",
                        )?)?,
                    };
                    super::insert_unique(
                        &mut registry.comet_clouds,
                        cloud.id,
                        cloud,
                        "comet cloud",
                    )?;
                }
                Some("D") => {
                    require_field_count(&fields, 7, "debris disk")?;
                    let disk = DebrisDisk {
                        id: DebrisDiskId::new(parse_number(fields[1], "debris-disk ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "debris-disk system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "debris-disk inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "debris-disk outer radius",
                        )?)?,
                        estimated_body_count: parse_number(fields[5], "debris-disk body count")?,
                        total_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "debris-disk mass",
                        )?)?,
                    };
                    super::insert_unique(&mut registry.debris_disks, disk.id, disk, "debris disk")?;
                }
                Some("Q") => {
                    require_field_count(&fields, 9, "protoplanetary disk")?;
                    let disk = ProtoplanetaryDisk {
                        id: ProtoplanetaryDiskId::new(parse_number(
                            fields[1],
                            "protoplanetary-disk ID",
                        )?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "protoplanetary-disk system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "protoplanetary-disk inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "protoplanetary-disk outer radius",
                        )?)?,
                        initial_mass: Mass::from_kilograms(parse_float_bits(
                            fields[5],
                            "protoplanetary-disk initial mass",
                        )?)?,
                        remaining_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "protoplanetary-disk remaining mass",
                        )?)?,
                        formed_at: SimDuration::from_seconds(parse_number(
                            fields[7],
                            "protoplanetary-disk formation time",
                        )?),
                        dissipates_at: SimDuration::from_seconds(parse_number(
                            fields[8],
                            "protoplanetary-disk dissipation time",
                        )?),
                    };
                    super::insert_unique(
                        &mut registry.protoplanetary_disks,
                        disk.id,
                        disk,
                        "protoplanetary disk",
                    )?;
                }
                Some(tag) => return Err(ObjectError(format!("unknown object record '{tag}'"))),
                None => return Err(ObjectError("empty object record".to_owned())),
            }
        }
        registry.validate()?;
        Ok(registry)
    }
}

fn encode_ids<T>(ids: &[T]) -> String
where
    T: Copy,
    u64: From<T>,
{
    ids.iter()
        .map(|id| u64::from(*id).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn float_bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn encode_optional_orbit(orbit: Option<OrbitalElements>) -> String {
    let Some(orbit) = orbit else {
        return "-".to_owned();
    };
    [
        float_bits(orbit.semi_major_axis().meters()),
        float_bits(orbit.eccentricity()),
        float_bits(orbit.inclination().radians()),
        float_bits(orbit.longitude_of_ascending_node().radians()),
        float_bits(orbit.argument_of_periapsis().radians()),
        float_bits(orbit.mean_anomaly_at_epoch().radians()),
        orbit.epoch().seconds().to_string(),
    ]
    .join(",")
}

fn parse_optional_orbit(value: &str) -> Result<Option<OrbitalElements>, ObjectError> {
    if value == "-" {
        return Ok(None);
    }
    let fields: Vec<_> = value.split(',').collect();
    if fields.len() != 7 {
        return Err(ObjectError("malformed companion orbit".to_owned()));
    }
    Ok(Some(OrbitalElements::new(
        Length::from_meters(parse_float_bits(fields[0], "companion semi-major axis")?)?,
        parse_float_bits(fields[1], "companion eccentricity")?,
        Angle::from_radians(parse_float_bits(fields[2], "companion inclination")?)?,
        Angle::from_radians(parse_float_bits(fields[3], "companion ascending node")?)?,
        Angle::from_radians(parse_float_bits(fields[4], "companion periapsis")?)?,
        Angle::from_radians(parse_float_bits(fields[5], "companion mean anomaly")?)?,
        SimDuration::from_seconds(parse_number(fields[6], "companion orbit epoch")?),
    )?))
}

fn parse_orbit(
    fields: &[&str],
    start: usize,
    prefix: &str,
) -> Result<OrbitalElements, ObjectError> {
    OrbitalElements::new(
        Length::from_meters(parse_float_bits(fields[start], "semi-major axis")?)?,
        parse_float_bits(fields[start + 1], "eccentricity")?,
        Angle::from_radians(parse_float_bits(fields[start + 2], "inclination")?)?,
        Angle::from_radians(parse_float_bits(fields[start + 3], "ascending node")?)?,
        Angle::from_radians(parse_float_bits(fields[start + 4], "periapsis argument")?)?,
        Angle::from_radians(parse_float_bits(fields[start + 5], "mean anomaly")?)?,
        SimDuration::from_seconds(parse_number(
            fields[start + 6],
            &format!("{prefix} orbit epoch"),
        )?),
    )
    .map_err(Into::into)
}

fn parse_float_bits(value: &str, field: &str) -> Result<f64, ObjectError> {
    u64::from_str_radix(value, 16)
        .map(f64::from_bits)
        .map_err(|_| ObjectError(format!("invalid {field}")))
}

fn parse_number<T>(value: &str, field: &str) -> Result<T, ObjectError>
where
    T: FromStr,
{
    value
        .parse()
        .map_err(|_| ObjectError(format!("invalid {field}")))
}

fn parse_position(fields: &[&str]) -> Result<Position, ObjectError> {
    Position::from_meters(
        parse_float_bits(fields[0], "position x")?,
        parse_float_bits(fields[1], "position y")?,
        parse_float_bits(fields[2], "position z")?,
    )
    .map_err(ObjectError::from)
}

fn parse_velocity(fields: &[&str]) -> Result<Velocity, ObjectError> {
    Velocity::from_meters_per_second(
        parse_float_bits(fields[0], "velocity x")?,
        parse_float_bits(fields[1], "velocity y")?,
        parse_float_bits(fields[2], "velocity z")?,
    )
    .map_err(ObjectError::from)
}

fn parse_ids<T>(value: &str, constructor: fn(u64) -> T) -> Result<Vec<T>, ObjectError> {
    if value.is_empty() {
        return Ok(Vec::new());
    }
    value
        .split(',')
        .map(|value| parse_number(value, "relationship ID").map(constructor))
        .collect()
}

fn require_field_count(fields: &[&str], expected: usize, kind: &str) -> Result<(), ObjectError> {
    if fields.len() == expected {
        Ok(())
    } else {
        Err(ObjectError(format!(
            "{kind} record has {} fields; expected {expected}",
            fields.len()
        )))
    }
}
