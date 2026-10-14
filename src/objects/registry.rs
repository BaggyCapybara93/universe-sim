use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, AsteroidId, CometCloudId, CometId, DebrisDiskId, GalaxyId, MoonId, PlanetId,
    ProtoplanetaryDiskId, StarId, StarSystemId,
};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::SimDuration;
use crate::units::{METERS_PER_PARSEC, Mass, Position, Velocity};

use super::{
    Asteroid, AsteroidBelt, Comet, CometCloud, DebrisDisk, DynamicsMode, GalacticComponent, Galaxy,
    GalaxyCatalogEntry, GalaxyMorphology, Moon, ObjectCounts, ObjectError, ObjectId, ObjectRef,
    PARALLEL_SYSTEM_THRESHOLD, Planet, ProtoplanetaryDisk, Star, StarSystem,
    StellarInteractionKind, StellarStage, StellarTransition, generation, insert_unique,
    sample_stellar_mass, tick_at_or_after,
};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ObjectRegistry {
    pub(crate) galaxies: BTreeMap<GalaxyId, Galaxy>,
    pub(crate) star_systems: BTreeMap<StarSystemId, StarSystem>,
    pub(crate) stars: BTreeMap<StarId, Star>,
    pub(crate) planets: BTreeMap<PlanetId, Planet>,
    pub(crate) moons: BTreeMap<MoonId, Moon>,
    pub(crate) asteroids: BTreeMap<AsteroidId, Asteroid>,
    pub(crate) comets: BTreeMap<CometId, Comet>,
    pub(crate) asteroid_belts: BTreeMap<AsteroidBeltId, AsteroidBelt>,
    pub(crate) comet_clouds: BTreeMap<CometCloudId, CometCloud>,
    pub(crate) debris_disks: BTreeMap<DebrisDiskId, DebrisDisk>,
    pub(crate) protoplanetary_disks: BTreeMap<ProtoplanetaryDiskId, ProtoplanetaryDisk>,
    pub(crate) gas_reservoir_solar_masses: f64,
    pub(crate) gas_composition: GasComposition,
    pub(crate) formation_sequence: u64,
}

impl ObjectRegistry {
    pub fn generate(seed: u64, config: GenerationConfig) -> Result<Self, ObjectError> {
        let galaxy = GalaxyCatalogEntry::primary(seed, config.galaxy_radius_parsecs);
        Self::generate_in_galaxy(seed, config, galaxy)
    }

    pub(crate) fn generate_in_galaxy(
        seed: u64,
        config: GenerationConfig,
        galaxy: GalaxyCatalogEntry,
    ) -> Result<Self, ObjectError> {
        config
            .validate()
            .map_err(|error| ObjectError(error.to_string()))?;

        let galaxy_id = galaxy.id;
        let mut registry = Self::default();
        let mut system_ids = Vec::with_capacity(config.star_systems as usize);

        let generate =
            |system_index| generation::generate_system(seed, config, galaxy, system_index);
        let generated: Result<Vec<_>, ObjectError> =
            if config.star_systems >= PARALLEL_SYSTEM_THRESHOLD {
                (0..config.star_systems)
                    .into_par_iter()
                    .map(generate)
                    .collect()
            } else {
                (0..config.star_systems).map(generate).collect()
            };
        for generated in generated? {
            let system_id = generated.system.id;
            for star in generated.stars {
                insert_unique(&mut registry.stars, star.id, star, "star")?;
            }
            for planet in generated.planets {
                insert_unique(&mut registry.planets, planet.id, planet, "planet")?;
            }
            for moon in generated.moons {
                insert_unique(&mut registry.moons, moon.id, moon, "moon")?;
            }
            for asteroid in generated.asteroids {
                insert_unique(&mut registry.asteroids, asteroid.id, asteroid, "asteroid")?;
            }
            for comet in generated.comets {
                insert_unique(&mut registry.comets, comet.id, comet, "comet")?;
            }
            for belt in generated.asteroid_belts {
                insert_unique(&mut registry.asteroid_belts, belt.id, belt, "asteroid belt")?;
            }
            for cloud in generated.comet_clouds {
                insert_unique(&mut registry.comet_clouds, cloud.id, cloud, "comet cloud")?;
            }
            for disk in generated.debris_disks {
                insert_unique(&mut registry.debris_disks, disk.id, disk, "debris disk")?;
            }
            for disk in generated.protoplanetary_disks {
                insert_unique(
                    &mut registry.protoplanetary_disks,
                    disk.id,
                    disk,
                    "protoplanetary disk",
                )?;
            }
            insert_unique(
                &mut registry.star_systems,
                system_id,
                generated.system,
                "star system",
            )?;
            system_ids.push(system_id);
        }

        let gas_factor = match galaxy.morphology {
            GalaxyMorphology::Spiral => 2.0,
            GalaxyMorphology::Elliptical => 0.25,
            GalaxyMorphology::Irregular => 3.0,
            GalaxyMorphology::Dwarf => 1.0,
        };
        let gas_metallicity = match galaxy.morphology {
            GalaxyMorphology::Spiral => 0.014,
            GalaxyMorphology::Elliptical => 0.020,
            GalaxyMorphology::Irregular => 0.006,
            GalaxyMorphology::Dwarf => 0.002,
        };
        let galaxy = Galaxy {
            id: galaxy_id,
            position: galaxy.position,
            star_system_ids: system_ids,
        };
        insert_unique(&mut registry.galaxies, galaxy_id, galaxy, "galaxy")?;
        registry.gas_reservoir_solar_masses = f64::from(config.star_systems) * gas_factor;
        registry.gas_composition = GasComposition::from_metallicity(gas_metallicity)
            .map_err(|error| ObjectError(error.to_string()))?;
        registry.validate()?;
        Ok(registry)
    }

    pub fn counts(&self) -> ObjectCounts {
        ObjectCounts {
            galaxies: self.galaxies.len(),
            star_systems: self.star_systems.len(),
            stars: self.stars.len(),
            planets: self.planets.len(),
            moons: self.moons.len(),
            asteroids: self.asteroids.len(),
            comets: self.comets.len(),
            asteroid_belts: self.asteroid_belts.len(),
            comet_clouds: self.comet_clouds.len(),
            debris_disks: self.debris_disks.len(),
            protoplanetary_disks: self.protoplanetary_disks.len(),
            gas_reservoirs: usize::from(self.gas_reservoir_solar_masses > 0.0),
        }
    }

    pub const fn gas_reservoir_solar_masses(&self) -> f64 {
        self.gas_reservoir_solar_masses
    }

    pub const fn gas_composition(&self) -> GasComposition {
        self.gas_composition
    }

    pub(super) fn return_gas(
        &mut self,
        mass: f64,
        composition: GasComposition,
    ) -> Result<(), ObjectError> {
        if mass == 0.0 {
            return Ok(());
        }
        self.gas_composition = self
            .gas_composition
            .mixed_with(self.gas_reservoir_solar_masses, composition, mass)
            .map_err(|error| ObjectError(error.to_string()))?;
        self.gas_reservoir_solar_masses += mass;
        Ok(())
    }

    pub(crate) fn form_stars_between(
        &mut self,
        seed: u64,
        config: GenerationConfig,
        from: SimDuration,
        to: SimDuration,
        tick_duration: SimDuration,
    ) -> Result<Vec<StellarTransition>, ObjectError> {
        const FORMATION_INTERVAL_YEARS: u128 = 10_000_000;
        const STAR_FORMATION_EFFICIENCY: f64 = 0.30;
        let interval = SimDuration::from_years(FORMATION_INTERVAL_YEARS)
            .map_err(|error| ObjectError(error.to_string()))?;
        let galaxy_id = self
            .galaxies
            .keys()
            .next()
            .copied()
            .ok_or_else(|| ObjectError("star formation requires a galaxy".to_owned()))?;
        let mut events = Vec::new();
        self.formation_sequence = self.formation_sequence.max(
            u64::try_from(from.seconds() / interval.seconds())
                .map_err(|_| ObjectError("star-formation sequence overflow".to_owned()))?,
        );
        loop {
            let sequence = self.formation_sequence.saturating_add(1);
            let formation_time = interval
                .checked_mul(sequence)
                .ok_or_else(|| ObjectError("star-formation time overflow".to_owned()))?;
            if formation_time > to {
                break;
            }
            self.formation_sequence = sequence;
            if formation_time <= from {
                continue;
            }
            // Feedback at or before this collapse must precede the new birth.
            events.extend(
                self.advance_stars_to(formation_time, tick_duration)
                    .map_err(|error| ObjectError(error.to_string()))?,
            );
            let mut rng = StableRng::new(derive_object_seed(
                seed,
                sequence,
                "dynamic-star-formation-v1",
            ));
            let available_stellar_mass =
                self.gas_reservoir_solar_masses * STAR_FORMATION_EFFICIENCY;
            let maximum_mass = config.star_mass_solar_max.min(available_stellar_mass);
            if maximum_mass < config.star_mass_solar_min {
                // Skip empty collapse opportunities, but revisit after the next
                // stellar transition, which may replenish an exhausted reservoir.
                let next = self
                    .next_stellar_transition_after(formation_time)
                    .map_err(|error| ObjectError(error.to_string()))?;
                let skip_until = next.map_or(to.seconds(), |time| {
                    time.seconds().saturating_sub(1).min(to.seconds())
                });
                self.formation_sequence = self.formation_sequence.max(
                    u64::try_from(skip_until / interval.seconds())
                        .map_err(|_| ObjectError("star-formation sequence overflow".to_owned()))?,
                );
                continue;
            }
            let mass_solar =
                sample_stellar_mass(&mut rng, config.star_mass_solar_min, maximum_mass);
            let system_id = StarSystemId::new(derive_object_seed(
                seed,
                sequence,
                "formed-star-system-id-v1",
            ));
            let star_id = StarId::new(derive_object_seed(seed, sequence, "formed-star-id-v1"));
            if self.star_systems.contains_key(&system_id) || self.stars.contains_key(&star_id) {
                return Err(ObjectError(
                    "dynamic star-formation ID collision".to_owned(),
                ));
            }
            let radius = config.galaxy_radius_parsecs * METERS_PER_PARSEC * rng.next_f64().sqrt();
            let angle = std::f64::consts::TAU * rng.next_f64();
            let (sin_angle, cos_angle) = angle.sin_cos();
            let [gx, gy, gz] = self
                .galaxies
                .get(&galaxy_id)
                .expect("formation galaxy exists")
                .position
                .meters();
            let position =
                Position::from_meters(gx + radius * cos_angle, gy + radius * sin_angle, gz)?;
            let velocity = Velocity::from_meters_per_second(
                -220_000.0 * sin_angle,
                220_000.0 * cos_angle,
                0.0,
            )?;
            let mut star = generation::generate_star_with_birth(
                star_id,
                system_id,
                formation_time,
                self.gas_composition,
                mass_solar,
            )?;
            star.position = position;
            star.velocity = velocity;
            self.stars.insert(star_id, star);
            let disk_id = ProtoplanetaryDiskId::new(derive_object_seed(
                seed,
                sequence,
                "formed-protoplanetary-disk-id-v1",
            ));
            let mut disk = super::planetary::generate_protoplanetary_disk(
                disk_id,
                system_id,
                Mass::from_solar_masses(mass_solar)?,
                SimDuration::ZERO,
                &mut rng,
            )?;
            disk.formed_at = formation_time;
            disk.initial_mass = Mass::from_solar_masses(
                disk.initial_mass
                    .solar_masses()
                    .min((self.gas_reservoir_solar_masses - mass_solar).max(0.0)),
            )?;
            disk.remaining_mass = disk.initial_mass;
            disk.dissipates_at = formation_time
                .checked_add(
                    SimDuration::from_years(10_000_000)
                        .map_err(|error| ObjectError(error.to_string()))?,
                )
                .ok_or_else(|| ObjectError("protoplanetary disk lifetime overflow".to_owned()))?;
            // Efficiency limits the collapsing mass; unused gas remains in the
            // reservoir instead of disappearing. Fund the disk as well as the star.
            let gas_consumed = mass_solar + disk.initial_mass.solar_masses();
            self.protoplanetary_disks.insert(disk_id, disk);
            self.star_systems.insert(
                system_id,
                StarSystem {
                    id: system_id,
                    galaxy_id,
                    position,
                    velocity,
                    star_id,
                    companion_star_ids: Vec::new(),
                    planet_ids: Vec::new(),
                    moon_ids: Vec::new(),
                    asteroid_ids: Vec::new(),
                    comet_ids: Vec::new(),
                    asteroid_belt_ids: Vec::new(),
                    comet_cloud_ids: Vec::new(),
                    debris_disk_ids: Vec::new(),
                    protoplanetary_disk_ids: vec![disk_id],
                    component: GalacticComponent::Disk,
                },
            );
            self.galaxies
                .get_mut(&galaxy_id)
                .expect("formation galaxy exists")
                .star_system_ids
                .push(system_id);
            self.gas_reservoir_solar_masses =
                (self.gas_reservoir_solar_masses - gas_consumed).max(0.0);
            events.push(StellarTransition {
                tick: tick_at_or_after(formation_time, tick_duration)
                    .map_err(|error| ObjectError(error.to_string()))?,
                star_id,
                from: StellarStage::Protostar,
                to: StellarStage::Protostar,
                impacts: Vec::new(),
                natal_kick_meters_per_second: 0,
                interaction: Some(StellarInteractionKind::StarFormation),
                related_star_id: None,
            });
        }
        self.validate()?;
        Ok(events)
    }

    /// Forms planets from active disks after a deterministic accretion delay.
    pub(crate) fn advance_planetary_formation_to(
        &mut self,
        seed: u64,
        config: GenerationConfig,
        elapsed: SimDuration,
    ) -> Result<(), ObjectError> {
        let accretion_delay =
            SimDuration::from_years(5_000_000).map_err(|error| ObjectError(error.to_string()))?;
        let system_ids: Vec<_> = self.star_systems.keys().copied().collect();
        for system_id in system_ids {
            let system = self.star_systems[&system_id].clone();
            if !system.planet_ids.is_empty() {
                continue;
            }
            let Some(disk_id) = system.protoplanetary_disk_ids.first().copied() else {
                continue;
            };
            let disk = self.protoplanetary_disks[&disk_id].clone();
            let formation_time = disk
                .formed_at
                .checked_add(accretion_delay)
                .ok_or_else(|| ObjectError("planet formation time overflow".to_owned()))?;
            if elapsed < formation_time {
                continue;
            }
            if elapsed >= disk.dissipates_at {
                if let Some(disk) = self.protoplanetary_disks.get_mut(&disk_id) {
                    disk.remaining_mass = Mass::from_kilograms(0.0)?;
                }
                continue;
            }

            let star = self.stars[&system.star_id].clone();
            let mut rng = StableRng::new(derive_object_seed(
                seed,
                disk_id.value(),
                "planet-formation-v1",
            ));
            let count = super::planetary::generate_planet_count(
                &mut rng,
                config,
                star.current_mass,
                star.metallicity,
            );
            let radii = super::planetary::stable_orbital_radii(&mut rng, count, config);
            let material = star.birth_composition;
            let mut planet_ids = Vec::new();
            let mut moon_ids = Vec::new();
            for (index, orbit_au) in radii.into_iter().enumerate() {
                let index = u64::try_from(index).expect("planet index fits in u64");
                let planet_id = PlanetId::new(derive_object_seed(
                    seed,
                    disk_id.value() ^ index,
                    "formed-planet-id-v1",
                ));
                let planet = super::planetary::generate_planet(
                    seed,
                    planet_id,
                    system_id,
                    star.id,
                    system.position,
                    system.velocity,
                    star.current_mass,
                    star.luminosity,
                    orbit_au,
                    material,
                    formation_time,
                    config,
                )?;
                let moon_count = super::planetary::generate_moon_count(&mut rng, planet.mass);
                for moon_index in 0..moon_count {
                    let moon_id = MoonId::new(derive_object_seed(
                        seed,
                        disk_id.value() ^ (index << 8) ^ u64::from(moon_index),
                        "formed-moon-id-v1",
                    ));
                    if let Some(moon) = super::planetary::generate_moon(
                        moon_id,
                        planet.id,
                        system_id,
                        system.position,
                        system.velocity,
                        planet.mass,
                        star.current_mass,
                        star.luminosity,
                        orbit_au,
                        formation_time,
                        &mut rng,
                    )? {
                        self.moons.insert(moon_id, moon);
                        moon_ids.push(moon_id);
                    }
                }
                self.planets.insert(planet_id, planet);
                planet_ids.push(planet_id);
            }
            if let Some(system) = self.star_systems.get_mut(&system_id) {
                system.planet_ids = planet_ids;
                system.moon_ids.extend(moon_ids);
            }
            if let Some(disk) = self.protoplanetary_disks.get_mut(&disk_id) {
                disk.remaining_mass = Mass::from_kilograms(0.0)?;
            }
        }
        self.validate()?;
        Ok(())
    }

    pub fn galaxy(&self, id: GalaxyId) -> Option<&Galaxy> {
        self.galaxies.get(&id)
    }

    pub fn star_system(&self, id: StarSystemId) -> Option<&StarSystem> {
        self.star_systems.get(&id)
    }

    pub fn star(&self, id: StarId) -> Option<&Star> {
        self.stars.get(&id)
    }

    pub fn planet(&self, id: PlanetId) -> Option<&Planet> {
        self.planets.get(&id)
    }

    pub fn moon(&self, id: MoonId) -> Option<&Moon> {
        self.moons.get(&id)
    }

    pub fn asteroid(&self, id: AsteroidId) -> Option<&Asteroid> {
        self.asteroids.get(&id)
    }

    pub fn comet(&self, id: CometId) -> Option<&Comet> {
        self.comets.get(&id)
    }

    pub fn asteroid_belt(&self, id: AsteroidBeltId) -> Option<&AsteroidBelt> {
        self.asteroid_belts.get(&id)
    }

    pub fn comet_cloud(&self, id: CometCloudId) -> Option<&CometCloud> {
        self.comet_clouds.get(&id)
    }

    pub fn debris_disk(&self, id: DebrisDiskId) -> Option<&DebrisDisk> {
        self.debris_disks.get(&id)
    }

    pub fn protoplanetary_disk(&self, id: ProtoplanetaryDiskId) -> Option<&ProtoplanetaryDisk> {
        self.protoplanetary_disks.get(&id)
    }

    pub fn find(&self, id: ObjectId) -> Option<ObjectRef<'_>> {
        match id {
            ObjectId::Galaxy(id) => self.galaxy(id).map(ObjectRef::Galaxy),
            ObjectId::StarSystem(id) => self.star_system(id).map(ObjectRef::StarSystem),
            ObjectId::Star(id) => self.star(id).map(ObjectRef::Star),
            ObjectId::Planet(id) => self.planet(id).map(ObjectRef::Planet),
            ObjectId::Moon(id) => self.moon(id).map(ObjectRef::Moon),
            ObjectId::Asteroid(id) => self.asteroid(id).map(ObjectRef::Asteroid),
            ObjectId::Comet(id) => self.comet(id).map(ObjectRef::Comet),
            ObjectId::AsteroidBelt(id) => self.asteroid_belt(id).map(ObjectRef::AsteroidBelt),
            ObjectId::CometCloud(id) => self.comet_cloud(id).map(ObjectRef::CometCloud),
            ObjectId::DebrisDisk(id) => self.debris_disk(id).map(ObjectRef::DebrisDisk),
            ObjectId::ProtoplanetaryDisk(id) => self
                .protoplanetary_disk(id)
                .map(ObjectRef::ProtoplanetaryDisk),
        }
    }

    pub fn galaxies(&self) -> impl Iterator<Item = &Galaxy> {
        self.galaxies.values()
    }

    pub fn star_systems(&self) -> impl Iterator<Item = &StarSystem> {
        self.star_systems.values()
    }

    pub fn stars(&self) -> impl Iterator<Item = &Star> {
        self.stars.values()
    }

    pub fn planets(&self) -> impl Iterator<Item = &Planet> {
        self.planets.values()
    }

    pub fn moons(&self) -> impl Iterator<Item = &Moon> {
        self.moons.values()
    }

    pub fn asteroids(&self) -> impl Iterator<Item = &Asteroid> {
        self.asteroids.values()
    }

    pub fn comets(&self) -> impl Iterator<Item = &Comet> {
        self.comets.values()
    }

    pub fn asteroid_belts(&self) -> impl Iterator<Item = &AsteroidBelt> {
        self.asteroid_belts.values()
    }

    pub fn comet_clouds(&self) -> impl Iterator<Item = &CometCloud> {
        self.comet_clouds.values()
    }

    pub fn debris_disks(&self) -> impl Iterator<Item = &DebrisDisk> {
        self.debris_disks.values()
    }

    pub fn protoplanetary_disks(&self) -> impl Iterator<Item = &ProtoplanetaryDisk> {
        self.protoplanetary_disks.values()
    }

    pub fn dynamics_mode(&self, system_id: StarSystemId) -> Option<DynamicsMode> {
        self.star_systems.get(&system_id).map(|system| {
            if self.system_requires_nbody(system) {
                DynamicsMode::Interacting
            } else {
                DynamicsMode::Keplerian
            }
        })
    }
}
