use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, CometCloudId, DebrisDiskId, GalaxyId, MoonId, PlanetId, StarId, StarSystemId,
};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::SimDuration;
use crate::units::{METERS_PER_PARSEC, Position, Velocity};

use super::{
    AsteroidBelt, CometCloud, DebrisDisk, DynamicsMode, GalacticComponent, Galaxy,
    GalaxyCatalogEntry, GalaxyMorphology, Moon, ObjectCounts, ObjectError, ObjectId, ObjectRef,
    PARALLEL_SYSTEM_THRESHOLD, Planet, Star, StarSystem, StellarInteractionKind, StellarStage,
    StellarTransition, generation, insert_unique, sample_stellar_mass, tick_at_or_after,
};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ObjectRegistry {
    pub(crate) galaxies: BTreeMap<GalaxyId, Galaxy>,
    pub(crate) star_systems: BTreeMap<StarSystemId, StarSystem>,
    pub(crate) stars: BTreeMap<StarId, Star>,
    pub(crate) planets: BTreeMap<PlanetId, Planet>,
    pub(crate) moons: BTreeMap<MoonId, Moon>,
    pub(crate) asteroid_belts: BTreeMap<AsteroidBeltId, AsteroidBelt>,
    pub(crate) comet_clouds: BTreeMap<CometCloudId, CometCloud>,
    pub(crate) debris_disks: BTreeMap<DebrisDiskId, DebrisDisk>,
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
            for belt in generated.asteroid_belts {
                insert_unique(&mut registry.asteroid_belts, belt.id, belt, "asteroid belt")?;
            }
            for cloud in generated.comet_clouds {
                insert_unique(&mut registry.comet_clouds, cloud.id, cloud, "comet cloud")?;
            }
            for disk in generated.debris_disks {
                insert_unique(&mut registry.debris_disks, disk.id, disk, "debris disk")?;
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
            asteroid_belts: self.asteroid_belts.len(),
            comet_clouds: self.comet_clouds.len(),
            debris_disks: self.debris_disks.len(),
            gas_reservoirs: usize::from(self.gas_reservoir_solar_masses > 0.0),
        }
    }

    pub const fn gas_reservoir_solar_masses(&self) -> f64 {
        self.gas_reservoir_solar_masses
    }

    pub const fn gas_composition(&self) -> GasComposition {
        self.gas_composition
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
            let mut rng = StableRng::new(derive_object_seed(
                seed,
                sequence,
                "dynamic-star-formation-v1",
            ));
            let available_stellar_mass =
                self.gas_reservoir_solar_masses * STAR_FORMATION_EFFICIENCY;
            let maximum_mass = config.star_mass_solar_max.min(available_stellar_mass);
            if maximum_mass < config.star_mass_solar_min {
                break;
            }
            let mass_solar =
                sample_stellar_mass(&mut rng, config.star_mass_solar_min, maximum_mass);
            let gas_consumed = mass_solar / STAR_FORMATION_EFFICIENCY;
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
            let metallicity = (0.004 + 0.020 * rng.next_f64()).clamp(0.0001, 0.04);
            let mut star = generation::generate_star_with_birth(
                star_id,
                system_id,
                formation_time,
                metallicity,
                mass_solar,
            )?;
            star.position = position;
            star.velocity = velocity;
            self.stars.insert(star_id, star);
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
                    asteroid_belt_ids: Vec::new(),
                    comet_cloud_ids: Vec::new(),
                    debris_disk_ids: Vec::new(),
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

    pub fn asteroid_belt(&self, id: AsteroidBeltId) -> Option<&AsteroidBelt> {
        self.asteroid_belts.get(&id)
    }

    pub fn comet_cloud(&self, id: CometCloudId) -> Option<&CometCloud> {
        self.comet_clouds.get(&id)
    }

    pub fn debris_disk(&self, id: DebrisDiskId) -> Option<&DebrisDisk> {
        self.debris_disks.get(&id)
    }

    pub fn find(&self, id: ObjectId) -> Option<ObjectRef<'_>> {
        match id {
            ObjectId::Galaxy(id) => self.galaxy(id).map(ObjectRef::Galaxy),
            ObjectId::StarSystem(id) => self.star_system(id).map(ObjectRef::StarSystem),
            ObjectId::Star(id) => self.star(id).map(ObjectRef::Star),
            ObjectId::Planet(id) => self.planet(id).map(ObjectRef::Planet),
            ObjectId::Moon(id) => self.moon(id).map(ObjectRef::Moon),
            ObjectId::AsteroidBelt(id) => self.asteroid_belt(id).map(ObjectRef::AsteroidBelt),
            ObjectId::CometCloud(id) => self.comet_cloud(id).map(ObjectRef::CometCloud),
            ObjectId::DebrisDisk(id) => self.debris_disk(id).map(ObjectRef::DebrisDisk),
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

    pub fn asteroid_belts(&self) -> impl Iterator<Item = &AsteroidBelt> {
        self.asteroid_belts.values()
    }

    pub fn comet_clouds(&self) -> impl Iterator<Item = &CometCloud> {
        self.comet_clouds.values()
    }

    pub fn debris_disks(&self) -> impl Iterator<Item = &DebrisDisk> {
        self.debris_disks.values()
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
