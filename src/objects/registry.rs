use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, AsteroidId, CometCloudId, CometId, DebrisDiskId, GalaxyId, GasCloudId, MoonId,
    PlanetId, ProtoplanetaryDiskId, StarClusterId, StarId, StarSystemId,
};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration};
use crate::units::{METERS_PER_PARSEC, Mass, Position, Velocity};

use super::{
    Asteroid, AsteroidBelt, CloudShock, Comet, CometCloud, DebrisDisk, DynamicsMode,
    GalacticComponent, Galaxy, GalaxyCatalogEntry, GalaxyMorphology, GasCloud, GasCloudKind, Moon,
    ObjectCounts, ObjectError, ObjectId, ObjectRef, PARALLEL_SYSTEM_THRESHOLD, Planet,
    ProtoplanetaryDisk, Star, StarCluster, StarSystem, StellarInteractionKind, StellarStage,
    StellarTransition, generation, insert_unique, sample_stellar_mass, tick_at_or_after,
};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ObjectRegistry {
    pub(crate) galaxies: BTreeMap<GalaxyId, Galaxy>,
    pub(crate) gas_clouds: BTreeMap<GasCloudId, GasCloud>,
    pub(crate) star_clusters: BTreeMap<StarClusterId, StarCluster>,
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
        registry.seed_clouds_and_ancestry(seed, galaxy_id, config.star_systems)?;
        registry.validate()?;
        Ok(registry)
    }

    pub fn counts(&self) -> ObjectCounts {
        ObjectCounts {
            galaxies: self.galaxies.len(),
            gas_clouds: self.gas_clouds.len(),
            star_clusters: self.star_clusters.len(),
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
            gas_reservoirs: self
                .gas_clouds
                .values()
                .filter(|cloud| cloud.mass_solar > 0.0)
                .count(),
        }
    }

    pub const fn gas_reservoir_solar_masses(&self) -> f64 {
        self.gas_reservoir_solar_masses
    }

    pub const fn gas_composition(&self) -> GasComposition {
        self.gas_composition
    }

    fn refresh_gas_summary(&mut self) -> Result<(), ObjectError> {
        let mut total = 0.0;
        let mut mixture = self.gas_composition;
        for cloud in self.gas_clouds.values() {
            if cloud.mass_solar > 0.0 {
                mixture = mixture
                    .mixed_with(total, cloud.composition, cloud.mass_solar)
                    .map_err(|error| ObjectError(error.to_string()))?;
                total += cloud.mass_solar;
            }
        }
        self.gas_reservoir_solar_masses = total;
        if total > 0.0 {
            self.gas_composition = mixture;
        }
        Ok(())
    }

    fn seed_clouds_and_ancestry(
        &mut self,
        seed: u64,
        galaxy_id: GalaxyId,
        systems: u32,
    ) -> Result<(), ObjectError> {
        let cloud_count = (systems.div_ceil(8) as usize).clamp(1, 16);
        let total = self.gas_reservoir_solar_masses;
        for index in 0..cloud_count {
            let id = GasCloudId::new(derive_object_seed(seed, index as u64, "gas-cloud-id-v1"));
            let system = self.star_systems.values().nth(index);
            let position =
                system.map_or(self.galaxies[&galaxy_id].position, |system| system.position);
            let kind = if cloud_count > 1 && index % 4 == 3 {
                GasCloudKind::Diffuse
            } else {
                GasCloudKind::Molecular
            };
            let mass_solar = if index + 1 == cloud_count {
                total - total * (index as f64 / cloud_count as f64)
            } else {
                total / cloud_count as f64
            };
            insert_unique(
                &mut self.gas_clouds,
                id,
                GasCloud {
                    id,
                    galaxy_id,
                    kind,
                    position,
                    mass_solar,
                    composition: self.gas_composition,
                    cluster_ids: Vec::new(),
                    pending_shock: None,
                },
                "gas cloud",
            )?;
        }
        let ids: Vec<_> = self.star_systems.keys().copied().collect();
        for (index, system_id) in ids.into_iter().enumerate() {
            let cloud_id = *self
                .gas_clouds
                .keys()
                .nth(index % cloud_count)
                .expect("cloud exists");
            let cluster_id = StarClusterId::new(derive_object_seed(
                seed,
                system_id.value(),
                "initial-star-cluster-id-v1",
            ));
            let position = self.star_systems[&system_id].position;
            insert_unique(
                &mut self.star_clusters,
                cluster_id,
                StarCluster {
                    id: cluster_id,
                    cloud_id,
                    formed_at: SimDuration::ZERO,
                    position,
                    system_ids: vec![system_id],
                    triggering_star_id: None,
                },
                "star cluster",
            )?;
            self.gas_clouds
                .get_mut(&cloud_id)
                .expect("cloud exists")
                .cluster_ids
                .push(cluster_id);
            self.star_systems
                .get_mut(&system_id)
                .expect("system exists")
                .cluster_id = Some(cluster_id);
        }
        for cloud in self.gas_clouds.values_mut() {
            cloud.cluster_ids.sort_unstable();
        }
        Ok(())
    }

    pub(crate) fn attach_legacy_ancestry(&mut self) -> Result<(), ObjectError> {
        if !self.gas_clouds.is_empty() {
            return Ok(());
        }
        let galaxy_id = *self
            .galaxies
            .keys()
            .next()
            .ok_or_else(|| ObjectError("legacy objects have no galaxy".to_owned()))?;
        let cloud_id = GasCloudId::new(derive_object_seed(
            galaxy_id.value(),
            0,
            "legacy-gas-cloud-id-v1",
        ));
        let mut cloud = GasCloud {
            id: cloud_id,
            galaxy_id,
            kind: GasCloudKind::Molecular,
            position: self.galaxies[&galaxy_id].position,
            mass_solar: self.gas_reservoir_solar_masses,
            composition: self.gas_composition,
            cluster_ids: Vec::new(),
            pending_shock: None,
        };
        for system in self.star_systems.values_mut() {
            let cluster_id = StarClusterId::new(derive_object_seed(
                system.id.value(),
                0,
                "legacy-star-cluster-id-v1",
            ));
            if self.star_clusters.contains_key(&cluster_id) {
                return Err(ObjectError("legacy cluster ID collision".to_owned()));
            }
            self.star_clusters.insert(
                cluster_id,
                StarCluster {
                    id: cluster_id,
                    cloud_id,
                    formed_at: SimDuration::ZERO,
                    position: system.position,
                    system_ids: vec![system.id],
                    triggering_star_id: None,
                },
            );
            cloud.cluster_ids.push(cluster_id);
            system.cluster_id = Some(cluster_id);
        }
        self.gas_clouds.insert(cloud_id, cloud);
        Ok(())
    }

    pub(super) fn return_gas(
        &mut self,
        star_id: StarId,
        mass: f64,
        composition: GasComposition,
    ) -> Result<(), ObjectError> {
        if mass == 0.0 {
            return Ok(());
        }
        let system_id = self
            .stars
            .get(&star_id)
            .ok_or_else(|| ObjectError(format!("missing gas-source star {star_id}")))?
            .system_id;
        let cluster_id = self.star_systems[&system_id]
            .cluster_id
            .ok_or_else(|| ObjectError("star has no source cluster".to_owned()))?;
        let cloud_id = self.star_clusters[&cluster_id].cloud_id;
        let cloud = self
            .gas_clouds
            .get_mut(&cloud_id)
            .expect("source cloud exists");
        cloud.composition = cloud
            .composition
            .mixed_with(cloud.mass_solar, composition, mass)
            .map_err(|error| ObjectError(error.to_string()))?;
        cloud.mass_solar += mass;
        self.refresh_gas_summary()
    }

    pub(super) fn apply_supernova_shock(
        &mut self,
        star_id: StarId,
        event_time: SimDuration,
        ejecta_solar_masses: f64,
    ) -> Result<(), ObjectError> {
        if !ejecta_solar_masses.is_finite() || ejecta_solar_masses < 0.0 {
            return Err(ObjectError(
                "supernova ejecta mass must be finite and nonnegative".to_owned(),
            ));
        }
        if ejecta_solar_masses == 0.0 {
            return Ok(());
        }
        let star = self
            .stars
            .get(&star_id)
            .ok_or_else(|| ObjectError(format!("missing shock-source star {star_id}")))?;
        let cluster_id = self.star_systems[&star.system_id]
            .cluster_id
            .ok_or_else(|| ObjectError("shock-source star has no cluster".to_owned()))?;
        let source_cloud_id = self.star_clusters[&cluster_id].cloud_id;
        let source_position = star.position.meters();
        const SHOCK_RADIUS_METERS: f64 = 100.0 * METERS_PER_PARSEC;
        for cloud in self.gas_clouds.values_mut() {
            let [x, y, z] = cloud.position.meters();
            let distance = ((x - source_position[0]).powi(2)
                + (y - source_position[1]).powi(2)
                + (z - source_position[2]).powi(2))
            .sqrt();
            if cloud.id != source_cloud_id && distance > SHOCK_RADIUS_METERS {
                continue;
            }
            let attenuation = if cloud.id == source_cloud_id {
                1.0
            } else {
                (1.0 - distance / SHOCK_RADIUS_METERS).powi(2)
            };
            let strength = (attenuation * ejecta_solar_masses
                / (ejecta_solar_masses + cloud.mass_solar))
                .clamp(0.0, 1.0);
            if strength <= 0.0 {
                continue;
            }
            let combined = cloud.pending_shock.map_or(strength, |previous| {
                1.0 - (1.0 - previous.strength) * (1.0 - strength)
            });
            cloud.pending_shock = Some(CloudShock {
                source_star_id: star_id,
                arrived_at: event_time,
                strength: combined,
            });
            cloud.kind = GasCloudKind::Molecular;
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn form_stars_between(
        &mut self,
        seed: u64,
        config: GenerationConfig,
        from: SimDuration,
        to: SimDuration,
        tick_duration: SimDuration,
    ) -> Result<Vec<StellarTransition>, ObjectError> {
        self.form_stars_between_with_encounters(seed, config, from, to, tick_duration, &[])
    }

    pub(crate) fn form_stars_between_with_encounters(
        &mut self,
        seed: u64,
        config: GenerationConfig,
        from: SimDuration,
        to: SimDuration,
        tick_duration: SimDuration,
        encounters: &[crate::events::GalaxyEncounter],
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
            let tidal_strength = encounters
                .iter()
                .filter(|encounter| {
                    encounter.involves(galaxy_id)
                        && encounter.occurred_at() <= formation_time
                        && formation_time.seconds() - encounter.occurred_at().seconds()
                            <= 200_000_000 * SECONDS_PER_YEAR
                })
                .map(|encounter| match encounter.kind() {
                    crate::events::GalaxyEncounterKind::ClosePass => 0.5_f64,
                    crate::events::GalaxyEncounterKind::Merger => 1.0_f64,
                })
                .fold(0.0_f64, f64::max);
            let tidal_multiplier = 1.0 + tidal_strength;
            let source = self
                .gas_clouds
                .values()
                .filter(|cloud| cloud.kind == GasCloudKind::Molecular)
                .max_by(|a, b| {
                    let budget = |cloud: &GasCloud| {
                        cloud.mass_solar
                            * tidal_multiplier
                            * (STAR_FORMATION_EFFICIENCY
                                + cloud
                                    .pending_shock
                                    .map_or(0.0, |shock| 0.30 * shock.strength))
                    };
                    budget(a)
                        .total_cmp(&budget(b))
                        .then_with(|| a.id.cmp(&b.id))
                })
                .cloned();
            let available_stellar_mass = source.as_ref().map_or(0.0, |cloud| {
                cloud.mass_solar
                    * tidal_multiplier
                    * (STAR_FORMATION_EFFICIENCY
                        + cloud
                            .pending_shock
                            .map_or(0.0, |shock| 0.30 * shock.strength))
            });
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
            let source = source.expect("eligible source cloud exists");
            let cluster_id = StarClusterId::new(derive_object_seed(
                seed,
                sequence,
                "formed-star-cluster-id-v1",
            ));
            if self.star_clusters.contains_key(&cluster_id) {
                return Err(ObjectError("dynamic star-cluster ID collision".to_owned()));
            }
            let [cx, cy, cz] = source.position.meters();
            let [gx, gy, _gz] = self
                .galaxies
                .get(&galaxy_id)
                .expect("formation galaxy exists")
                .position
                .meters();
            let angle = (cy - gy).atan2(cx - gx);
            let (sin_angle, cos_angle) = angle.sin_cos();
            let velocity = Velocity::from_meters_per_second(
                -220_000.0 * sin_angle,
                220_000.0 * cos_angle,
                0.0,
            )?;
            let mut cluster = StarCluster {
                id: cluster_id,
                cloud_id: source.id,
                formed_at: formation_time,
                position: source.position,
                system_ids: Vec::new(),
                triggering_star_id: source.pending_shock.map(|shock| shock.source_star_id),
            };
            let mut remaining_budget = available_stellar_mass;
            let member_limit = if source.pending_shock.is_some() || tidal_strength > 0.0 {
                8
            } else {
                4
            };
            for member in 0..member_limit {
                let remaining_cloud = self.gas_clouds[&source.id].mass_solar;
                let maximum_mass = config
                    .star_mass_solar_max
                    .min(remaining_budget)
                    .min(remaining_cloud);
                if maximum_mass < config.star_mass_solar_min {
                    break;
                }
                let key = sequence.wrapping_mul(8).wrapping_add(member);
                let mut rng =
                    StableRng::new(derive_object_seed(seed, key, "dynamic-star-formation-v3"));
                let mass_solar =
                    sample_stellar_mass(&mut rng, config.star_mass_solar_min, maximum_mass);
                let system_id =
                    StarSystemId::new(derive_object_seed(seed, key, "formed-star-system-id-v3"));
                let star_id = StarId::new(derive_object_seed(seed, key, "formed-star-id-v3"));
                let disk_id = ProtoplanetaryDiskId::new(derive_object_seed(
                    seed,
                    key,
                    "formed-protoplanetary-disk-id-v3",
                ));
                if self.star_systems.contains_key(&system_id)
                    || self.stars.contains_key(&star_id)
                    || self.protoplanetary_disks.contains_key(&disk_id)
                {
                    return Err(ObjectError(
                        "dynamic star-formation ID collision".to_owned(),
                    ));
                }
                let offset = 0.1 * METERS_PER_PARSEC * rng.next_f64().sqrt();
                let direction = std::f64::consts::TAU * rng.next_f64();
                let position = Position::from_meters(
                    cx + offset * direction.cos(),
                    cy + offset * direction.sin(),
                    cz,
                )?;
                let mut star = generation::generate_star_with_birth(
                    star_id,
                    system_id,
                    formation_time,
                    source.composition,
                    mass_solar,
                )?;
                star.position = position;
                star.velocity = velocity;
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
                        .min((remaining_cloud - mass_solar).max(0.0)),
                )?;
                disk.remaining_mass = disk.initial_mass;
                disk.dissipates_at = formation_time
                    .checked_add(
                        SimDuration::from_years(10_000_000)
                            .map_err(|error| ObjectError(error.to_string()))?,
                    )
                    .ok_or_else(|| {
                        ObjectError("protoplanetary disk lifetime overflow".to_owned())
                    })?;
                let gas_consumed = mass_solar + disk.initial_mass.solar_masses();
                self.gas_clouds
                    .get_mut(&source.id)
                    .expect("source cloud exists")
                    .mass_solar = (remaining_cloud - gas_consumed).max(0.0);
                self.refresh_gas_summary()?;
                remaining_budget = (remaining_budget - gas_consumed).max(0.0);
                self.stars.insert(star_id, star);
                self.protoplanetary_disks.insert(disk_id, disk);
                self.star_systems.insert(
                    system_id,
                    StarSystem {
                        id: system_id,
                        galaxy_id,
                        cluster_id: Some(cluster_id),
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
                cluster.system_ids.push(system_id);
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
            let cloud = self
                .gas_clouds
                .get_mut(&source.id)
                .expect("source cloud exists");
            cloud.pending_shock = None;
            cloud.cluster_ids.push(cluster_id);
            cloud.cluster_ids.sort_unstable();
            self.star_clusters.insert(cluster_id, cluster);
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

    pub fn gas_cloud(&self, id: GasCloudId) -> Option<&GasCloud> {
        self.gas_clouds.get(&id)
    }
    pub fn star_cluster(&self, id: StarClusterId) -> Option<&StarCluster> {
        self.star_clusters.get(&id)
    }
    pub fn gas_clouds(&self) -> impl Iterator<Item = &GasCloud> {
        self.gas_clouds.values()
    }
    pub fn star_clusters(&self) -> impl Iterator<Item = &StarCluster> {
        self.star_clusters.values()
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
            ObjectId::GasCloud(id) => self.gas_cloud(id).map(ObjectRef::GasCloud),
            ObjectId::StarCluster(id) => self.star_cluster(id).map(ObjectRef::StarCluster),
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
