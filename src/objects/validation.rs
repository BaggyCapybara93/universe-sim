use std::collections::BTreeSet;
use std::fmt;

use super::{ObjectError, ObjectRegistry, planetary};

impl ObjectRegistry {
    pub fn validate(&self) -> Result<(), ObjectError> {
        if !self.gas_reservoir_solar_masses.is_finite() || self.gas_reservoir_solar_masses < 0.0 {
            return Err(ObjectError(
                "gas reservoir mass must be finite and nonnegative".to_owned(),
            ));
        }
        self.gas_composition
            .validate()
            .map_err(|error| ObjectError(error.to_string()))?;
        let mut cloud_mass = 0.0;
        let mut cloud_mix = self.gas_composition;
        for (cloud_id, cloud) in &self.gas_clouds {
            if *cloud_id != cloud.id || !cloud.mass_solar.is_finite() || cloud.mass_solar < 0.0 {
                return Err(ObjectError(format!("invalid gas cloud {cloud_id}")));
            }
            cloud
                .composition
                .validate()
                .map_err(|error| ObjectError(error.to_string()))?;
            if !self.galaxies.contains_key(&cloud.galaxy_id) {
                return Err(ObjectError(format!(
                    "{cloud_id} references missing galaxy {}",
                    cloud.galaxy_id
                )));
            }
            ensure_unique(&cloud.cluster_ids, "gas-cloud cluster relationship")?;
            for cluster_id in &cloud.cluster_ids {
                let cluster = self.star_clusters.get(cluster_id).ok_or_else(|| {
                    ObjectError(format!("{cloud_id} references missing {cluster_id}"))
                })?;
                if cluster.cloud_id != *cloud_id {
                    return Err(ObjectError(format!(
                        "{cluster_id} does not point back to {cloud_id}"
                    )));
                }
            }
            if cloud.mass_solar > 0.0 {
                cloud_mix = cloud_mix
                    .mixed_with(cloud_mass, cloud.composition, cloud.mass_solar)
                    .map_err(|error| ObjectError(error.to_string()))?;
                cloud_mass += cloud.mass_solar;
            }
        }
        if (cloud_mass - self.gas_reservoir_solar_masses).abs() > 1e-9 * cloud_mass.max(1.0) {
            return Err(ObjectError(
                "gas-cloud mass differs from aggregate gas mass".to_owned(),
            ));
        }
        let expected = [
            cloud_mix.hydrogen_fraction(),
            cloud_mix.helium_fraction(),
            cloud_mix.carbon_fraction(),
            cloud_mix.oxygen_fraction(),
            cloud_mix.iron_fraction(),
            cloud_mix.other_metals_fraction(),
        ];
        let actual = [
            self.gas_composition.hydrogen_fraction(),
            self.gas_composition.helium_fraction(),
            self.gas_composition.carbon_fraction(),
            self.gas_composition.oxygen_fraction(),
            self.gas_composition.iron_fraction(),
            self.gas_composition.other_metals_fraction(),
        ];
        if cloud_mass > 0.0
            && expected
                .iter()
                .zip(actual)
                .any(|(a, b)| (a - b).abs() > 1e-9)
        {
            return Err(ObjectError(
                "gas-cloud composition differs from aggregate gas composition".to_owned(),
            ));
        }
        for (cluster_id, cluster) in &self.star_clusters {
            if *cluster_id != cluster.id {
                return Err(ObjectError(
                    "star-cluster registry key does not match its ID".to_owned(),
                ));
            }
            let cloud = self.gas_clouds.get(&cluster.cloud_id).ok_or_else(|| {
                ObjectError(format!(
                    "{cluster_id} references missing gas cloud {}",
                    cluster.cloud_id
                ))
            })?;
            if !cloud.cluster_ids.contains(cluster_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {cluster_id}",
                    cluster.cloud_id
                )));
            }
            ensure_unique(&cluster.system_ids, "star-cluster system relationship")?;
            if cluster.system_ids.is_empty() {
                return Err(ObjectError(format!("{cluster_id} has no systems")));
            }
            for system_id in &cluster.system_ids {
                let system = self.star_systems.get(system_id).ok_or_else(|| {
                    ObjectError(format!("{cluster_id} references missing {system_id}"))
                })?;
                if system.cluster_id != Some(*cluster_id) || system.galaxy_id != cloud.galaxy_id {
                    return Err(ObjectError(format!(
                        "{system_id} has inconsistent cluster ancestry"
                    )));
                }
            }
        }
        for (galaxy_id, galaxy) in &self.galaxies {
            if *galaxy_id != galaxy.id {
                return Err(ObjectError(
                    "galaxy registry key does not match its ID".to_owned(),
                ));
            }
            ensure_unique(&galaxy.star_system_ids, "galaxy star-system relationship")?;
            for system_id in &galaxy.star_system_ids {
                let system = self.star_systems.get(system_id).ok_or_else(|| {
                    ObjectError(format!("galaxy {galaxy_id} references missing {system_id}"))
                })?;
                if system.galaxy_id != *galaxy_id {
                    return Err(ObjectError(format!(
                        "{system_id} does not point back to {galaxy_id}"
                    )));
                }
            }
        }

        for (system_id, system) in &self.star_systems {
            if *system_id != system.id {
                return Err(ObjectError(
                    "star-system registry key does not match its ID".to_owned(),
                ));
            }
            let cluster_id = system
                .cluster_id
                .ok_or_else(|| ObjectError(format!("{system_id} has no source cluster")))?;
            if !self
                .star_clusters
                .get(&cluster_id)
                .is_some_and(|cluster| cluster.system_ids.contains(system_id))
            {
                return Err(ObjectError(format!(
                    "{system_id} references missing source cluster {cluster_id}"
                )));
            }
            ensure_unique(&system.planet_ids, "star-system planet relationship")?;
            let galaxy = self.galaxies.get(&system.galaxy_id).ok_or_else(|| {
                ObjectError(format!(
                    "{system_id} references missing galaxy {}",
                    system.galaxy_id
                ))
            })?;
            if !galaxy.star_system_ids.contains(system_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {system_id}",
                    system.galaxy_id
                )));
            }
            let star = self.stars.get(&system.star_id).ok_or_else(|| {
                ObjectError(format!(
                    "{system_id} references missing star {}",
                    system.star_id
                ))
            })?;
            if star.system_id != *system_id {
                return Err(ObjectError(format!(
                    "{} does not point back to {system_id}",
                    system.star_id
                )));
            }
            ensure_unique(
                &system.companion_star_ids,
                "star-system companion relationship",
            )?;
            for companion_id in &system.companion_star_ids {
                let companion = self.stars.get(companion_id).ok_or_else(|| {
                    ObjectError(format!(
                        "{system_id} references missing companion {companion_id}"
                    ))
                })?;
                if companion.system_id != *system_id || *companion_id == system.star_id {
                    return Err(ObjectError(format!(
                        "{companion_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for planet_id in &system.planet_ids {
                let planet = self.planets.get(planet_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {planet_id}"))
                })?;
                if planet.system_id != *system_id || planet.star_id != system.star_id {
                    return Err(ObjectError(format!(
                        "{planet_id} has inconsistent parent relationships"
                    )));
                }
            }
            ensure_unique(&system.moon_ids, "star-system moon relationship")?;
            ensure_unique(&system.asteroid_ids, "star-system asteroid relationship")?;
            ensure_unique(&system.comet_ids, "star-system comet relationship")?;
            for moon_id in &system.moon_ids {
                let moon = self.moons.get(moon_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {moon_id}"))
                })?;
                if moon.system_id != *system_id || !system.planet_ids.contains(&moon.planet_id) {
                    return Err(ObjectError(format!(
                        "{moon_id} has inconsistent parent relationships"
                    )));
                }
            }
            ensure_unique(
                &system.asteroid_belt_ids,
                "star-system asteroid-belt relationship",
            )?;
            ensure_unique(
                &system.comet_cloud_ids,
                "star-system comet-cloud relationship",
            )?;
            ensure_unique(
                &system.debris_disk_ids,
                "star-system debris-disk relationship",
            )?;
            ensure_unique(
                &system.protoplanetary_disk_ids,
                "star-system protoplanetary-disk relationship",
            )?;
            for belt_id in &system.asteroid_belt_ids {
                let belt = self.asteroid_belts.get(belt_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {belt_id}"))
                })?;
                if belt.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{belt_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for cloud_id in &system.comet_cloud_ids {
                let cloud = self.comet_clouds.get(cloud_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {cloud_id}"))
                })?;
                if cloud.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{cloud_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for disk_id in &system.debris_disk_ids {
                let disk = self.debris_disks.get(disk_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {disk_id}"))
                })?;
                if disk.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{disk_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for disk_id in &system.protoplanetary_disk_ids {
                let disk = self.protoplanetary_disks.get(disk_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {disk_id}"))
                })?;
                if disk.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{disk_id} has an inconsistent parent relationship"
                    )));
                }
            }
        }

        for (star_id, star) in &self.stars {
            star.birth_composition
                .validate()
                .map_err(|error| ObjectError(error.to_string()))?;
            if !star.metallicity.is_finite()
                || (star.birth_composition.metallicity() - star.metallicity).abs() > 1e-12
            {
                return Err(ObjectError(format!(
                    "{star_id} metallicity differs from its birth composition"
                )));
            }
            if *star_id != star.id {
                return Err(ObjectError(
                    "star registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&star.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{star_id} references missing system {}",
                    star.system_id
                ))
            })?;
            if system.star_id != *star_id && !system.companion_star_ids.contains(star_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {star_id}",
                    star.system_id
                )));
            }
        }

        for (planet_id, planet) in &self.planets {
            if *planet_id != planet.id {
                return Err(ObjectError(
                    "planet registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&planet.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{planet_id} references missing system {}",
                    planet.system_id
                ))
            })?;
            if !system.planet_ids.contains(planet_id) || system.star_id != planet.star_id {
                return Err(ObjectError(format!(
                    "{} has inconsistent child relationship for {planet_id}",
                    planet.system_id
                )));
            }
        }
        for (moon_id, moon) in &self.moons {
            if *moon_id != moon.id {
                return Err(ObjectError(
                    "moon registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&moon.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{moon_id} references missing system {}",
                    moon.system_id
                ))
            })?;
            let planet = self.planets.get(&moon.planet_id).ok_or_else(|| {
                ObjectError(format!(
                    "{moon_id} references missing planet {}",
                    moon.planet_id
                ))
            })?;
            if !system.moon_ids.contains(moon_id)
                || planet.system_id != moon.system_id
                || planet.star_id != system.star_id
            {
                return Err(ObjectError(format!(
                    "{moon_id} has inconsistent child relationship"
                )));
            }
        }
        for (belt_id, belt) in &self.asteroid_belts {
            if *belt_id != belt.id {
                return Err(ObjectError(
                    "asteroid-belt registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&belt.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{belt_id} references missing system {}",
                    belt.system_id
                ))
            })?;
            if !system.asteroid_belt_ids.contains(belt_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {belt_id}",
                    belt.system_id
                )));
            }
            planetary::validate_small_body_population(
                belt.inner_radius,
                belt.outer_radius,
                belt.estimated_body_count,
                "asteroid belt",
            )?;
        }
        for (asteroid_id, asteroid) in &self.asteroids {
            if *asteroid_id != asteroid.id {
                return Err(ObjectError(
                    "asteroid registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&asteroid.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{asteroid_id} references missing system {}",
                    asteroid.system_id
                ))
            })?;
            let belt = self.asteroid_belts.get(&asteroid.belt_id).ok_or_else(|| {
                ObjectError(format!(
                    "{asteroid_id} references missing belt {}",
                    asteroid.belt_id
                ))
            })?;
            if belt.system_id != asteroid.system_id || !system.asteroid_ids.contains(asteroid_id) {
                return Err(ObjectError(format!(
                    "{asteroid_id} has an inconsistent parent relationship"
                )));
            }
        }
        for (cloud_id, cloud) in &self.comet_clouds {
            if *cloud_id != cloud.id {
                return Err(ObjectError(
                    "comet-cloud registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&cloud.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{cloud_id} references missing system {}",
                    cloud.system_id
                ))
            })?;
            if !system.comet_cloud_ids.contains(cloud_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {cloud_id}",
                    cloud.system_id
                )));
            }
            planetary::validate_small_body_population(
                cloud.inner_radius,
                cloud.outer_radius,
                cloud.estimated_body_count,
                "comet cloud",
            )?;
        }
        for (comet_id, comet) in &self.comets {
            if *comet_id != comet.id {
                return Err(ObjectError(
                    "comet registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&comet.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{comet_id} references missing system {}",
                    comet.system_id
                ))
            })?;
            let cloud = self.comet_clouds.get(&comet.cloud_id).ok_or_else(|| {
                ObjectError(format!(
                    "{comet_id} references missing cloud {}",
                    comet.cloud_id
                ))
            })?;
            if cloud.system_id != comet.system_id || !system.comet_ids.contains(comet_id) {
                return Err(ObjectError(format!(
                    "{comet_id} has an inconsistent parent relationship"
                )));
            }
        }
        for (disk_id, disk) in &self.debris_disks {
            if *disk_id != disk.id {
                return Err(ObjectError(
                    "debris-disk registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&disk.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{disk_id} references missing system {}",
                    disk.system_id
                ))
            })?;
            if !system.debris_disk_ids.contains(disk_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {disk_id}",
                    disk.system_id
                )));
            }
            planetary::validate_small_body_population(
                disk.inner_radius,
                disk.outer_radius,
                disk.estimated_body_count,
                "debris disk",
            )?;
        }
        for (disk_id, disk) in &self.protoplanetary_disks {
            if *disk_id != disk.id {
                return Err(ObjectError(
                    "protoplanetary-disk registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&disk.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{disk_id} references missing system {}",
                    disk.system_id
                ))
            })?;
            if !system.protoplanetary_disk_ids.contains(disk_id)
                || disk.inner_radius >= disk.outer_radius
                || disk.initial_mass.kilograms() < disk.remaining_mass.kilograms()
                || disk.formed_at > disk.dissipates_at
            {
                return Err(ObjectError(format!("{disk_id} has invalid disk state")));
            }
        }
        Ok(())
    }
}

fn ensure_unique<T>(ids: &[T], relationship: &str) -> Result<(), ObjectError>
where
    T: Ord + Copy + fmt::Display,
{
    let mut seen = BTreeSet::new();
    for id in ids {
        if !seen.insert(*id) {
            return Err(ObjectError(format!("duplicate {id} in {relationship}")));
        }
    }
    Ok(())
}
