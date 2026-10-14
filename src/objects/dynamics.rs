use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;

use crate::id::{PlanetId, StarId, StarSystemId};
use crate::orbital::{
    GRAVITATIONAL_CONSTANT, GravityBody, OrbitalElements, leapfrog_step, propagate_kepler,
};
use crate::time::{SECONDS_PER_YEAR, SimDuration, TimeError};
use crate::units::{Length, METERS_PER_ASTRONOMICAL_UNIT, Mass, Position, Temperature, Velocity};

use super::{
    Moon, ObjectError, ObjectRegistry, Planet, PlanetAtmosphere, PlanetImpact, PlanetImpactKind,
    PlanetStatus, StarSystem, StellarInteractionKind, StellarStage, StellarTransition,
};

impl ObjectRegistry {
    pub(super) fn refresh_planet_temperatures(&mut self) -> Result<(), TimeError> {
        let mut updates = Vec::new();
        for planet in self.planets.values() {
            let Some(star) = self.stars.get(&planet.star_id) else {
                continue;
            };
            let temperature =
                Temperature::from_kelvin(super::planetary::equilibrium_temperature_kelvin(
                    star.luminosity,
                    planet.orbital_radius.astronomical_units(),
                ))
                .map_err(|_| TimeError::Overflow)?;
            updates.push((planet.id, temperature));
        }
        for (planet_id, temperature) in updates {
            self.planets
                .get_mut(&planet_id)
                .expect("refreshed planet exists")
                .equilibrium_temperature = temperature;
        }
        Ok(())
    }

    pub(super) fn apply_direct_planet_collisions(
        &mut self,
        star_id: StarId,
        stage: StellarStage,
    ) -> Vec<PlanetImpact> {
        let Some(star) = self.stars.get(&star_id) else {
            return Vec::new();
        };
        let collision_radius_au = super::stellar::stellar_radius_solar(star.current_mass, stage)
            * super::SOLAR_RADIUS_METERS
            / METERS_PER_ASTRONOMICAL_UNIT;
        let mut impacts = Vec::new();
        for planet in self
            .planets
            .values_mut()
            .filter(|planet| planet.star_id == star_id && planet.status != PlanetStatus::Destroyed)
        {
            let periapsis_au = planet
                .orbital_elements
                .semi_major_axis()
                .astronomical_units()
                * (1.0 - planet.orbital_elements.eccentricity());
            if periapsis_au <= collision_radius_au {
                planet.status = PlanetStatus::Destroyed;
                impacts.push(PlanetImpact {
                    planet_id: planet.id,
                    kind: PlanetImpactKind::Destroyed,
                });
            }
        }
        impacts
    }

    pub(super) fn resolve_binary_interactions(
        &mut self,
        trigger_ticks: &BTreeMap<StarId, u64>,
        elapsed: SimDuration,
    ) -> Result<Vec<StellarTransition>, TimeError> {
        let mut candidates = Vec::new();
        let affected_systems: BTreeSet<_> = trigger_ticks
            .keys()
            .filter_map(|id| self.stars.get(id).map(|star| star.system_id))
            .collect();
        for system_id in affected_systems {
            let system = &self.star_systems[&system_id];
            for companion_id in &system.companion_star_ids {
                let Some(companion) = self.stars.get(companion_id) else {
                    continue;
                };
                let Some(orbit) = companion.orbital_elements else {
                    continue;
                };
                let tick = trigger_ticks
                    .get(&system.star_id)
                    .or_else(|| trigger_ticks.get(companion_id));
                if let Some(tick) = tick {
                    candidates.push((system.id, system.star_id, *companion_id, orbit, *tick));
                }
            }
        }
        let mut outcomes = Vec::new();
        for (system_id, primary_id, companion_id, orbit, tick) in candidates {
            let (Some(primary), Some(companion)) =
                (self.stars.get(&primary_id), self.stars.get(&companion_id))
            else {
                continue;
            };
            if primary.current_mass.kilograms() == 0.0 || companion.current_mass.kilograms() == 0.0
            {
                continue;
            }
            let separation_au =
                orbit.semi_major_axis().astronomical_units() * (1.0 - orbit.eccentricity());
            let primary_radius_au =
                super::stellar::stellar_radius_solar(primary.current_mass, primary.stage)
                    * super::SOLAR_RADIUS_METERS
                    / METERS_PER_ASTRONOMICAL_UNIT;
            let companion_radius_au =
                super::stellar::stellar_radius_solar(companion.current_mass, companion.stage)
                    * super::SOLAR_RADIUS_METERS
                    / METERS_PER_ASTRONOMICAL_UNIT;
            let primary_lobe = super::stellar::roche_lobe_radius(
                separation_au,
                primary.current_mass.solar_masses() / companion.current_mass.solar_masses(),
            );
            let companion_lobe = super::stellar::roche_lobe_radius(
                separation_au,
                companion.current_mass.solar_masses() / primary.current_mass.solar_masses(),
            );
            let primary_overflow = primary_radius_au > primary_lobe;
            let companion_overflow = companion_radius_au > companion_lobe;
            if !primary_overflow
                && !companion_overflow
                && primary_radius_au + companion_radius_au < separation_au
            {
                continue;
            }
            let mass_ratio = primary
                .current_mass
                .solar_masses()
                .max(companion.current_mass.solar_masses())
                / primary
                    .current_mass
                    .solar_masses()
                    .min(companion.current_mass.solar_masses());
            if primary_radius_au + companion_radius_au >= separation_au
                || (primary_overflow || companion_overflow) && mass_ratio > 3.0
            {
                let previous_stage = primary.stage;
                let total_mass =
                    primary.current_mass.solar_masses() + companion.current_mass.solar_masses();
                let composition = primary
                    .birth_composition
                    .mixed_with(
                        primary.current_mass.solar_masses(),
                        companion.birth_composition,
                        companion.current_mass.solar_masses(),
                    )
                    .map_err(|_| TimeError::Overflow)?;
                let merged_mass = 0.9
                    * (primary.current_mass.solar_masses() + companion.current_mass.solar_masses());
                let merged_mass =
                    Mass::from_solar_masses(merged_mass).map_err(|_| TimeError::Overflow)?;
                let primary = self
                    .stars
                    .get_mut(&primary_id)
                    .expect("merger primary exists");
                primary.initial_mass = merged_mass;
                primary.current_mass = merged_mass;
                primary.birth_composition = composition;
                primary.metallicity = composition.metallicity();
                primary.epoch_age = SimDuration::ZERO;
                primary.formed_at = elapsed;
                primary.age = SimDuration::ZERO;
                primary.stage = StellarStage::MainSequence;
                let (_, luminosity, temperature) = super::stellar::stellar_properties(
                    merged_mass,
                    primary.metallicity,
                    StellarStage::MainSequence,
                )
                .map_err(|_| TimeError::Overflow)?;
                primary.luminosity = luminosity;
                primary.temperature = temperature;
                self.stars.remove(&companion_id);
                self.return_gas(total_mass - merged_mass.solar_masses(), composition)
                    .map_err(|_| TimeError::Overflow)?;
                if let Some(system) = self.star_systems.get_mut(&system_id) {
                    system.companion_star_ids.retain(|id| *id != companion_id);
                }
                outcomes.push(StellarTransition {
                    tick,
                    star_id: primary_id,
                    from: previous_stage,
                    to: StellarStage::MainSequence,
                    impacts: Vec::new(),
                    natal_kick_meters_per_second: 0,
                    interaction: Some(StellarInteractionKind::StellarMerger),
                    related_star_id: Some(companion_id),
                });
                continue;
            }
            let (donor_id, accretor_id) = if primary_overflow {
                (primary_id, companion_id)
            } else {
                (companion_id, primary_id)
            };
            let donor = self.stars.get(&donor_id).expect("binary donor exists");
            let donor_stage = donor.stage;
            let donor_old_solar = donor.current_mass.solar_masses();
            let donor_composition = donor.birth_composition;
            let accretor_old_solar = self
                .stars
                .get(&accretor_id)
                .expect("binary accretor exists")
                .current_mass
                .solar_masses();
            let transferred_solar = 0.05 * donor_old_solar;
            let donor_mass = Mass::from_solar_masses(donor_old_solar - transferred_solar)
                .map_err(|_| TimeError::Overflow)?;
            let accretor_mass =
                Mass::from_solar_masses(accretor_old_solar + 0.70 * transferred_solar)
                    .map_err(|_| TimeError::Overflow)?;
            self.stars
                .get_mut(&donor_id)
                .expect("binary donor exists")
                .current_mass = donor_mass;
            self.stars
                .get_mut(&accretor_id)
                .expect("binary accretor exists")
                .current_mass = accretor_mass;
            self.return_gas(
                donor_old_solar + accretor_old_solar
                    - donor_mass.solar_masses()
                    - accretor_mass.solar_masses(),
                donor_composition,
            )
            .map_err(|_| TimeError::Overflow)?;
            let old_product = donor_old_solar * accretor_old_solar;
            let new_product = donor_mass.solar_masses() * accretor_mass.solar_masses();
            let axis_scale = (old_product / new_product).powi(2).clamp(0.25, 4.0);
            if let Some(companion) = self.stars.get_mut(&companion_id)
                && let Some(existing) = companion.orbital_elements
            {
                companion.orbital_elements = Some(
                    OrbitalElements::new(
                        Length::from_meters(existing.semi_major_axis().meters() * axis_scale)
                            .map_err(|_| TimeError::Overflow)?,
                        existing.eccentricity() * 0.5,
                        existing.inclination(),
                        existing.longitude_of_ascending_node(),
                        existing.argument_of_periapsis(),
                        existing.mean_anomaly_at_epoch(),
                        existing.epoch(),
                    )
                    .map_err(|_| TimeError::Overflow)?,
                );
            }
            outcomes.push(StellarTransition {
                tick,
                star_id: donor_id,
                from: donor_stage,
                to: donor_stage,
                impacts: Vec::new(),
                natal_kick_meters_per_second: 0,
                interaction: Some(StellarInteractionKind::MassTransfer),
                related_star_id: Some(accretor_id),
            });
        }
        Ok(outcomes)
    }

    pub(super) fn apply_stellar_impacts(
        &mut self,
        star_id: StarId,
        previous_mass: Mass,
        current_mass: Mass,
        stage: StellarStage,
    ) -> Vec<PlanetImpact> {
        let retained_mass = current_mass.kilograms() / previous_mass.kilograms();
        let mut impacts = Vec::new();
        for planet in self
            .planets
            .values_mut()
            .filter(|planet| planet.star_id == star_id)
        {
            if matches!(
                planet.status,
                PlanetStatus::Unbound | PlanetStatus::Destroyed
            ) {
                continue;
            }
            let periapsis_au = planet
                .orbital_elements
                .semi_major_axis()
                .astronomical_units()
                * (1.0 - planet.orbital_elements.eccentricity());
            let outcome = match stage {
                StellarStage::AsymptoticGiantBranch if periapsis_au <= 2.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::AsymptoticGiantBranch if periapsis_au <= 5.0 => {
                    Some(PlanetImpactKind::Irradiated)
                }
                StellarStage::Giant if periapsis_au <= 0.5 => Some(PlanetImpactKind::Destroyed),
                StellarStage::Giant if periapsis_au <= 2.0 => Some(PlanetImpactKind::Irradiated),
                StellarStage::Supergiant if periapsis_au <= 3.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::Supergiant if periapsis_au <= 10.0 => {
                    Some(PlanetImpactKind::Irradiated)
                }
                StellarStage::PulsationalPairInstability if periapsis_au <= 20.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::PulsationalPairInstability => Some(PlanetImpactKind::Irradiated),
                StellarStage::CoreCollapseSupernova if periapsis_au <= 20.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::CoreCollapseSupernova => Some(PlanetImpactKind::Irradiated),
                StellarStage::Disrupted if periapsis_au <= 20.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::Disrupted => Some(PlanetImpactKind::Unbound),
                StellarStage::NeutronStar | StellarStage::BlackHole if periapsis_au <= 5.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::PlanetaryNebula
                | StellarStage::WhiteDwarf
                | StellarStage::NeutronStar
                | StellarStage::BlackHole
                    if retained_mass < 0.5 =>
                {
                    Some(PlanetImpactKind::Unbound)
                }
                StellarStage::NeutronStar | StellarStage::BlackHole => {
                    Some(PlanetImpactKind::Irradiated)
                }
                _ => None,
            };
            let Some(kind) = outcome else {
                continue;
            };
            let new_status = match kind {
                PlanetImpactKind::Irradiated => PlanetStatus::Irradiated,
                PlanetImpactKind::Unbound => PlanetStatus::Unbound,
                PlanetImpactKind::Destroyed => PlanetStatus::Destroyed,
            };
            if planet.status == new_status {
                continue;
            }
            planet.status = new_status;
            impacts.push(PlanetImpact {
                planet_id: planet.id,
                kind,
            });
        }
        impacts
    }

    pub(crate) fn advance_orbits_to(&mut self, elapsed: SimDuration) -> Result<(), ObjectError> {
        self.evolve_planetary_systems(elapsed)?;
        let system_ids: Vec<_> = self.star_systems.keys().copied().collect();
        let mut integrated_systems = BTreeSet::new();
        for system_id in system_ids {
            if self.integrate_system_from_epoch(system_id, elapsed)? {
                integrated_systems.insert(system_id);
            }
        }
        let companion_inputs: Vec<_> = self
            .stars
            .values()
            .filter_map(|star| {
                (!integrated_systems.contains(&star.system_id))
                    .then_some(star.orbital_elements)
                    .flatten()
                    .map(|orbit| (star.id, star.system_id, orbit))
            })
            .collect();
        let propagate_companion =
            |(star_id, system_id, orbit): &(StarId, StarSystemId, OrbitalElements)| {
                let system = &self.star_systems[system_id];
                let primary = &self.stars[&system.star_id];
                let companion = &self.stars[star_id];
                let interior_mass: f64 = system
                    .companion_star_ids
                    .iter()
                    .filter_map(|id| self.stars.get(id))
                    .filter_map(|star| star.orbital_elements.map(|candidate| (star, candidate)))
                    .filter(|(_, candidate)| {
                        candidate.semi_major_axis().meters() < orbit.semi_major_axis().meters()
                    })
                    .map(|(star, _)| star.current_mass.kilograms())
                    .sum();
                let orbital_mass = Mass::from_kilograms(
                    primary.current_mass.kilograms()
                        + interior_mass
                        + companion.current_mass.kilograms(),
                )?;
                if orbital_mass.kilograms() == 0.0 {
                    return Ok((*star_id, companion.position, companion.velocity));
                }
                let (position, velocity) = propagate_kepler(
                    *orbit,
                    orbital_mass,
                    elapsed,
                    system.position,
                    system.velocity,
                )?;
                Ok((*star_id, position, velocity))
            };
        let companion_updates: Result<Vec<_>, ObjectError> =
            if companion_inputs.len() >= super::PARALLEL_BODY_THRESHOLD {
                companion_inputs
                    .par_iter()
                    .map(propagate_companion)
                    .collect()
            } else {
                companion_inputs.iter().map(propagate_companion).collect()
            };
        for (star_id, position, velocity) in companion_updates? {
            let star = self
                .stars
                .get_mut(&star_id)
                .expect("propagated star exists");
            star.position = position;
            star.velocity = velocity;
        }

        let planets: Vec<_> = self
            .planets
            .values()
            .filter(|planet| {
                !integrated_systems.contains(&planet.system_id)
                    && !matches!(
                        planet.status,
                        PlanetStatus::Unbound | PlanetStatus::Destroyed
                    )
            })
            .collect();
        let propagate = |planet: &&Planet| -> Result<_, ObjectError> {
            let system = &self.star_systems[&planet.system_id];
            let star = &self.stars[&planet.star_id];
            let (position, velocity) = propagate_kepler(
                planet.orbital_elements,
                star.current_mass,
                elapsed,
                system.position,
                system.velocity,
            )?;
            Ok((planet.id, position, velocity))
        };
        let updates: Result<Vec<_>, ObjectError> =
            if planets.len() >= super::PARALLEL_BODY_THRESHOLD {
                planets.par_iter().map(propagate).collect()
            } else {
                planets.iter().map(propagate).collect()
            };
        for (planet_id, position, velocity) in updates? {
            let planet = self
                .planets
                .get_mut(&planet_id)
                .expect("propagated planet exists");
            planet.position = position;
            planet.velocity = velocity;
        }
        let moons: Vec<_> = self
            .moons
            .values()
            .filter_map(|moon| {
                let planet = self.planets.get(&moon.planet_id)?;
                (!matches!(
                    planet.status,
                    PlanetStatus::Unbound | PlanetStatus::Destroyed
                ))
                .then_some((moon, planet.position, planet.velocity, planet.mass))
            })
            .collect();
        let propagate_moon = |(moon, planet_position, planet_velocity, planet_mass): &(
            &Moon,
            Position,
            Velocity,
            Mass,
        )| {
            let (position, velocity) = propagate_kepler(
                moon.orbital_elements,
                *planet_mass,
                elapsed,
                *planet_position,
                *planet_velocity,
            )?;
            Ok((moon.id, position, velocity))
        };
        let moon_updates: Result<Vec<_>, ObjectError> =
            if moons.len() >= super::PARALLEL_BODY_THRESHOLD {
                moons.par_iter().map(propagate_moon).collect()
            } else {
                moons.iter().map(propagate_moon).collect()
            };
        for (moon_id, position, velocity) in moon_updates? {
            let moon = self
                .moons
                .get_mut(&moon_id)
                .expect("propagated moon exists");
            moon.position = position;
            moon.velocity = velocity;
        }
        for asteroid in self.asteroids.values_mut() {
            let system = &self.star_systems[&asteroid.system_id];
            let star = &self.stars[&system.star_id];
            let (position, velocity) = propagate_kepler(
                asteroid.orbital_elements,
                star.current_mass,
                elapsed,
                system.position,
                system.velocity,
            )?;
            asteroid.position = position;
            asteroid.velocity = velocity;
        }
        for comet in self.comets.values_mut() {
            let system = &self.star_systems[&comet.system_id];
            let star = &self.stars[&system.star_id];
            let (position, velocity) = propagate_kepler(
                comet.orbital_elements,
                star.current_mass,
                elapsed,
                system.position,
                system.velocity,
            )?;
            comet.position = position;
            comet.velocity = velocity;
        }
        Ok(())
    }

    fn evolve_planetary_systems(&mut self, elapsed: SimDuration) -> Result<(), ObjectError> {
        for planet in self.planets.values_mut() {
            let age_seconds = elapsed.seconds().saturating_sub(planet.formed_at.seconds()) as f64;
            let age_years = age_seconds / SECONDS_PER_YEAR as f64;
            let migrated_au = (planet.initial_orbital_radius.astronomical_units()
                + planet.migration_rate_au_per_year * age_years)
                .clamp(0.03, 10_000.0);
            planet.orbital_elements = planet
                .orbital_elements
                .with_semi_major_axis(Length::from_astronomical_units(migrated_au)?)?;
            planet.orbital_radius = Length::from_astronomical_units(migrated_au)?;
            if let Some(star) = self.stars.get(&planet.star_id) {
                let flux = star.luminosity.solar_luminosities() / migrated_au.powi(2).max(0.001);
                planet.atmosphere_loss_fraction = (age_years / 1.0e9 * flux * 0.08).clamp(0.0, 1.0);
                planet.atmosphere = if planet.atmosphere_loss_fraction >= 0.90 {
                    PlanetAtmosphere::None
                } else if planet.atmosphere_loss_fraction >= 0.45 {
                    PlanetAtmosphere::Thin
                } else {
                    planet.atmosphere
                };
                planet.equilibrium_temperature = Temperature::from_kelvin(
                    super::planetary::equilibrium_temperature_kelvin(star.luminosity, migrated_au),
                )?;
                planet.in_habitable_zone =
                    super::planetary::is_in_habitable_zone(star.luminosity, migrated_au);
            }
        }

        let system_ids: Vec<_> = self.star_systems.keys().copied().collect();
        for system_id in system_ids {
            let mut ids = self.star_systems[&system_id].planet_ids.clone();
            ids.retain(|id| self.planets.contains_key(id));
            ids.sort_by(|left, right| {
                self.planets[left]
                    .orbital_radius
                    .meters()
                    .total_cmp(&self.planets[right].orbital_radius.meters())
            });
            for pair in ids.windows(2) {
                let (inner, outer) = (&self.planets[&pair[0]], &self.planets[&pair[1]]);
                let inner_apo = inner.orbital_elements.semi_major_axis().meters()
                    * (1.0 + inner.orbital_elements.eccentricity());
                let outer_peri = outer.orbital_elements.semi_major_axis().meters()
                    * (1.0 - outer.orbital_elements.eccentricity());
                if outer.orbital_radius.meters() / inner.orbital_radius.meters() >= 1.5
                    && inner_apo < outer_peri
                {
                    continue;
                }
                let inner_axis = inner.orbital_radius.meters() * 0.99;
                let outer_axis = outer.orbital_radius.meters() * 1.01;
                let inner_e = (inner.orbital_elements.eccentricity() + 0.03).min(0.85);
                let outer_e = (outer.orbital_elements.eccentricity() + 0.03).min(0.85);
                let inner_orbit = inner
                    .orbital_elements
                    .with_semi_major_axis(Length::from_meters(inner_axis)?)?
                    .with_eccentricity(inner_e)?;
                let outer_orbit = outer
                    .orbital_elements
                    .with_semi_major_axis(Length::from_meters(outer_axis)?)?
                    .with_eccentricity(outer_e)?;
                let inner_planet = self.planets.get_mut(&pair[0]).expect("inner planet exists");
                inner_planet.orbital_elements = inner_orbit;
                inner_planet.orbital_radius = Length::from_meters(inner_axis)?;
                inner_planet.scattering_count = inner_planet.scattering_count.saturating_add(1);
                let outer_planet = self.planets.get_mut(&pair[1]).expect("outer planet exists");
                outer_planet.orbital_elements = outer_orbit;
                outer_planet.orbital_radius = Length::from_meters(outer_axis)?;
                outer_planet.scattering_count = outer_planet.scattering_count.saturating_add(1);
            }
        }

        for moon in self.moons.values_mut() {
            let Some(planet) = self.planets.get(&moon.planet_id) else {
                continue;
            };
            let age_seconds = elapsed
                .seconds()
                .saturating_sub(moon.formation_history.formed_at().seconds())
                as f64;
            let age_years = age_seconds / SECONDS_PER_YEAR as f64;
            let target_axis = (moon.initial_orbital_radius.meters()
                + moon.tidal_migration_rate_m_per_year * age_years)
                .max(planet.radius.meters() * 2.0);
            moon.orbital_elements = moon
                .orbital_elements
                .with_semi_major_axis(Length::from_meters(target_axis)?)?;
            moon.orbital_radius = Length::from_meters(target_axis)?;
            if let Some(star) = self.stars.get(&planet.star_id) {
                moon.equilibrium_temperature =
                    Temperature::from_kelvin(super::planetary::equilibrium_temperature_kelvin(
                        star.luminosity,
                        planet.orbital_radius.astronomical_units(),
                    ))?;
            }
        }
        Ok(())
    }

    pub(super) fn system_requires_nbody(&self, system: &StarSystem) -> bool {
        let mut companion_orbits: Vec<_> = system
            .companion_star_ids
            .iter()
            .filter_map(|id| self.stars[id].orbital_elements)
            .collect();
        companion_orbits.sort_by(|left, right| {
            left.semi_major_axis()
                .meters()
                .total_cmp(&right.semi_major_axis().meters())
        });
        let max_planet_apoapsis = system
            .planet_ids
            .iter()
            .filter_map(|id| self.planets.get(id))
            .filter(|planet| {
                !matches!(
                    planet.status,
                    PlanetStatus::Destroyed | PlanetStatus::Unbound
                )
            })
            .map(|planet| {
                planet.orbital_elements.semi_major_axis().meters()
                    * (1.0 + planet.orbital_elements.eccentricity())
            })
            .fold(0.0_f64, f64::max);
        let mut planet_orbits: Vec<_> = system
            .planet_ids
            .iter()
            .filter_map(|id| self.planets.get(id))
            .filter(|planet| {
                !matches!(
                    planet.status,
                    PlanetStatus::Destroyed | PlanetStatus::Unbound
                )
            })
            .collect();
        planet_orbits.sort_by(|left, right| {
            left.orbital_radius
                .meters()
                .total_cmp(&right.orbital_radius.meters())
        });
        if planet_orbits.windows(2).any(|pair| {
            pair[1].orbital_radius.meters() / pair[0].orbital_radius.meters() < 1.5
                || pair[0].orbital_elements.eccentricity() + pair[1].orbital_elements.eccentricity()
                    > 0.5
        }) {
            return true;
        }
        if max_planet_apoapsis > 0.0
            && companion_orbits.iter().any(|orbit| {
                orbit.semi_major_axis().meters() * (1.0 - orbit.eccentricity())
                    < 5.0 * max_planet_apoapsis
            })
        {
            return true;
        }
        companion_orbits.windows(2).any(|orbits| {
            let inner_apoapsis =
                orbits[0].semi_major_axis().meters() * (1.0 + orbits[0].eccentricity());
            let outer_periapsis =
                orbits[1].semi_major_axis().meters() * (1.0 - orbits[1].eccentricity());
            outer_periapsis < 5.0 * inner_apoapsis
        })
    }

    pub(super) fn integrate_system_from_epoch(
        &mut self,
        system_id: StarSystemId,
        elapsed: SimDuration,
    ) -> Result<bool, ObjectError> {
        const STEPS_PER_SHORTEST_ORBIT: f64 = 32.0;
        const MAX_NBODY_STEPS: usize = 4_096;
        let system = &self.star_systems[&system_id];
        if !self.system_requires_nbody(system) || elapsed == SimDuration::ZERO {
            return Ok(false);
        }
        let primary = &self.stars[&system.star_id];
        let active_planets: Vec<_> = system
            .planet_ids
            .iter()
            .filter_map(|id| self.planets.get(id))
            .filter(|planet| {
                !matches!(
                    planet.status,
                    PlanetStatus::Destroyed | PlanetStatus::Unbound
                )
            })
            .collect();
        let companions: Vec<_> = system
            .companion_star_ids
            .iter()
            .filter_map(|id| self.stars.get(id))
            .filter(|star| star.orbital_elements.is_some())
            .collect();
        let primary_mass = primary.current_mass.kilograms();
        if primary_mass <= 0.0 {
            return Ok(false);
        }
        let shortest_axis = companions
            .iter()
            .filter_map(|star| star.orbital_elements)
            .map(|orbit| orbit.semi_major_axis().meters())
            .chain(
                active_planets
                    .iter()
                    .map(|planet| planet.orbital_elements.semi_major_axis().meters()),
            )
            .fold(f64::INFINITY, f64::min);
        if !shortest_axis.is_finite() {
            return Ok(false);
        }
        let shortest_period = std::f64::consts::TAU
            * (shortest_axis.powi(3) / (GRAVITATIONAL_CONSTANT * primary_mass)).sqrt();
        let desired_step = shortest_period / STEPS_PER_SHORTEST_ORBIT;
        let elapsed_seconds = elapsed.seconds() as f64;
        let step_count = (elapsed_seconds / desired_step).ceil() as usize;
        if step_count == 0 || step_count > MAX_NBODY_STEPS {
            return Ok(false);
        }
        let step_seconds = elapsed_seconds / step_count as f64;
        let mut keys = Vec::with_capacity(1 + companions.len() + active_planets.len());
        let mut bodies = Vec::with_capacity(keys.capacity());
        keys.push(NBodyObject::Star(primary.id));
        bodies.push(GravityBody {
            mass: primary.current_mass,
            position: Position::default(),
            velocity: Velocity::default(),
        });
        for companion in companions {
            let orbit = companion
                .orbital_elements
                .expect("filtered companion orbit");
            let interior_mass: f64 = self.star_systems[&system_id]
                .companion_star_ids
                .iter()
                .filter_map(|id| self.stars.get(id))
                .filter_map(|star| star.orbital_elements.map(|candidate| (star, candidate)))
                .filter(|(_, candidate)| {
                    candidate.semi_major_axis().meters() < orbit.semi_major_axis().meters()
                })
                .map(|(star, _)| star.current_mass.kilograms())
                .sum();
            let total_mass = Mass::from_kilograms(
                primary.current_mass.kilograms()
                    + interior_mass
                    + companion.current_mass.kilograms(),
            )?;
            let (position, velocity) = propagate_kepler(
                orbit,
                total_mass,
                SimDuration::ZERO,
                Position::default(),
                Velocity::default(),
            )?;
            keys.push(NBodyObject::Star(companion.id));
            bodies.push(GravityBody {
                mass: companion.current_mass,
                position,
                velocity,
            });
        }
        for planet in active_planets {
            let (position, velocity) = propagate_kepler(
                planet.orbital_elements,
                primary.current_mass,
                SimDuration::ZERO,
                Position::default(),
                Velocity::default(),
            )?;
            keys.push(NBodyObject::Planet(planet.id));
            bodies.push(GravityBody {
                mass: planet.mass,
                position,
                velocity,
            });
        }
        for _ in 0..step_count {
            leapfrog_step(&mut bodies, step_seconds)?;
        }
        let [sx, sy, sz] = system.position.meters();
        let [svx, svy, svz] = system.velocity.meters_per_second();
        for (key, body) in keys.into_iter().zip(bodies) {
            let [x, y, z] = body.position.meters();
            let [vx, vy, vz] = body.velocity.meters_per_second();
            let position = Position::from_meters(sx + x, sy + y, sz + z)?;
            let velocity = Velocity::from_meters_per_second(svx + vx, svy + vy, svz + vz)?;
            match key {
                NBodyObject::Star(id) => {
                    let star = self.stars.get_mut(&id).expect("integrated star exists");
                    star.position = position;
                    star.velocity = velocity;
                }
                NBodyObject::Planet(id) => {
                    let planet = self.planets.get_mut(&id).expect("integrated planet exists");
                    planet.position = position;
                    planet.velocity = velocity;
                }
            }
        }
        Ok(true)
    }
}

#[derive(Debug, Clone, Copy)]
enum NBodyObject {
    Star(StarId),
    Planet(PlanetId),
}
