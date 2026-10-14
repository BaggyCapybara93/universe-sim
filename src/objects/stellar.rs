use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};

use rayon::prelude::*;

use crate::chemistry::GasComposition;
use crate::id::StarId;
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration, TimeError};
use crate::units::{Luminosity, Mass, Temperature, Velocity};

use super::{
    ObjectError, ObjectRegistry, PARALLEL_BODY_THRESHOLD, Star, StellarInteractionKind,
    StellarStage, StellarTransition, random_unit_vector, tick_at_or_after,
};

pub(super) fn stellar_schedule(
    initial_mass: Mass,
    metallicity: f64,
) -> Result<Vec<(SimDuration, StellarStage)>, TimeError> {
    let solar_masses = initial_mass.solar_masses();
    let metallicity_factor = (metallicity / 0.014).clamp(0.05, 3.0).powf(0.12);
    let main_sequence_years = (10_000_000_000.0 * solar_masses.powf(-2.5) * metallicity_factor)
        .clamp(2_500_000.0, 200_000_000_000.0);
    let pre_main_sequence_years =
        (30_000_000.0 * solar_masses.powf(-2.5)).clamp(50_000.0, 2_000_000_000.0);
    let protostar_years = (pre_main_sequence_years * 0.15).clamp(10_000.0, 500_000.0);
    let main_sequence_end = pre_main_sequence_years + main_sequence_years;
    let mut schedule = vec![
        (
            years_to_duration(protostar_years)?,
            StellarStage::PreMainSequence,
        ),
        (
            years_to_duration(pre_main_sequence_years)?,
            StellarStage::MainSequence,
        ),
    ];
    match stellar_fate(initial_mass, metallicity) {
        StellarFate::WhiteDwarf => {
            schedule.push((
                years_to_duration(main_sequence_end)?,
                StellarStage::Subgiant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.05 * main_sequence_years)?,
                StellarStage::Giant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.10 * main_sequence_years)?,
                StellarStage::CoreHeliumBurning,
            ));
            let agb_start = main_sequence_end + 0.11 * main_sequence_years;
            let agb_years = (1_000_000.0 * solar_masses.powf(-1.5)).clamp(30_000.0, 2_000_000.0);
            schedule.push((
                years_to_duration(agb_start)?,
                StellarStage::AsymptoticGiantBranch,
            ));
            schedule.push((
                years_to_duration(agb_start + agb_years)?,
                StellarStage::PlanetaryNebula,
            ));
            schedule.push((
                years_to_duration(agb_start + agb_years + 30_000.0)?,
                StellarStage::WhiteDwarf,
            ));
        }
        StellarFate::NeutronStar | StellarFate::BlackHole => {
            schedule_massive_post_main_sequence(
                &mut schedule,
                main_sequence_end,
                main_sequence_years,
                solar_masses,
                metallicity,
            )?;
            let collapse = main_sequence_end + 0.05 * main_sequence_years;
            schedule.push((
                years_to_duration(collapse)?,
                StellarStage::CoreCollapseSupernova,
            ));
            let remnant = if matches!(
                stellar_fate(initial_mass, metallicity),
                StellarFate::NeutronStar
            ) {
                StellarStage::NeutronStar
            } else {
                StellarStage::BlackHole
            };
            schedule.push((years_to_duration(collapse + 1.0 / 365.25)?, remnant));
        }
        StellarFate::DirectCollapse => {
            schedule_massive_post_main_sequence(
                &mut schedule,
                main_sequence_end,
                main_sequence_years,
                solar_masses,
                metallicity,
            )?;
            schedule.push((
                years_to_duration(main_sequence_end + 0.05 * main_sequence_years)?,
                StellarStage::BlackHole,
            ));
        }
        StellarFate::PulsationalPairInstability => {
            schedule.push((
                years_to_duration(main_sequence_end)?,
                StellarStage::Supergiant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.03 * main_sequence_years)?,
                StellarStage::PulsationalPairInstability,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.08 * main_sequence_years)?,
                StellarStage::BlackHole,
            ));
        }
        StellarFate::PairInstability => {
            schedule.push((
                years_to_duration(main_sequence_end)?,
                StellarStage::Supergiant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.05 * main_sequence_years)?,
                StellarStage::Disrupted,
            ));
        }
    }
    Ok(schedule)
}

fn schedule_massive_post_main_sequence(
    schedule: &mut Vec<(SimDuration, StellarStage)>,
    main_sequence_end: f64,
    main_sequence_years: f64,
    solar_masses: f64,
    metallicity: f64,
) -> Result<(), TimeError> {
    schedule.push((
        years_to_duration(main_sequence_end)?,
        StellarStage::Supergiant,
    ));
    if solar_masses >= 25.0 && metallicity >= 0.008 {
        schedule.push((
            years_to_duration(main_sequence_end + 0.025 * main_sequence_years)?,
            StellarStage::WolfRayet,
        ));
    }
    Ok(())
}

fn years_to_duration(years: f64) -> Result<SimDuration, TimeError> {
    let seconds = years * SECONDS_PER_YEAR as f64;
    if !seconds.is_finite() || seconds < 0.0 || seconds > u128::MAX as f64 {
        return Err(TimeError::Overflow);
    }
    Ok(SimDuration::from_seconds(seconds.round() as u128))
}

pub(super) fn stellar_stage_at(
    initial_mass: Mass,
    metallicity: f64,
    age: SimDuration,
) -> Result<StellarStage, ObjectError> {
    let schedule = stellar_schedule(initial_mass, metallicity)
        .map_err(|error| ObjectError(error.to_string()))?;
    let mut stage = StellarStage::Protostar;
    for (threshold, scheduled_stage) in schedule {
        if age < threshold {
            break;
        }
        stage = scheduled_stage;
    }
    Ok(stage)
}

pub(super) fn stellar_properties(
    initial_mass: Mass,
    metallicity: f64,
    stage: StellarStage,
) -> Result<(Mass, Luminosity, Temperature), ObjectError> {
    let solar_masses = initial_mass.solar_masses();
    let retained = wind_retained_fraction(solar_masses, metallicity);
    let pre_collapse_mass = solar_masses * retained;
    let main_sequence_luminosity = if solar_masses <= 20.0 {
        solar_masses.powf(3.5)
    } else {
        20.0_f64.powf(3.5) * (solar_masses / 20.0).powf(1.5)
    };
    let (current_mass_solar, luminosity_solar, temperature_kelvin) = match stage {
        StellarStage::Protostar => (
            solar_masses,
            (2.0 * solar_masses.powf(2.0)).max(0.001),
            (3_000.0 * solar_masses.powf(0.08)).clamp(2_500.0, 5_000.0),
        ),
        StellarStage::PreMainSequence => (
            solar_masses * 0.995,
            (1.5 * solar_masses.powf(2.5)).max(0.0005),
            (4_200.0 * solar_masses.powf(0.15)).clamp(2_700.0, 15_000.0),
        ),
        StellarStage::MainSequence => (
            solar_masses * (1.0 - 0.5 * (1.0 - retained)),
            main_sequence_luminosity.max(0.0001),
            (5_772.0 * solar_masses.powf(0.505)).min(60_000.0),
        ),
        StellarStage::Subgiant => (
            solar_masses * 0.97,
            (3.0 * main_sequence_luminosity).max(0.001),
            (5_200.0 * solar_masses.powf(0.08)).clamp(3_800.0, 8_000.0),
        ),
        StellarStage::Giant => (
            solar_masses * 0.85,
            1_000.0 * solar_masses.powf(2.0),
            3_800.0,
        ),
        StellarStage::CoreHeliumBurning => {
            (solar_masses * 0.80, 50.0 * solar_masses.powf(2.0), 5_000.0)
        }
        StellarStage::AsymptoticGiantBranch => (
            solar_masses * 0.65,
            10_000.0 * solar_masses.powf(1.2),
            3_200.0,
        ),
        StellarStage::PlanetaryNebula => (
            (0.109 * solar_masses + 0.394).clamp(0.45, 1.38),
            1_000.0,
            100_000.0,
        ),
        StellarStage::Supergiant => (
            pre_collapse_mass,
            30_000.0 * solar_masses.powf(1.5),
            4_500.0,
        ),
        StellarStage::WolfRayet => (
            pre_collapse_mass,
            50_000.0 * solar_masses.powf(1.3),
            80_000.0,
        ),
        StellarStage::PulsationalPairInstability => {
            (pre_collapse_mass.min(60.0), 500_000.0, 12_000.0)
        }
        StellarStage::CoreCollapseSupernova => (pre_collapse_mass, 1_000_000_000.0, 10_000.0),
        StellarStage::WhiteDwarf => (
            (0.109 * solar_masses + 0.394).clamp(0.45, 1.38),
            0.01,
            12_000.0,
        ),
        StellarStage::NeutronStar => (
            (1.25 + 0.05 * (solar_masses - 8.0)).clamp(1.2, 2.2),
            0.0001,
            600_000.0,
        ),
        StellarStage::BlackHole => {
            let remnant = match stellar_fate(initial_mass, metallicity) {
                StellarFate::PulsationalPairInstability => pre_collapse_mass.min(50.0) * 0.9,
                StellarFate::DirectCollapse => pre_collapse_mass * 0.9,
                _ => (3.0 + 0.25 * solar_masses).min(pre_collapse_mass).max(3.0),
            };
            (remnant, 0.0, 0.0)
        }
        StellarStage::Disrupted => (0.0, 0.0, 0.0),
    };
    Ok((
        Mass::from_solar_masses(current_mass_solar)?,
        Luminosity::from_solar_luminosities(luminosity_solar)?,
        Temperature::from_kelvin(temperature_kelvin)?,
    ))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum StellarFate {
    WhiteDwarf,
    NeutronStar,
    BlackHole,
    PulsationalPairInstability,
    PairInstability,
    DirectCollapse,
}

fn wind_retained_fraction(initial_mass_solar: f64, metallicity: f64) -> f64 {
    if initial_mass_solar < 25.0 {
        return 1.0;
    }
    let metallicity_ratio = (metallicity / 0.014).clamp(0.001, 3.0);
    let wind_loss = (0.03
        + 0.45 * metallicity_ratio.powf(0.7) * (initial_mass_solar / 100.0).sqrt())
    .clamp(0.03, 0.80);
    1.0 - wind_loss
}

pub(super) fn stellar_fate(initial_mass: Mass, metallicity: f64) -> StellarFate {
    let solar_masses = initial_mass.solar_masses();
    if solar_masses < 8.0 {
        return StellarFate::WhiteDwarf;
    }
    if solar_masses < 25.0 {
        return StellarFate::NeutronStar;
    }
    let helium_core = 0.55 * solar_masses * wind_retained_fraction(solar_masses, metallicity);
    if helium_core < 32.0 {
        StellarFate::BlackHole
    } else if helium_core < 64.0 {
        StellarFate::PulsationalPairInstability
    } else if helium_core < 135.0 {
        StellarFate::PairInstability
    } else {
        StellarFate::DirectCollapse
    }
}

pub(super) fn natal_kick_speed(star: &Star, stage: StellarStage, tick: u64) -> u32 {
    let mut rng = StableRng::new(derive_object_seed(
        star.id.value(),
        tick,
        "natal-kick-speed-v1",
    ));
    match stage {
        StellarStage::NeutronStar => 100_000 + rng.next_u64_below(400_001).unwrap_or(0) as u32,
        StellarStage::BlackHole => match stellar_fate(star.initial_mass, star.metallicity) {
            StellarFate::DirectCollapse => 0,
            StellarFate::PulsationalPairInstability => {
                10_000 + rng.next_u64_below(40_001).unwrap_or(0) as u32
            }
            _ => 20_000 + rng.next_u64_below(80_001).unwrap_or(0) as u32,
        },
        _ => 0,
    }
}

pub(super) fn stellar_radius_solar(mass: Mass, stage: StellarStage) -> f64 {
    let solar_masses = mass.solar_masses().max(0.01);
    match stage {
        StellarStage::Protostar => 3.0 * solar_masses.powf(0.5),
        StellarStage::PreMainSequence => 2.0 * solar_masses.powf(0.6),
        StellarStage::MainSequence => {
            solar_masses.powf(if solar_masses < 1.0 { 0.8 } else { 0.57 })
        }
        StellarStage::Subgiant => 3.0 * solar_masses.powf(0.6),
        StellarStage::Giant => 100.0 * solar_masses.powf(0.5),
        StellarStage::CoreHeliumBurning => 10.0 * solar_masses.powf(0.5),
        StellarStage::AsymptoticGiantBranch => 200.0 * solar_masses.powf(0.5),
        StellarStage::PlanetaryNebula => 0.02,
        StellarStage::Supergiant => 500.0 * (solar_masses / 10.0).sqrt(),
        StellarStage::WolfRayet => 3.0 * solar_masses.powf(0.4),
        StellarStage::PulsationalPairInstability => 100.0,
        StellarStage::CoreCollapseSupernova => 500.0,
        StellarStage::WhiteDwarf => 0.012,
        StellarStage::NeutronStar => 0.000_02,
        StellarStage::BlackHole | StellarStage::Disrupted => 0.0,
    }
}

pub(super) fn roche_lobe_radius(separation: f64, mass_ratio: f64) -> f64 {
    let q_one_third = mass_ratio.max(1e-12).cbrt();
    let q_two_thirds = q_one_third * q_one_third;
    separation * 0.49 * q_two_thirds / (0.6 * q_two_thirds + (1.0 + q_one_third).ln())
}

// Use physical seconds for ordering even when several transitions share a tick.
// Birth time and initial mass invalidate queued entries after rejuvenating mergers.
type ScheduledTransition = (u128, StarId, StellarStage, u128, u64);

fn scheduled_transitions(
    star: &Star,
    until: SimDuration,
) -> Result<Vec<ScheduledTransition>, TimeError> {
    let target_age = star
        .epoch_age
        .checked_add(SimDuration::from_seconds(
            until.seconds().saturating_sub(star.formed_at.seconds()),
        ))
        .ok_or(TimeError::Overflow)?;
    let mut items = Vec::new();
    for (threshold, stage) in stellar_schedule(star.initial_mass, star.metallicity)? {
        if star.age < threshold && threshold <= target_age && stage > star.stage {
            let seconds = star
                .formed_at
                .seconds()
                .checked_add(threshold.seconds().saturating_sub(star.epoch_age.seconds()))
                .ok_or(TimeError::Overflow)?;
            items.push((
                seconds,
                star.id,
                stage,
                star.formed_at.seconds(),
                star.initial_mass.kilograms().to_bits(),
            ));
        }
    }
    Ok(items)
}

fn ejecta_composition(
    birth: GasComposition,
    from: StellarStage,
    to: StellarStage,
) -> GasComposition {
    match to {
        StellarStage::NeutronStar | StellarStage::Disrupted => birth.processed(0.35, 0.20),
        StellarStage::BlackHole if from == StellarStage::CoreCollapseSupernova => {
            birth.processed(0.35, 0.20)
        }
        StellarStage::Supergiant
        | StellarStage::WolfRayet
        | StellarStage::PulsationalPairInstability => birth.processed(0.20, 0.08),
        StellarStage::Giant
        | StellarStage::CoreHeliumBurning
        | StellarStage::AsymptoticGiantBranch
        | StellarStage::PlanetaryNebula
        | StellarStage::WhiteDwarf => birth.processed(0.10, 0.05),
        _ => birth,
    }
}

impl ObjectRegistry {
    pub(super) fn next_stellar_transition_after(
        &self,
        elapsed: SimDuration,
    ) -> Result<Option<SimDuration>, TimeError> {
        let mut next = None;
        for star in self.stars.values() {
            for (threshold, stage) in stellar_schedule(star.initial_mass, star.metallicity)? {
                if stage <= star.stage {
                    continue;
                }
                let seconds = star
                    .formed_at
                    .seconds()
                    .checked_add(threshold.seconds().saturating_sub(star.epoch_age.seconds()))
                    .ok_or(TimeError::Overflow)?;
                if seconds > elapsed.seconds() {
                    next = Some(next.map_or(seconds, |old: u128| old.min(seconds)));
                }
            }
        }
        Ok(next.map(SimDuration::from_seconds))
    }

    pub(crate) fn advance_stars_to(
        &mut self,
        elapsed: SimDuration,
        tick_duration: SimDuration,
    ) -> Result<Vec<StellarTransition>, TimeError> {
        let mut scheduled = BinaryHeap::new();
        let stars: Vec<_> = self.stars.values().collect();
        let schedule_star = |star: &&Star| scheduled_transitions(star, elapsed);
        let schedules: Result<Vec<_>, TimeError> = if stars.len() >= PARALLEL_BODY_THRESHOLD {
            stars.par_iter().map(schedule_star).collect()
        } else {
            stars.iter().map(schedule_star).collect()
        };
        for item in schedules?.into_iter().flatten() {
            scheduled.push(Reverse(item));
        }

        let mut transitions = Vec::new();
        while let Some(Reverse((seconds, star_id, stage, birth_seconds, initial_mass_bits))) =
            scheduled.pop()
        {
            let Some(star) = self.stars.get_mut(&star_id) else {
                continue;
            };
            if star.formed_at.seconds() != birth_seconds
                || star.initial_mass.kilograms().to_bits() != initial_mass_bits
                || star.stage >= stage
            {
                continue;
            }
            let event_time = SimDuration::from_seconds(seconds);
            let tick = tick_at_or_after(event_time, tick_duration)?;
            let from = star.stage;
            let previous_mass = star.current_mass;
            star.stage = stage;
            let (prescribed_mass, luminosity, temperature) =
                stellar_properties(star.initial_mass, star.metallicity, stage)
                    .map_err(|_| TimeError::Overflow)?;
            // Stage prescriptions must never recreate previously lost/transferred mass.
            let mass =
                Mass::from_kilograms(prescribed_mass.kilograms().min(previous_mass.kilograms()))
                    .map_err(|_| TimeError::Overflow)?;
            star.current_mass = mass;
            star.luminosity = luminosity;
            star.temperature = temperature;
            let composition = ejecta_composition(star.birth_composition, from, stage);
            let natal_kick_meters_per_second = natal_kick_speed(star, stage, tick);
            if natal_kick_meters_per_second > 0 {
                let mut rng = StableRng::new(derive_object_seed(
                    star.id.value(),
                    tick,
                    "natal-kick-direction-v1",
                ));
                let direction = random_unit_vector(&mut rng);
                let [vx, vy, vz] = star.velocity.meters_per_second();
                let speed = f64::from(natal_kick_meters_per_second);
                star.velocity = Velocity::from_meters_per_second(
                    vx + direction[0] * speed,
                    vy + direction[1] * speed,
                    vz + direction[2] * speed,
                )
                .map_err(|_| TimeError::Overflow)?;
            }
            self.return_gas(
                star_id,
                previous_mass.solar_masses() - mass.solar_masses(),
                composition,
            )
            .map_err(|_| TimeError::Overflow)?;
            let collision_impacts = self.apply_direct_planet_collisions(star_id, stage);
            let impacts = self.apply_stellar_impacts(star_id, previous_mass, mass, stage);
            transitions.push(StellarTransition {
                tick,
                star_id,
                from,
                to: stage,
                impacts,
                natal_kick_meters_per_second,
                interaction: None,
                related_star_id: None,
            });
            if !collision_impacts.is_empty() {
                transitions.push(StellarTransition {
                    tick,
                    star_id,
                    from,
                    to: stage,
                    impacts: collision_impacts,
                    natal_kick_meters_per_second: 0,
                    interaction: Some(StellarInteractionKind::PlanetaryCollision),
                    related_star_id: None,
                });
            }
            let interactions =
                self.resolve_binary_interactions(&BTreeMap::from([(star_id, tick)]), event_time)?;
            for interaction in &interactions {
                if interaction.interaction == Some(StellarInteractionKind::StellarMerger) {
                    let merged = &self.stars[&interaction.star_id];
                    for item in scheduled_transitions(merged, elapsed)? {
                        scheduled.push(Reverse(item));
                    }
                }
            }
            transitions.extend(interactions);
        }
        for star in self.stars.values_mut() {
            star.age = star
                .epoch_age
                .checked_add(SimDuration::from_seconds(
                    elapsed.seconds().saturating_sub(star.formed_at.seconds()),
                ))
                .ok_or(TimeError::Overflow)?;
        }
        self.refresh_planet_temperatures()?;
        Ok(transitions)
    }
}
