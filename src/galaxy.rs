use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

use crate::events::{GalaxyEncounter, GalaxyEncounterKind};
use crate::id::GalaxyId;
use crate::objects::{CentralBlackHole, GalaxyCatalogEntry, GalaxyMorphology};
use crate::rng::derive_object_seed;
use crate::time::{SimDuration, TimeError};
use crate::units::{METERS_PER_PARSEC, Position, Velocity};

pub(crate) fn agn_eddington_ratio(
    galaxy_id: GalaxyId,
    encounters: &[GalaxyEncounter],
    elapsed: SimDuration,
) -> f64 {
    const STEADY_ACCRETION_RATIO: f64 = 1.0e-6;
    (STEADY_ACCRETION_RATIO
        + encounters
            .iter()
            .filter(|event| event.involves(galaxy_id) && event.occurred_at() <= elapsed)
            .map(|event| encounter_eddington_ratio(event, elapsed))
            .sum::<f64>())
    .clamp(0.0, 1.0)
}

pub(crate) fn estimate_central_black_hole(
    galaxies: &[GalaxyCatalogEntry],
    progenitors: &BTreeSet<GalaxyId>,
    encounters: &[GalaxyEncounter],
    elapsed: SimDuration,
) -> CentralBlackHole {
    let mut mass_solar: f64 = galaxies
        .iter()
        .filter(|galaxy| progenitors.contains(&galaxy.id()))
        .map(|galaxy| galaxy.black_hole_seed_mass_solar())
        .sum();
    let weighted_spin = galaxies
        .iter()
        .filter(|galaxy| progenitors.contains(&galaxy.id()))
        .map(|galaxy| {
            (
                galaxy.black_hole_seed_spin() * galaxy.black_hole_seed_mass_solar(),
                galaxy.black_hole_seed_mass_solar(),
            )
        })
        .fold((0.0, 0.0), |(spin_sum, mass_sum), (spin, mass)| {
            (spin_sum + spin, mass_sum + mass)
        });
    let integrated_eddington_years = 1.0e-6 * elapsed.seconds() as f64 / 31_557_600.0
        + encounters
            .iter()
            .filter(|event| {
                event.occurred_at() <= elapsed
                    && (progenitors.contains(&event.first())
                        || progenitors.contains(&event.second()))
            })
            .map(|event| {
                let (peak, timescale_myr) = encounter_parameters(event.kind());
                let age_seconds = elapsed.seconds() - event.occurred_at().seconds();
                let timescale_years = timescale_myr * 1.0e6;
                peak * timescale_years
                    * (1.0 - (-(age_seconds as f64) / (timescale_years * 31_557_600.0)).exp())
            })
            .sum::<f64>();
    // At 10% radiative efficiency, an Eddington-limited black hole grows on
    // an approximately 45-million-year e-folding timescale.
    mass_solar *= (integrated_eddington_years / 45.0e6).exp();
    let spin = if weighted_spin.1 > 0.0 {
        weighted_spin.0 / weighted_spin.1
    } else {
        0.5
    };
    let ratio = 1.0e-6
        + encounters
            .iter()
            .filter(|event| {
                event.occurred_at() <= elapsed
                    && (progenitors.contains(&event.first())
                        || progenitors.contains(&event.second()))
            })
            .map(|event| encounter_eddington_ratio(event, elapsed))
            .sum::<f64>()
            .clamp(0.0, 1.0);
    CentralBlackHole::from_activity(mass_solar, spin, ratio)
}

fn encounter_eddington_ratio(event: &GalaxyEncounter, elapsed: SimDuration) -> f64 {
    let (peak, timescale_myr) = encounter_parameters(event.kind());
    let age_seconds = elapsed.seconds() - event.occurred_at().seconds();
    peak * (-(age_seconds as f64) / (timescale_myr * 1.0e6 * 31_557_600.0)).exp()
}

fn encounter_parameters(kind: GalaxyEncounterKind) -> (f64, f64) {
    match kind {
        GalaxyEncounterKind::ClosePass => (0.03, 50.0),
        GalaxyEncounterKind::Merger => (0.5, 100.0),
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct GalaxyMotion {
    pub(crate) position_at_epoch: Position,
    pub(crate) velocity: Velocity,
    pub(crate) epoch: SimDuration,
    pub(crate) morphology: GalaxyMorphology,
    pub(crate) radius_parsecs: f64,
    pub(crate) merged_into: Option<GalaxyId>,
}

impl GalaxyMotion {
    pub(crate) fn from_catalog(galaxy: GalaxyCatalogEntry, epoch: SimDuration) -> Self {
        Self {
            position_at_epoch: galaxy.position(),
            velocity: galaxy.velocity(),
            epoch,
            morphology: galaxy.morphology(),
            radius_parsecs: galaxy.radius_parsecs(),
            merged_into: None,
        }
    }

    pub(crate) fn position_at(self, time: SimDuration) -> Result<Position, TimeError> {
        let seconds = time.seconds().saturating_sub(self.epoch.seconds()) as f64;
        let start = self.position_at_epoch.meters();
        let velocity = self.velocity.meters_per_second();
        Position::from_meters(
            start[0] + velocity[0] * seconds,
            start[1] + velocity[1] * seconds,
            start[2] + velocity[2] * seconds,
        )
        .map_err(|_| TimeError::Overflow)
    }

    pub(crate) fn snapshot(
        self,
        catalog: GalaxyCatalogEntry,
        time: SimDuration,
    ) -> Result<GalaxyCatalogEntry, TimeError> {
        Ok(catalog.with_dynamics(
            self.position_at(time)?,
            self.velocity,
            self.morphology,
            self.radius_parsecs,
            self.merged_into,
        ))
    }
}

pub(crate) fn encounter_id(
    seed: u64,
    first: GalaxyId,
    second: GalaxyId,
    kind: GalaxyEncounterKind,
    seconds: u128,
) -> u64 {
    let namespace = format!("galaxy-encounter-v1-{}-{seconds}", kind.code());
    derive_object_seed(seed ^ first.value(), second.value(), &namespace)
}

pub(crate) fn advance_encounters(
    motions: &mut BTreeMap<GalaxyId, GalaxyMotion>,
    encounters: &mut Vec<GalaxyEncounter>,
    seed: u64,
    primary_id: GalaxyId,
    from: SimDuration,
    to: SimDuration,
    tick_duration: SimDuration,
) -> Result<(), TimeError> {
    let mut seen: BTreeSet<_> = encounters
        .iter()
        .map(|event| (event.first(), event.second(), event.kind()))
        .collect();
    let ids: Vec<_> = motions
        .iter()
        .filter_map(|(id, motion)| motion.merged_into.is_none().then_some(*id))
        .collect();
    let mut queue = BinaryHeap::new();
    for (index, first) in ids.iter().enumerate() {
        for second in &ids[index + 1..] {
            enqueue_pair(&mut queue, motions, &seen, *first, *second, from, to)?;
        }
    }
    let mut cursor = from;
    while let Some(Reverse((seconds, kind, first_id, second_id))) = queue.pop() {
        if seen.contains(&(first_id, second_id, kind))
            || motions[&first_id].merged_into.is_some()
            || motions[&second_id].merged_into.is_some()
        {
            continue;
        }
        let first = motions[&first_id];
        let second = motions[&second_id];
        let threshold = encounter_threshold(first, second, kind);
        if entry_time(first, second, threshold, cursor)? != Some(seconds) {
            continue;
        }
        let at = SimDuration::from_seconds(seconds);
        let tick = u64::try_from(
            seconds
                .checked_add(tick_duration.seconds() - 1)
                .ok_or(TimeError::Overflow)?
                / tick_duration.seconds(),
        )
        .map_err(|_| TimeError::Overflow)?;
        let event = GalaxyEncounter::new(
            encounter_id(seed, first_id, second_id, kind, seconds),
            at,
            tick,
            first_id,
            second_id,
            kind,
        );
        encounters.push(event);
        seen.insert((first_id, second_id, kind));
        if kind == GalaxyEncounterKind::Merger {
            let survivor_id = if first_id == primary_id || second_id != primary_id {
                first_id
            } else {
                second_id
            };
            merge(motions, primary_id, first_id, second_id, at)?;
            let others: Vec<_> = motions
                .iter()
                .filter_map(|(id, motion)| {
                    (motion.merged_into.is_none() && *id != survivor_id).then_some(*id)
                })
                .collect();
            for other_id in others {
                let (a, b) = if survivor_id < other_id {
                    (survivor_id, other_id)
                } else {
                    (other_id, survivor_id)
                };
                enqueue_pair(&mut queue, motions, &seen, a, b, at, to)?;
            }
        }
        cursor = at;
    }
    for motion in motions.values() {
        motion.position_at(to)?;
    }
    encounters.sort_by_key(|event| (event.occurred_at(), event.id()));
    Ok(())
}

fn encounter_threshold(
    first: GalaxyMotion,
    second: GalaxyMotion,
    kind: GalaxyEncounterKind,
) -> f64 {
    let radius_multiplier = match kind {
        GalaxyEncounterKind::ClosePass => 3.0,
        GalaxyEncounterKind::Merger => 0.5,
    };
    radius_multiplier * (first.radius_parsecs + second.radius_parsecs) * METERS_PER_PARSEC
}

fn enqueue_pair(
    queue: &mut BinaryHeap<Reverse<(u128, GalaxyEncounterKind, GalaxyId, GalaxyId)>>,
    motions: &BTreeMap<GalaxyId, GalaxyMotion>,
    seen: &BTreeSet<(GalaxyId, GalaxyId, GalaxyEncounterKind)>,
    first_id: GalaxyId,
    second_id: GalaxyId,
    from: SimDuration,
    to: SimDuration,
) -> Result<(), TimeError> {
    let first = motions[&first_id];
    let second = motions[&second_id];
    for kind in [GalaxyEncounterKind::ClosePass, GalaxyEncounterKind::Merger] {
        if seen.contains(&(first_id, second_id, kind)) {
            continue;
        }
        let threshold = encounter_threshold(first, second, kind);
        if let Some(seconds) = entry_time(first, second, threshold, from)?
            && seconds <= to.seconds()
        {
            queue.push(Reverse((seconds, kind, first_id, second_id)));
        }
    }
    Ok(())
}

fn entry_time(
    first: GalaxyMotion,
    second: GalaxyMotion,
    threshold: f64,
    from: SimDuration,
) -> Result<Option<u128>, TimeError> {
    let position_at = |time: SimDuration| -> Result<[f64; 3], TimeError> {
        let a = first.position_at(time)?.meters();
        let b = second.position_at(time)?.meters();
        Ok([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
    };
    let distance_squared = |r: [f64; 3]| r.iter().map(|value| value * value).sum::<f64>();
    if distance_squared(position_at(from)?) <= threshold * threshold {
        return Ok(Some(from.seconds()));
    }
    let base = first.epoch.max(second.epoch);
    let r = position_at(base)?;
    let a = first.velocity.meters_per_second();
    let b = second.velocity.meters_per_second();
    let v = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let speed_squared = v.iter().map(|value| value * value).sum::<f64>();
    if speed_squared == 0.0 {
        return Ok(None);
    }
    let radial = r
        .iter()
        .zip(v)
        .map(|(position, speed)| position * speed)
        .sum::<f64>();
    let discriminant =
        radial * radial - speed_squared * (distance_squared(r) - threshold * threshold);
    if discriminant < 0.0 {
        return Ok(None);
    }
    let offset = (-radial - discriminant.sqrt()) / speed_squared;
    if !offset.is_finite() || offset < 0.0 || offset >= u128::MAX as f64 {
        return Ok(None);
    }
    let seconds = base.seconds().checked_add(offset.ceil() as u128);
    Ok(seconds.filter(|seconds| *seconds >= from.seconds()))
}

fn merge(
    motions: &mut BTreeMap<GalaxyId, GalaxyMotion>,
    primary_id: GalaxyId,
    first_id: GalaxyId,
    second_id: GalaxyId,
    at: SimDuration,
) -> Result<(), TimeError> {
    let (survivor_id, absorbed_id) = if first_id == primary_id {
        (first_id, second_id)
    } else if second_id == primary_id {
        (second_id, first_id)
    } else {
        (first_id, second_id)
    };
    let survivor = motions[&survivor_id];
    let absorbed = motions[&absorbed_id];
    let weights = [
        survivor.radius_parsecs.powi(3),
        absorbed.radius_parsecs.powi(3),
    ];
    let total = weights[0] + weights[1];
    let p1 = survivor.position_at(at)?.meters();
    let p2 = absorbed.position_at(at)?.meters();
    let v1 = survivor.velocity.meters_per_second();
    let v2 = absorbed.velocity.meters_per_second();
    let weighted = |left: [f64; 3], right: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|axis| (weights[0] * left[axis] + weights[1] * right[axis]) / total)
    };
    let position = weighted(p1, p2);
    let velocity = weighted(v1, v2);
    *motions.get_mut(&survivor_id).expect("survivor exists") = GalaxyMotion {
        position_at_epoch: Position::from_meters(position[0], position[1], position[2])
            .map_err(|_| TimeError::Overflow)?,
        velocity: Velocity::from_meters_per_second(velocity[0], velocity[1], velocity[2])
            .map_err(|_| TimeError::Overflow)?,
        epoch: at,
        morphology: GalaxyMorphology::Elliptical,
        radius_parsecs: total.cbrt(),
        merged_into: None,
    };
    *motions
        .get_mut(&absorbed_id)
        .expect("absorbed galaxy exists") = GalaxyMotion {
        position_at_epoch: absorbed.position_at(at)?,
        velocity: Velocity::default(),
        epoch: at,
        merged_into: Some(survivor_id),
        ..absorbed
    };
    Ok(())
}
