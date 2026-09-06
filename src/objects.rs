use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use crate::orbital::OrbitalError;
use crate::rng::StableRng;
use crate::time::{SimDuration, TimeError};
use crate::units::UnitError;

mod dynamics;
mod generation;
mod identity;
mod model;
mod persistence;
mod planetary;
mod registry;
mod stellar;
mod validation;

#[cfg(test)]
mod tests;

pub use identity::*;
pub use model::*;
pub use registry::ObjectRegistry;

const MAX_SYSTEM_SPEED_METERS_PER_SECOND: f64 = 250_000.0;
const PARALLEL_SYSTEM_THRESHOLD: u32 = 128;
const PARALLEL_BODY_THRESHOLD: usize = 256;
const SOLAR_RADIUS_METERS: f64 = 695_700_000.0;

fn sample_stellar_mass(rng: &mut StableRng, min: f64, max: f64) -> f64 {
    let split = 0.5;
    if max <= split {
        return sample_power_law(rng, min, max, 1.3);
    }
    if min >= split {
        return sample_power_law(rng, min, max, 2.3);
    }
    let low_weight = power_law_integral(min, split, 1.3);
    let high_weight = power_law_integral(split, max, 2.3);
    if rng.next_f64() < low_weight / (low_weight + high_weight) {
        sample_power_law(rng, min, split, 1.3)
    } else {
        sample_power_law(rng, split, max, 2.3)
    }
}

fn power_law_integral(min: f64, max: f64, alpha: f64) -> f64 {
    let exponent = 1.0 - alpha;
    (max.powf(exponent) - min.powf(exponent)) / exponent
}

fn sample_power_law(rng: &mut StableRng, min: f64, max: f64, alpha: f64) -> f64 {
    let exponent = 1.0 - alpha;
    let value = min.powf(exponent) + rng.next_f64() * (max.powf(exponent) - min.powf(exponent));
    value.powf(1.0 / exponent)
}

fn tick_at_or_after(threshold: SimDuration, tick_duration: SimDuration) -> Result<u64, TimeError> {
    let seconds = threshold.seconds();
    let tick_seconds = tick_duration.seconds();
    let tick = seconds
        .checked_add(tick_seconds - 1)
        .and_then(|value| value.checked_div(tick_seconds))
        .ok_or(TimeError::Overflow)?;
    u64::try_from(tick).map_err(|_| TimeError::Overflow)
}

fn random_range(rng: &mut StableRng, min: f64, max: f64) -> f64 {
    min + (max - min) * rng.next_f64()
}

fn random_signed(rng: &mut StableRng) -> f64 {
    rng.next_f64() * 2.0 - 1.0
}

fn random_unit_vector(rng: &mut StableRng) -> [f64; 3] {
    let z = random_signed(rng);
    let angle = std::f64::consts::TAU * rng.next_f64();
    let radial = (1.0 - z * z).sqrt();
    let (sin_angle, cos_angle) = angle.sin_cos();
    [radial * cos_angle, radial * sin_angle, z]
}

fn insert_unique<K, V>(
    map: &mut BTreeMap<K, V>,
    key: K,
    value: V,
    kind: &str,
) -> Result<(), ObjectError>
where
    K: Ord + Copy + fmt::Display,
{
    if map.insert(key, value).is_some() {
        Err(ObjectError(format!("duplicate {kind} ID {key}")))
    } else {
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ObjectError(String);

impl fmt::Display for ObjectError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for ObjectError {}

impl From<UnitError> for ObjectError {
    fn from(error: UnitError) -> Self {
        Self(error.to_string())
    }
}

impl From<OrbitalError> for ObjectError {
    fn from(error: OrbitalError) -> Self {
        Self(error.to_string())
    }
}
