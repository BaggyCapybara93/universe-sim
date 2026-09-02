use std::error::Error;
use std::f64::consts::TAU;
use std::fmt;

use crate::time::SimDuration;
use crate::units::{Angle, Length, Mass, Position, UnitError, Velocity};

pub const GRAVITATIONAL_CONSTANT: f64 = 6.674_30e-11;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitalElements {
    semi_major_axis: Length,
    eccentricity: f64,
    inclination: Angle,
    longitude_of_ascending_node: Angle,
    argument_of_periapsis: Angle,
    mean_anomaly_at_epoch: Angle,
    epoch: SimDuration,
}

impl OrbitalElements {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        semi_major_axis: Length,
        eccentricity: f64,
        inclination: Angle,
        longitude_of_ascending_node: Angle,
        argument_of_periapsis: Angle,
        mean_anomaly_at_epoch: Angle,
        epoch: SimDuration,
    ) -> Result<Self, OrbitalError> {
        if semi_major_axis.meters() <= 0.0 {
            return Err(OrbitalError::InvalidElements(
                "semi-major axis must be greater than zero".to_owned(),
            ));
        }
        if !eccentricity.is_finite() || !(0.0..1.0).contains(&eccentricity) {
            return Err(OrbitalError::InvalidElements(
                "eccentricity must be finite and in [0, 1)".to_owned(),
            ));
        }
        Ok(Self {
            semi_major_axis,
            eccentricity,
            inclination,
            longitude_of_ascending_node,
            argument_of_periapsis,
            mean_anomaly_at_epoch,
            epoch,
        })
    }

    pub const fn semi_major_axis(self) -> Length {
        self.semi_major_axis
    }

    pub const fn eccentricity(self) -> f64 {
        self.eccentricity
    }

    pub const fn inclination(self) -> Angle {
        self.inclination
    }

    pub const fn longitude_of_ascending_node(self) -> Angle {
        self.longitude_of_ascending_node
    }

    pub const fn argument_of_periapsis(self) -> Angle {
        self.argument_of_periapsis
    }

    pub const fn mean_anomaly_at_epoch(self) -> Angle {
        self.mean_anomaly_at_epoch
    }

    pub const fn epoch(self) -> SimDuration {
        self.epoch
    }
}

pub fn propagate_kepler(
    elements: OrbitalElements,
    primary_mass: Mass,
    elapsed: SimDuration,
    origin: Position,
    origin_velocity: Velocity,
) -> Result<(Position, Velocity), OrbitalError> {
    let mu = GRAVITATIONAL_CONSTANT * primary_mass.kilograms();
    if mu <= 0.0 || !mu.is_finite() {
        return Err(OrbitalError::InvalidPrimaryMass);
    }
    let time_seconds = if elapsed >= elements.epoch {
        (elapsed.seconds() - elements.epoch.seconds()) as f64
    } else {
        -((elements.epoch.seconds() - elapsed.seconds()) as f64)
    };
    let a = elements.semi_major_axis.meters();
    let eccentricity = elements.eccentricity;
    let mean_motion = (mu / a.powi(3)).sqrt();
    let mean_anomaly =
        (elements.mean_anomaly_at_epoch.radians() + mean_motion * time_seconds).rem_euclid(TAU);
    let eccentric_anomaly = solve_kepler(mean_anomaly, eccentricity)?;

    let cos_e = eccentric_anomaly.cos();
    let sin_e = eccentric_anomaly.sin();
    let one_minus_e_cos = 1.0 - eccentricity * cos_e;
    let sqrt_one_minus_e2 = (1.0 - eccentricity * eccentricity).sqrt();
    let plane_position = [
        a * (cos_e - eccentricity),
        a * sqrt_one_minus_e2 * sin_e,
        0.0,
    ];
    let velocity_scale = mean_motion * a / one_minus_e_cos;
    let plane_velocity = [
        -velocity_scale * sin_e,
        velocity_scale * sqrt_one_minus_e2 * cos_e,
        0.0,
    ];
    let position_delta = rotate_from_orbital_plane(plane_position, elements);
    let velocity_delta = rotate_from_orbital_plane(plane_velocity, elements);
    let [ox, oy, oz] = origin.meters();
    let [ovx, ovy, ovz] = origin_velocity.meters_per_second();
    Ok((
        Position::from_meters(
            ox + position_delta[0],
            oy + position_delta[1],
            oz + position_delta[2],
        )?,
        Velocity::from_meters_per_second(
            ovx + velocity_delta[0],
            ovy + velocity_delta[1],
            ovz + velocity_delta[2],
        )?,
    ))
}

fn solve_kepler(mean_anomaly: f64, eccentricity: f64) -> Result<f64, OrbitalError> {
    let mut eccentric_anomaly = if eccentricity < 0.8 {
        mean_anomaly
    } else {
        std::f64::consts::PI
    };
    for _ in 0..20 {
        let residual = eccentric_anomaly - eccentricity * eccentric_anomaly.sin() - mean_anomaly;
        let derivative = 1.0 - eccentricity * eccentric_anomaly.cos();
        let correction = residual / derivative;
        eccentric_anomaly -= correction;
        if correction.abs() < 1e-13 {
            return Ok(eccentric_anomaly);
        }
    }
    Err(OrbitalError::KeplerDidNotConverge)
}

fn rotate_from_orbital_plane(vector: [f64; 3], elements: OrbitalElements) -> [f64; 3] {
    let node = elements.longitude_of_ascending_node.radians();
    let inclination = elements.inclination.radians();
    let periapsis = elements.argument_of_periapsis.radians();
    let (sin_node, cos_node) = node.sin_cos();
    let (sin_inclination, cos_inclination) = inclination.sin_cos();
    let (sin_periapsis, cos_periapsis) = periapsis.sin_cos();
    let [x, y, _] = vector;
    [
        (cos_node * cos_periapsis - sin_node * sin_periapsis * cos_inclination) * x
            + (-cos_node * sin_periapsis - sin_node * cos_periapsis * cos_inclination) * y,
        (sin_node * cos_periapsis + cos_node * sin_periapsis * cos_inclination) * x
            + (-sin_node * sin_periapsis + cos_node * cos_periapsis * cos_inclination) * y,
        sin_periapsis * sin_inclination * x + cos_periapsis * sin_inclination * y,
    ]
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GravityBody {
    pub mass: Mass,
    pub position: Position,
    pub velocity: Velocity,
}

pub fn leapfrog_step(bodies: &mut [GravityBody], step_seconds: f64) -> Result<(), OrbitalError> {
    if !step_seconds.is_finite() || step_seconds <= 0.0 {
        return Err(OrbitalError::InvalidStep);
    }
    let initial_accelerations = accelerations(bodies)?;
    for (body, acceleration) in bodies.iter_mut().zip(&initial_accelerations) {
        let [vx, vy, vz] = body.velocity.meters_per_second();
        let half_velocity = [
            vx + 0.5 * step_seconds * acceleration[0],
            vy + 0.5 * step_seconds * acceleration[1],
            vz + 0.5 * step_seconds * acceleration[2],
        ];
        let [x, y, z] = body.position.meters();
        body.position = Position::from_meters(
            x + step_seconds * half_velocity[0],
            y + step_seconds * half_velocity[1],
            z + step_seconds * half_velocity[2],
        )?;
        body.velocity =
            Velocity::from_meters_per_second(half_velocity[0], half_velocity[1], half_velocity[2])?;
    }
    let final_accelerations = accelerations(bodies)?;
    for (body, acceleration) in bodies.iter_mut().zip(&final_accelerations) {
        let [vx, vy, vz] = body.velocity.meters_per_second();
        body.velocity = Velocity::from_meters_per_second(
            vx + 0.5 * step_seconds * acceleration[0],
            vy + 0.5 * step_seconds * acceleration[1],
            vz + 0.5 * step_seconds * acceleration[2],
        )?;
    }
    Ok(())
}

fn accelerations(bodies: &[GravityBody]) -> Result<Vec<[f64; 3]>, OrbitalError> {
    let mut values = vec![[0.0; 3]; bodies.len()];
    for first in 0..bodies.len() {
        for second in first + 1..bodies.len() {
            let first_position = bodies[first].position.meters();
            let second_position = bodies[second].position.meters();
            let delta = [
                second_position[0] - first_position[0],
                second_position[1] - first_position[1],
                second_position[2] - first_position[2],
            ];
            let distance_squared = dot(delta, delta);
            if distance_squared == 0.0 {
                return Err(OrbitalError::CoincidentBodies);
            }
            let inverse_distance_cubed = 1.0 / (distance_squared * distance_squared.sqrt());
            for axis in 0..3 {
                values[first][axis] += GRAVITATIONAL_CONSTANT
                    * bodies[second].mass.kilograms()
                    * delta[axis]
                    * inverse_distance_cubed;
                values[second][axis] -= GRAVITATIONAL_CONSTANT
                    * bodies[first].mass.kilograms()
                    * delta[axis]
                    * inverse_distance_cubed;
            }
        }
    }
    Ok(values)
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

#[derive(Debug)]
pub enum OrbitalError {
    InvalidElements(String),
    InvalidPrimaryMass,
    InvalidStep,
    CoincidentBodies,
    KeplerDidNotConverge,
    Units(UnitError),
}

impl fmt::Display for OrbitalError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidElements(message) => write!(formatter, "invalid orbit: {message}"),
            Self::InvalidPrimaryMass => formatter.write_str("primary mass must be positive"),
            Self::InvalidStep => {
                formatter.write_str("integration step must be positive and finite")
            }
            Self::CoincidentBodies => {
                formatter.write_str("gravity bodies occupy the same position")
            }
            Self::KeplerDidNotConverge => formatter.write_str("Kepler solver did not converge"),
            Self::Units(error) => error.fmt(formatter),
        }
    }
}

impl Error for OrbitalError {}

impl From<UnitError> for OrbitalError {
    fn from(error: UnitError) -> Self {
        Self::Units(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::SECONDS_PER_YEAR;
    use crate::units::{KILOGRAMS_PER_EARTH_MASS, KILOGRAMS_PER_SOLAR_MASS};

    #[test]
    fn circular_kepler_orbit_returns_near_its_start_after_one_year() {
        let elements = OrbitalElements::new(
            Length::from_astronomical_units(1.0).unwrap(),
            0.0167,
            Angle::from_degrees(0.0).unwrap(),
            Angle::from_degrees(0.0).unwrap(),
            Angle::from_degrees(0.0).unwrap(),
            Angle::from_degrees(0.0).unwrap(),
            SimDuration::ZERO,
        )
        .unwrap();
        let origin = Position::default();
        let velocity = Velocity::default();
        let (start, _) = propagate_kepler(
            elements,
            Mass::from_solar_masses(1.0).unwrap(),
            SimDuration::ZERO,
            origin,
            velocity,
        )
        .unwrap();
        let (end, _) = propagate_kepler(
            elements,
            Mass::from_solar_masses(1.0).unwrap(),
            SimDuration::from_seconds(SECONDS_PER_YEAR),
            origin,
            velocity,
        )
        .unwrap();
        let start = start.meters();
        let end = end.meters();
        let error = ((start[0] - end[0]).powi(2) + (start[1] - end[1]).powi(2)).sqrt();
        assert!(error / elements.semi_major_axis().meters() < 0.01);
    }

    #[test]
    fn leapfrog_conserves_two_body_energy_and_angular_momentum() {
        let solar_mass = Mass::from_kilograms(KILOGRAMS_PER_SOLAR_MASS).unwrap();
        let earth_mass = Mass::from_kilograms(KILOGRAMS_PER_EARTH_MASS).unwrap();
        let radius = Length::from_astronomical_units(1.0).unwrap().meters();
        let speed = (GRAVITATIONAL_CONSTANT * solar_mass.kilograms() / radius).sqrt();
        let mut bodies = [
            GravityBody {
                mass: solar_mass,
                position: Position::default(),
                velocity: Velocity::default(),
            },
            GravityBody {
                mass: earth_mass,
                position: Position::from_meters(radius, 0.0, 0.0).unwrap(),
                velocity: Velocity::from_meters_per_second(0.0, speed, 0.0).unwrap(),
            },
        ];
        let initial_energy = total_energy(&bodies);
        let initial_angular_momentum = total_angular_momentum(&bodies);
        for _ in 0..365 {
            leapfrog_step(&mut bodies, 86_400.0).unwrap();
        }
        let position = bodies[1].position.meters();
        let final_radius = dot(position, position).sqrt();
        assert!((final_radius - radius).abs() / radius < 0.001);
        let energy_error = ((total_energy(&bodies) - initial_energy) / initial_energy).abs();
        let angular_momentum_error = ((total_angular_momentum(&bodies) - initial_angular_momentum)
            / initial_angular_momentum)
            .abs();
        assert!(energy_error < 1e-6, "relative energy error: {energy_error}");
        assert!(
            angular_momentum_error < 1e-12,
            "relative angular momentum error: {angular_momentum_error}"
        );
    }

    fn total_energy(bodies: &[GravityBody; 2]) -> f64 {
        let kinetic: f64 = bodies
            .iter()
            .map(|body| {
                0.5 * body.mass.kilograms()
                    * dot(
                        body.velocity.meters_per_second(),
                        body.velocity.meters_per_second(),
                    )
            })
            .sum();
        let first = bodies[0].position.meters();
        let second = bodies[1].position.meters();
        let separation = dot(
            [
                second[0] - first[0],
                second[1] - first[1],
                second[2] - first[2],
            ],
            [
                second[0] - first[0],
                second[1] - first[1],
                second[2] - first[2],
            ],
        )
        .sqrt();
        kinetic
            - GRAVITATIONAL_CONSTANT * bodies[0].mass.kilograms() * bodies[1].mass.kilograms()
                / separation
    }

    fn total_angular_momentum(bodies: &[GravityBody; 2]) -> f64 {
        bodies
            .iter()
            .map(|body| {
                let [x, y, _] = body.position.meters();
                let [vx, vy, _] = body.velocity.meters_per_second();
                body.mass.kilograms() * (x * vy - y * vx)
            })
            .sum()
    }
}
