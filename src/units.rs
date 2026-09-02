use std::error::Error;
use std::fmt;

pub const KILOGRAMS_PER_EARTH_MASS: f64 = 5.972_2e24;
pub const KILOGRAMS_PER_SOLAR_MASS: f64 = 1.988_47e30;
pub const METERS_PER_ASTRONOMICAL_UNIT: f64 = 149_597_870_700.0;
pub const METERS_PER_PARSEC: f64 = 3.085_677_581_491_367e16;
pub const WATTS_PER_SOLAR_LUMINOSITY: f64 = 3.828e26;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Mass(f64);

impl Mass {
    pub fn from_kilograms(value: f64) -> Result<Self, UnitError> {
        validate_nonnegative_finite("mass", value).map(Self)
    }

    pub fn from_earth_masses(value: f64) -> Result<Self, UnitError> {
        Self::from_kilograms(value * KILOGRAMS_PER_EARTH_MASS)
    }

    pub fn from_solar_masses(value: f64) -> Result<Self, UnitError> {
        Self::from_kilograms(value * KILOGRAMS_PER_SOLAR_MASS)
    }

    pub const fn kilograms(self) -> f64 {
        self.0
    }

    pub fn solar_masses(self) -> f64 {
        self.0 / KILOGRAMS_PER_SOLAR_MASS
    }

    pub fn earth_masses(self) -> f64 {
        self.0 / KILOGRAMS_PER_EARTH_MASS
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Length(f64);

impl Length {
    pub fn from_meters(value: f64) -> Result<Self, UnitError> {
        validate_nonnegative_finite("length", value).map(Self)
    }

    pub fn from_astronomical_units(value: f64) -> Result<Self, UnitError> {
        Self::from_meters(value * METERS_PER_ASTRONOMICAL_UNIT)
    }

    pub fn from_parsecs(value: f64) -> Result<Self, UnitError> {
        Self::from_meters(value * METERS_PER_PARSEC)
    }

    pub const fn meters(self) -> f64 {
        self.0
    }

    pub fn astronomical_units(self) -> f64 {
        self.0 / METERS_PER_ASTRONOMICAL_UNIT
    }

    pub fn parsecs(self) -> f64 {
        self.0 / METERS_PER_PARSEC
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Temperature(f64);

impl Temperature {
    pub fn from_kelvin(value: f64) -> Result<Self, UnitError> {
        validate_nonnegative_finite("temperature", value).map(Self)
    }

    pub const fn kelvin(self) -> f64 {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Luminosity(f64);

impl Luminosity {
    pub fn from_watts(value: f64) -> Result<Self, UnitError> {
        validate_nonnegative_finite("luminosity", value).map(Self)
    }

    pub fn from_solar_luminosities(value: f64) -> Result<Self, UnitError> {
        Self::from_watts(value * WATTS_PER_SOLAR_LUMINOSITY)
    }

    pub const fn watts(self) -> f64 {
        self.0
    }

    pub fn solar_luminosities(self) -> f64 {
        self.0 / WATTS_PER_SOLAR_LUMINOSITY
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Position {
    x_meters: f64,
    y_meters: f64,
    z_meters: f64,
}

impl Position {
    pub fn from_meters(x: f64, y: f64, z: f64) -> Result<Self, UnitError> {
        validate_finite("position", x)?;
        validate_finite("position", y)?;
        validate_finite("position", z)?;
        Ok(Self {
            x_meters: x,
            y_meters: y,
            z_meters: z,
        })
    }

    pub const fn meters(self) -> [f64; 3] {
        [self.x_meters, self.y_meters, self.z_meters]
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Velocity {
    x_meters_per_second: f64,
    y_meters_per_second: f64,
    z_meters_per_second: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, PartialOrd)]
pub struct Angle(f64);

impl Angle {
    pub fn from_radians(value: f64) -> Result<Self, UnitError> {
        validate_finite("angle", value).map(Self)
    }

    pub fn from_degrees(value: f64) -> Result<Self, UnitError> {
        Self::from_radians(value.to_radians())
    }

    pub const fn radians(self) -> f64 {
        self.0
    }

    pub fn degrees(self) -> f64 {
        self.0.to_degrees()
    }
}

impl Velocity {
    pub fn from_meters_per_second(x: f64, y: f64, z: f64) -> Result<Self, UnitError> {
        validate_finite("velocity", x)?;
        validate_finite("velocity", y)?;
        validate_finite("velocity", z)?;
        Ok(Self {
            x_meters_per_second: x,
            y_meters_per_second: y,
            z_meters_per_second: z,
        })
    }

    pub const fn meters_per_second(self) -> [f64; 3] {
        [
            self.x_meters_per_second,
            self.y_meters_per_second,
            self.z_meters_per_second,
        ]
    }
}

fn validate_nonnegative_finite(name: &'static str, value: f64) -> Result<f64, UnitError> {
    validate_finite(name, value)?;
    if value < 0.0 {
        Err(UnitError::Negative(name))
    } else {
        Ok(value)
    }
}

fn validate_finite(name: &'static str, value: f64) -> Result<f64, UnitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(UnitError::NotFinite(name))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UnitError {
    Negative(&'static str),
    NotFinite(&'static str),
}

impl fmt::Display for UnitError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Negative(name) => write!(formatter, "{name} cannot be negative"),
            Self::NotFinite(name) => write!(formatter, "{name} must be finite"),
        }
    }
}

impl Error for UnitError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn astronomical_units_convert_to_si() {
        let distance = Length::from_astronomical_units(1.0).unwrap();
        assert_eq!(distance.meters(), METERS_PER_ASTRONOMICAL_UNIT);

        let mass = Mass::from_solar_masses(1.0).unwrap();
        assert_eq!(mass.kilograms(), KILOGRAMS_PER_SOLAR_MASS);

        let luminosity = Luminosity::from_solar_luminosities(1.0).unwrap();
        assert_eq!(luminosity.watts(), WATTS_PER_SOLAR_LUMINOSITY);
    }

    #[test]
    fn units_reject_invalid_physical_values() {
        assert_eq!(Mass::from_kilograms(-1.0), Err(UnitError::Negative("mass")));
        assert_eq!(
            Temperature::from_kelvin(f64::NAN),
            Err(UnitError::NotFinite("temperature"))
        );
    }
}
