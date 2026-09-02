use std::error::Error;
use std::fmt;
use std::fs;
use std::path::Path;

use crate::time::TickDuration;

pub const MAX_STAR_SYSTEMS: u32 = 100_000;
pub const MAX_GALAXIES: u16 = 1_024;
pub const MAX_PLANETS_PER_SYSTEM: u8 = 64;
pub const MAX_STAR_MASS_SOLAR: f64 = 300.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimulationConfig {
    pub seed: u64,
    pub tick_duration: TickDuration,
    pub generation: GenerationConfig,
}

impl SimulationConfig {
    pub const fn new(seed: u64, tick_duration: TickDuration) -> Self {
        Self {
            seed,
            tick_duration,
            generation: GenerationConfig::DEFAULT,
        }
    }

    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, ConfigLoadError> {
        let contents = fs::read_to_string(path).map_err(ConfigLoadError::Io)?;
        Self::parse(&contents).map_err(ConfigLoadError::Parse)
    }

    pub fn parse(contents: &str) -> Result<Self, ConfigParseError> {
        let mut section = None;
        let mut seed = None;
        let mut tick_duration = None;
        let mut generation = GenerationConfigBuilder::default();

        for (line_index, original_line) in contents.lines().enumerate() {
            let line_number = line_index + 1;
            let line = strip_comment(original_line).trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with('[') && line.ends_with(']') {
                let name = line[1..line.len() - 1].trim();
                if !matches!(name, "simulation" | "generation") {
                    return Err(ConfigParseError::new(
                        line_number,
                        format!("unknown section '[{name}]'"),
                    ));
                }
                section = Some(name);
                continue;
            }

            let (key, raw_value) = line.split_once('=').ok_or_else(|| {
                ConfigParseError::new(line_number, "expected a 'key = value' setting")
            })?;
            let key = key.trim();
            let raw_value = raw_value.trim();

            match section {
                Some("simulation") => match key {
                    "seed" => set_once(&mut seed, raw_value, line_number, key)?,
                    "tick_duration" => {
                        if tick_duration.is_some() {
                            return Err(ConfigParseError::new(
                                line_number,
                                "duplicate 'tick_duration'",
                            ));
                        }
                        let value = unquote(raw_value).ok_or_else(|| {
                            ConfigParseError::new(
                                line_number,
                                "tick_duration must be a quoted string",
                            )
                        })?;
                        tick_duration = Some(value.parse().map_err(|error| {
                            ConfigParseError::new(
                                line_number,
                                format!("invalid tick_duration: {error}"),
                            )
                        })?);
                    }
                    _ => return Err(unknown_setting(line_number, "simulation", key)),
                },
                Some("generation") => generation.set(key, raw_value, line_number)?,
                _ => {
                    return Err(ConfigParseError::new(
                        line_number,
                        "settings must appear under a section",
                    ));
                }
            }
        }

        let config = Self {
            seed: seed.ok_or_else(|| ConfigParseError::new(0, "missing 'simulation.seed'"))?,
            tick_duration: tick_duration
                .ok_or_else(|| ConfigParseError::new(0, "missing 'simulation.tick_duration'"))?,
            generation: generation.build(),
        };
        config
            .validate()
            .map_err(|error| ConfigParseError::new(0, error.to_string()))?;
        Ok(config)
    }

    pub fn validate(self) -> Result<(), ConfigValidationError> {
        self.generation.validate()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenerationConfig {
    pub galaxies: u16,
    pub star_systems: u32,
    pub planets_per_system_min: u8,
    pub planets_per_system_max: u8,
    pub star_mass_solar_min: f64,
    pub star_mass_solar_max: f64,
    pub planet_mass_earth_min: f64,
    pub planet_mass_earth_max: f64,
    pub orbit_au_min: f64,
    pub orbit_au_max: f64,
    pub galaxy_radius_parsecs: f64,
}

impl GenerationConfig {
    pub const DEFAULT: Self = Self {
        galaxies: 1,
        star_systems: 4,
        planets_per_system_min: 1,
        planets_per_system_max: 6,
        star_mass_solar_min: 0.3,
        star_mass_solar_max: MAX_STAR_MASS_SOLAR,
        planet_mass_earth_min: 0.1,
        planet_mass_earth_max: 300.0,
        orbit_au_min: 0.3,
        orbit_au_max: 30.0,
        galaxy_radius_parsecs: 25_000.0,
    };

    pub fn validate(self) -> Result<(), ConfigValidationError> {
        if self.galaxies == 0 || self.galaxies > MAX_GALAXIES {
            return Err(ConfigValidationError(format!(
                "generation.galaxies must be between 1 and {MAX_GALAXIES}"
            )));
        }
        if self.star_systems == 0 || self.star_systems > MAX_STAR_SYSTEMS {
            return Err(ConfigValidationError(format!(
                "generation.star_systems is the per-region population and must be between 1 and {MAX_STAR_SYSTEMS}; materialize additional regions for a larger universe-wide population"
            )));
        }
        if self.planets_per_system_min > self.planets_per_system_max {
            return Err(ConfigValidationError(
                "generation.planets_per_system_min cannot exceed planets_per_system_max".to_owned(),
            ));
        }
        if self.planets_per_system_max > MAX_PLANETS_PER_SYSTEM {
            return Err(ConfigValidationError(format!(
                "generation.planets_per_system_max cannot exceed {MAX_PLANETS_PER_SYSTEM}"
            )));
        }
        validate_range(
            "star_mass_solar",
            self.star_mass_solar_min,
            self.star_mass_solar_max,
        )?;
        if self.star_mass_solar_max > MAX_STAR_MASS_SOLAR {
            return Err(ConfigValidationError(format!(
                "generation.star_mass_solar_max cannot exceed {MAX_STAR_MASS_SOLAR} solar masses"
            )));
        }
        validate_range(
            "planet_mass_earth",
            self.planet_mass_earth_min,
            self.planet_mass_earth_max,
        )?;
        validate_range("orbit_au", self.orbit_au_min, self.orbit_au_max)?;
        if self.planets_per_system_max > 1 {
            let required_ratio = 1.35_f64.powi(i32::from(self.planets_per_system_max - 1));
            if self.orbit_au_max / self.orbit_au_min < required_ratio {
                return Err(ConfigValidationError(format!(
                    "generation orbit range is too narrow to stably space {} planets",
                    self.planets_per_system_max
                )));
            }
        }
        if !self.galaxy_radius_parsecs.is_finite() || self.galaxy_radius_parsecs <= 0.0 {
            return Err(ConfigValidationError(
                "generation.galaxy_radius_parsecs must be finite and greater than zero".to_owned(),
            ));
        }
        Ok(())
    }
}

impl Default for GenerationConfig {
    fn default() -> Self {
        Self::DEFAULT
    }
}

fn validate_range(name: &str, min: f64, max: f64) -> Result<(), ConfigValidationError> {
    if !min.is_finite() || !max.is_finite() || min <= 0.0 || max <= 0.0 {
        return Err(ConfigValidationError(format!(
            "generation.{name}_min and {name}_max must be finite and greater than zero"
        )));
    }
    if min > max {
        return Err(ConfigValidationError(format!(
            "generation.{name}_min cannot exceed {name}_max"
        )));
    }
    Ok(())
}

#[derive(Debug, Default)]
struct GenerationConfigBuilder {
    galaxies: Option<u16>,
    star_systems: Option<u32>,
    planets_per_system_min: Option<u8>,
    planets_per_system_max: Option<u8>,
    star_mass_solar_min: Option<f64>,
    star_mass_solar_max: Option<f64>,
    planet_mass_earth_min: Option<f64>,
    planet_mass_earth_max: Option<f64>,
    orbit_au_min: Option<f64>,
    orbit_au_max: Option<f64>,
    galaxy_radius_parsecs: Option<f64>,
}

impl GenerationConfigBuilder {
    fn set(&mut self, key: &str, raw_value: &str, line: usize) -> Result<(), ConfigParseError> {
        match key {
            "galaxies" => set_once(&mut self.galaxies, raw_value, line, key),
            "star_systems" => set_once(&mut self.star_systems, raw_value, line, key),
            "planets_per_system_min" => {
                set_once(&mut self.planets_per_system_min, raw_value, line, key)
            }
            "planets_per_system_max" => {
                set_once(&mut self.planets_per_system_max, raw_value, line, key)
            }
            "star_mass_solar_min" => set_once(&mut self.star_mass_solar_min, raw_value, line, key),
            "star_mass_solar_max" => set_once(&mut self.star_mass_solar_max, raw_value, line, key),
            "planet_mass_earth_min" => {
                set_once(&mut self.planet_mass_earth_min, raw_value, line, key)
            }
            "planet_mass_earth_max" => {
                set_once(&mut self.planet_mass_earth_max, raw_value, line, key)
            }
            "orbit_au_min" => set_once(&mut self.orbit_au_min, raw_value, line, key),
            "orbit_au_max" => set_once(&mut self.orbit_au_max, raw_value, line, key),
            "galaxy_radius_parsecs" => {
                set_once(&mut self.galaxy_radius_parsecs, raw_value, line, key)
            }
            _ => Err(unknown_setting(line, "generation", key)),
        }
    }

    fn build(self) -> GenerationConfig {
        let defaults = GenerationConfig::DEFAULT;
        GenerationConfig {
            galaxies: self.galaxies.unwrap_or(defaults.galaxies),
            star_systems: self.star_systems.unwrap_or(defaults.star_systems),
            planets_per_system_min: self
                .planets_per_system_min
                .unwrap_or(defaults.planets_per_system_min),
            planets_per_system_max: self
                .planets_per_system_max
                .unwrap_or(defaults.planets_per_system_max),
            star_mass_solar_min: self
                .star_mass_solar_min
                .unwrap_or(defaults.star_mass_solar_min),
            star_mass_solar_max: self
                .star_mass_solar_max
                .unwrap_or(defaults.star_mass_solar_max),
            planet_mass_earth_min: self
                .planet_mass_earth_min
                .unwrap_or(defaults.planet_mass_earth_min),
            planet_mass_earth_max: self
                .planet_mass_earth_max
                .unwrap_or(defaults.planet_mass_earth_max),
            orbit_au_min: self.orbit_au_min.unwrap_or(defaults.orbit_au_min),
            orbit_au_max: self.orbit_au_max.unwrap_or(defaults.orbit_au_max),
            galaxy_radius_parsecs: self
                .galaxy_radius_parsecs
                .unwrap_or(defaults.galaxy_radius_parsecs),
        }
    }
}

fn set_once<T>(
    slot: &mut Option<T>,
    raw_value: &str,
    line: usize,
    key: &str,
) -> Result<(), ConfigParseError>
where
    T: std::str::FromStr,
{
    if slot.is_some() {
        return Err(ConfigParseError::new(line, format!("duplicate '{key}'")));
    }
    *slot = Some(raw_value.parse().map_err(|_| {
        ConfigParseError::new(line, format!("'{key}' has an invalid numeric value"))
    })?);
    Ok(())
}

fn unknown_setting(line: usize, section: &str, key: &str) -> ConfigParseError {
    ConfigParseError::new(line, format!("unknown {section} setting '{key}'"))
}

fn strip_comment(line: &str) -> &str {
    let mut in_quotes = false;
    for (index, character) in line.char_indices() {
        match character {
            '"' => in_quotes = !in_quotes,
            '#' if !in_quotes => return &line[..index],
            _ => {}
        }
    }
    line
}

fn unquote(value: &str) -> Option<&str> {
    value.strip_prefix('"')?.strip_suffix('"')
}

#[derive(Debug)]
pub enum ConfigLoadError {
    Io(std::io::Error),
    Parse(ConfigParseError),
}

impl fmt::Display for ConfigLoadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "could not read configuration: {error}"),
            Self::Parse(error) => error.fmt(formatter),
        }
    }
}

impl Error for ConfigLoadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Parse(error) => Some(error),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConfigValidationError(String);

impl fmt::Display for ConfigValidationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for ConfigValidationError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConfigParseError {
    line: usize,
    message: String,
}

impl ConfigParseError {
    fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for ConfigParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(formatter, "invalid configuration: {}", self.message)
        } else {
            write!(
                formatter,
                "invalid configuration at line {}: {}",
                self.line, self.message
            )
        }
    }
}

impl Error for ConfigParseError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::SimDuration;

    #[test]
    fn parses_supported_configuration() {
        let config = SimulationConfig::parse(
            r#"
                # Initial deterministic universe
                [simulation]
                seed = 918273645
                tick_duration = "1000 years"

                [generation]
                galaxies = 3
                star_systems = 8
                planets_per_system_min = 0
                planets_per_system_max = 4
                orbit_au_max = 50.0
            "#,
        )
        .unwrap();

        assert_eq!(config.seed, 918_273_645);
        assert_eq!(
            config.tick_duration.duration(),
            SimDuration::from_years(1_000).unwrap()
        );
        assert_eq!(config.generation.star_systems, 8);
        assert_eq!(config.generation.galaxies, 3);
        assert_eq!(config.generation.planets_per_system_min, 0);
        assert_eq!(config.generation.orbit_au_max, 50.0);
    }

    #[test]
    fn generation_settings_are_optional() {
        let config = SimulationConfig::parse(
            r#"
                [simulation]
                seed = 1
                tick_duration = "1 year"
            "#,
        )
        .unwrap();
        assert_eq!(config.generation, GenerationConfig::DEFAULT);
    }

    #[test]
    fn rejects_unknown_settings_to_catch_typos() {
        let error = SimulationConfig::parse(
            r#"
                [simulation]
                seed = 1
                tick_duration = "1 year"
                tick_lenght = "2 years"
            "#,
        )
        .unwrap_err();
        assert!(error.to_string().contains("tick_lenght"));
    }

    #[test]
    fn rejects_invalid_generation_ranges() {
        let error = SimulationConfig::parse(
            r#"
                [simulation]
                seed = 1
                tick_duration = "1 year"
                [generation]
                orbit_au_min = 10
                orbit_au_max = 1
            "#,
        )
        .unwrap_err();
        assert!(error.to_string().contains("orbit_au_min"));
    }

    #[test]
    fn enforces_the_observed_upper_stellar_mass_limit() {
        let mut generation = GenerationConfig::DEFAULT;
        generation.star_mass_solar_max = MAX_STAR_MASS_SOLAR + 1.0;
        let error = generation.validate().unwrap_err();
        assert!(error.to_string().contains("300"));

        generation.star_mass_solar_max = MAX_STAR_MASS_SOLAR;
        generation.validate().unwrap();
    }

    #[test]
    fn validates_sparse_galaxy_catalog_size() {
        let mut generation = GenerationConfig::DEFAULT;
        generation.galaxies = 0;
        assert!(generation.validate().is_err());
        generation.galaxies = MAX_GALAXIES;
        generation.validate().unwrap();
    }

    #[test]
    fn excessive_region_population_explains_sparse_scaling() {
        let mut generation = GenerationConfig::DEFAULT;
        generation.star_systems = MAX_STAR_SYSTEMS + 1;
        let error = generation.validate().unwrap_err().to_string();
        assert!(error.contains("per-region population"));
        assert!(error.contains("materialize additional regions"));
    }
}
