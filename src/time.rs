use std::error::Error;
use std::fmt;
use std::str::FromStr;

pub const SECONDS_PER_MINUTE: u128 = 60;
pub const SECONDS_PER_HOUR: u128 = 60 * SECONDS_PER_MINUTE;
pub const SECONDS_PER_DAY: u128 = 24 * SECONDS_PER_HOUR;
pub const SECONDS_PER_YEAR: u128 = 31_557_600;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SimDuration {
    seconds: u128,
}

impl SimDuration {
    pub const ZERO: Self = Self { seconds: 0 };

    pub const fn from_seconds(seconds: u128) -> Self {
        Self { seconds }
    }

    pub fn from_years(years: u128) -> Result<Self, TimeError> {
        years
            .checked_mul(SECONDS_PER_YEAR)
            .map(Self::from_seconds)
            .ok_or(TimeError::Overflow)
    }

    pub const fn seconds(self) -> u128 {
        self.seconds
    }

    pub fn checked_mul(self, factor: u64) -> Option<Self> {
        self.seconds
            .checked_mul(u128::from(factor))
            .map(Self::from_seconds)
    }

    pub fn checked_add(self, other: Self) -> Option<Self> {
        self.seconds
            .checked_add(other.seconds)
            .map(Self::from_seconds)
    }
}

impl FromStr for SimDuration {
    type Err = DurationParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim();
        let number_end = value
            .find(|character: char| !character.is_ascii_digit() && character != '_')
            .unwrap_or(value.len());
        if number_end == 0 {
            return Err(DurationParseError::MissingNumber);
        }

        let number: u128 = value[..number_end]
            .replace('_', "")
            .parse()
            .map_err(|_| DurationParseError::InvalidNumber)?;
        let unit: String = value[number_end..]
            .chars()
            .filter(|character| !character.is_whitespace())
            .flat_map(char::to_lowercase)
            .collect();

        let multiplier = match unit.as_str() {
            "s" | "sec" | "secs" | "second" | "seconds" => 1,
            "min" | "mins" | "minute" | "minutes" => SECONDS_PER_MINUTE,
            "h" | "hr" | "hrs" | "hour" | "hours" => SECONDS_PER_HOUR,
            "d" | "day" | "days" => SECONDS_PER_DAY,
            "y" | "yr" | "yrs" | "year" | "years" => SECONDS_PER_YEAR,
            "kyr" | "kiloyear" | "kiloyears" | "thousandyears" => 1_000 * SECONDS_PER_YEAR,
            "myr" | "megayear" | "megayears" | "millionyears" => 1_000_000 * SECONDS_PER_YEAR,
            "gyr" | "gigayear" | "gigayears" | "billionyears" => 1_000_000_000 * SECONDS_PER_YEAR,
            "" => return Err(DurationParseError::MissingUnit),
            _ => return Err(DurationParseError::UnknownUnit(unit)),
        };

        number
            .checked_mul(multiplier)
            .map(Self::from_seconds)
            .ok_or(DurationParseError::Overflow)
    }
}

impl fmt::Display for SimDuration {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        const UNITS: [(&str, u128); 7] = [
            ("Gyr", 1_000_000_000 * SECONDS_PER_YEAR),
            ("Myr", 1_000_000 * SECONDS_PER_YEAR),
            ("kyr", 1_000 * SECONDS_PER_YEAR),
            ("years", SECONDS_PER_YEAR),
            ("days", SECONDS_PER_DAY),
            ("hours", SECONDS_PER_HOUR),
            ("minutes", SECONDS_PER_MINUTE),
        ];

        for (unit, seconds) in UNITS {
            if self.seconds >= seconds && self.seconds.is_multiple_of(seconds) {
                return write!(formatter, "{} {unit}", self.seconds / seconds);
            }
        }
        write!(formatter, "{} seconds", self.seconds)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TickDuration(SimDuration);

impl TickDuration {
    pub fn new(duration: SimDuration) -> Result<Self, TimeError> {
        if duration == SimDuration::ZERO {
            Err(TimeError::ZeroTickDuration)
        } else {
            Ok(Self(duration))
        }
    }

    pub const fn duration(self) -> SimDuration {
        self.0
    }
}

impl FromStr for TickDuration {
    type Err = TickDurationParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let duration = value.parse().map_err(TickDurationParseError::Duration)?;
        Self::new(duration).map_err(TickDurationParseError::Time)
    }
}

impl fmt::Display for TickDuration {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SimulationClock {
    tick: u64,
    elapsed: SimDuration,
    tick_duration: TickDuration,
}

impl SimulationClock {
    pub const fn new(tick_duration: TickDuration) -> Self {
        Self {
            tick: 0,
            elapsed: SimDuration::ZERO,
            tick_duration,
        }
    }

    pub fn from_parts(
        tick: u64,
        elapsed: SimDuration,
        tick_duration: TickDuration,
    ) -> Result<Self, TimeError> {
        let expected = tick_duration
            .duration()
            .checked_mul(tick)
            .ok_or(TimeError::Overflow)?;
        if expected != elapsed {
            return Err(TimeError::InconsistentClock);
        }
        Ok(Self {
            tick,
            elapsed,
            tick_duration,
        })
    }

    pub const fn tick(self) -> u64 {
        self.tick
    }

    pub const fn elapsed(self) -> SimDuration {
        self.elapsed
    }

    pub const fn tick_duration(self) -> TickDuration {
        self.tick_duration
    }

    pub fn advance(&mut self) -> Result<(), TimeError> {
        self.advance_by(1)
    }

    pub fn advance_by(&mut self, ticks: u64) -> Result<(), TimeError> {
        let tick = self.tick.checked_add(ticks).ok_or(TimeError::Overflow)?;
        let increment = self
            .tick_duration
            .duration()
            .checked_mul(ticks)
            .ok_or(TimeError::Overflow)?;
        let elapsed = self
            .elapsed
            .checked_add(increment)
            .ok_or(TimeError::Overflow)?;
        self.tick = tick;
        self.elapsed = elapsed;
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DurationParseError {
    MissingNumber,
    InvalidNumber,
    MissingUnit,
    UnknownUnit(String),
    Overflow,
}

impl fmt::Display for DurationParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingNumber => formatter.write_str("duration is missing a number"),
            Self::InvalidNumber => formatter.write_str("duration has an invalid number"),
            Self::MissingUnit => formatter.write_str("duration is missing a unit"),
            Self::UnknownUnit(unit) => write!(formatter, "unknown duration unit '{unit}'"),
            Self::Overflow => formatter.write_str("duration is too large"),
        }
    }
}

impl Error for DurationParseError {}

#[derive(Debug)]
pub enum TickDurationParseError {
    Duration(DurationParseError),
    Time(TimeError),
}

impl fmt::Display for TickDurationParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Duration(error) => error.fmt(formatter),
            Self::Time(error) => error.fmt(formatter),
        }
    }
}

impl Error for TickDurationParseError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeError {
    ZeroTickDuration,
    InconsistentClock,
    Overflow,
}

impl fmt::Display for TimeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroTickDuration => {
                formatter.write_str("tick duration must be greater than zero")
            }
            Self::InconsistentClock => {
                formatter.write_str("elapsed time does not match tick and tick duration")
            }
            Self::Overflow => formatter.write_str("simulation time overflow"),
        }
    }
}

impl Error for TimeError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_astronomical_durations() {
        assert_eq!(
            "1000 years".parse::<SimDuration>().unwrap(),
            SimDuration::from_years(1_000).unwrap()
        );
        assert_eq!(
            "1 Myr".parse::<SimDuration>().unwrap(),
            SimDuration::from_years(1_000_000).unwrap()
        );
        assert_eq!(
            "24h".parse::<SimDuration>().unwrap(),
            SimDuration::from_seconds(SECONDS_PER_DAY)
        );
    }

    #[test]
    fn clock_advances_canonically() {
        let tick_duration = "10 kyr".parse().unwrap();
        let mut clock = SimulationClock::new(tick_duration);
        clock.advance_by(1_000).unwrap();
        assert_eq!(clock.tick(), 1_000);
        assert_eq!(
            clock.elapsed(),
            SimDuration::from_years(10_000_000).unwrap()
        );
    }

    #[test]
    fn tick_duration_cannot_be_zero() {
        assert!(matches!(
            "0s".parse::<TickDuration>(),
            Err(TickDurationParseError::Time(TimeError::ZeroTickDuration))
        ));
    }
}
