use std::error::Error;
use std::fmt;
use std::str::FromStr;

macro_rules! stable_id {
    ($name:ident, $prefix:literal) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u64);

        impl $name {
            pub const fn new(value: u64) -> Self {
                Self(value)
            }

            pub const fn value(self) -> u64 {
                self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(formatter, "{}-{}", $prefix, self.0)
            }
        }

        impl FromStr for $name {
            type Err = IdParseError;

            fn from_str(value: &str) -> Result<Self, Self::Err> {
                let number = value
                    .strip_prefix(concat!($prefix, "-"))
                    .ok_or(IdParseError::WrongPrefix { expected: $prefix })?;
                let number = number
                    .parse()
                    .map_err(|_| IdParseError::InvalidNumber(value.to_owned()))?;
                Ok(Self(number))
            }
        }
    };
}

stable_id!(UniverseId, "U");
stable_id!(GalaxyId, "G");
stable_id!(StarSystemId, "SYS");
stable_id!(StarId, "S");
stable_id!(PlanetId, "P");
stable_id!(MoonId, "M");
stable_id!(AsteroidId, "AST");
stable_id!(CometId, "COM");
stable_id!(AsteroidBeltId, "AB");
stable_id!(CometCloudId, "CC");
stable_id!(DebrisDiskId, "DD");
stable_id!(ProtoplanetaryDiskId, "PD");

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdParseError {
    WrongPrefix { expected: &'static str },
    InvalidNumber(String),
}

impl fmt::Display for IdParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongPrefix { expected } => {
                write!(formatter, "expected an ID beginning with '{expected}-'")
            }
            Self::InvalidNumber(value) => write!(formatter, "invalid numeric ID '{value}'"),
        }
    }
}

impl Error for IdParseError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_have_stable_text_representations() {
        let id = StarSystemId::new(42);
        assert_eq!(id.to_string(), "SYS-42");
        assert_eq!("SYS-42".parse::<StarSystemId>(), Ok(id));
        assert!("S-42".parse::<StarSystemId>().is_err());
    }
}
