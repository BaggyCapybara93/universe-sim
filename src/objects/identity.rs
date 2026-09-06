use std::error::Error;
use std::fmt;
use std::str::FromStr;

use super::{AsteroidBelt, CometCloud, DebrisDisk, Galaxy, Moon, Planet, Star, StarSystem};
use crate::id::{
    AsteroidBeltId, CometCloudId, DebrisDiskId, GalaxyId, IdParseError, MoonId, PlanetId, StarId,
    StarSystemId,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjectId {
    Galaxy(GalaxyId),
    StarSystem(StarSystemId),
    Star(StarId),
    Planet(PlanetId),
    Moon(MoonId),
    AsteroidBelt(AsteroidBeltId),
    CometCloud(CometCloudId),
    DebrisDisk(DebrisDiskId),
}

impl fmt::Display for ObjectId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Galaxy(id) => id.fmt(formatter),
            Self::StarSystem(id) => id.fmt(formatter),
            Self::Star(id) => id.fmt(formatter),
            Self::Planet(id) => id.fmt(formatter),
            Self::Moon(id) => id.fmt(formatter),
            Self::AsteroidBelt(id) => id.fmt(formatter),
            Self::CometCloud(id) => id.fmt(formatter),
            Self::DebrisDisk(id) => id.fmt(formatter),
        }
    }
}

impl FromStr for ObjectId {
    type Err = ObjectIdParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.starts_with("SYS-") {
            value
                .parse()
                .map(Self::StarSystem)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("G-") {
            value
                .parse()
                .map(Self::Galaxy)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("S-") {
            value
                .parse()
                .map(Self::Star)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("P-") {
            value
                .parse()
                .map(Self::Planet)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("M-") {
            value
                .parse()
                .map(Self::Moon)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("AB-") {
            value
                .parse()
                .map(Self::AsteroidBelt)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("CC-") {
            value
                .parse()
                .map(Self::CometCloud)
                .map_err(ObjectIdParseError::Id)
        } else if value.starts_with("DD-") {
            value
                .parse()
                .map(Self::DebrisDisk)
                .map_err(ObjectIdParseError::Id)
        } else {
            Err(ObjectIdParseError::UnknownPrefix(value.to_owned()))
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ObjectIdParseError {
    UnknownPrefix(String),
    Id(IdParseError),
}

impl fmt::Display for ObjectIdParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownPrefix(value) => write!(
                formatter,
                "unknown object ID '{value}'; expected G-, SYS-, S-, P-, M-, AB-, CC-, or DD-"
            ),
            Self::Id(error) => error.fmt(formatter),
        }
    }
}

impl Error for ObjectIdParseError {}

#[derive(Debug, Clone, Copy)]
pub enum ObjectRef<'a> {
    Galaxy(&'a Galaxy),
    StarSystem(&'a StarSystem),
    Star(&'a Star),
    Planet(&'a Planet),
    Moon(&'a Moon),
    AsteroidBelt(&'a AsteroidBelt),
    CometCloud(&'a CometCloud),
    DebrisDisk(&'a DebrisDisk),
}

impl From<GalaxyId> for u64 {
    fn from(id: GalaxyId) -> Self {
        id.value()
    }
}

impl From<StarSystemId> for u64 {
    fn from(id: StarSystemId) -> Self {
        id.value()
    }
}

impl From<StarId> for u64 {
    fn from(id: StarId) -> Self {
        id.value()
    }
}

impl From<PlanetId> for u64 {
    fn from(id: PlanetId) -> Self {
        id.value()
    }
}

impl From<MoonId> for u64 {
    fn from(id: MoonId) -> Self {
        id.value()
    }
}

impl From<AsteroidBeltId> for u64 {
    fn from(id: AsteroidBeltId) -> Self {
        id.value()
    }
}

impl From<CometCloudId> for u64 {
    fn from(id: CometCloudId) -> Self {
        id.value()
    }
}

impl From<DebrisDiskId> for u64 {
    fn from(id: DebrisDiskId) -> Self {
        id.value()
    }
}
