use std::error::Error;
use std::fmt;
use std::str::FromStr;

use crate::id::GalaxyId;
use crate::objects::ObjectRegistry;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RegionCoord {
    pub x: i64,
    pub y: i64,
    pub z: i64,
}

impl RegionCoord {
    pub const ORIGIN: Self = Self { x: 0, y: 0, z: 0 };

    pub const fn new(x: i64, y: i64, z: i64) -> Self {
        Self { x, y, z }
    }
}

impl fmt::Display for RegionCoord {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{},{},{}", self.x, self.y, self.z)
    }
}

impl FromStr for RegionCoord {
    type Err = RegionCoordParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let mut components = value.split(',');
        let x = parse_component(components.next(), "x")?;
        let y = parse_component(components.next(), "y")?;
        let z = parse_component(components.next(), "z")?;
        if components.next().is_some() {
            return Err(RegionCoordParseError(
                "region must contain exactly three coordinates".to_owned(),
            ));
        }
        Ok(Self { x, y, z })
    }
}

fn parse_component(value: Option<&str>, name: &str) -> Result<i64, RegionCoordParseError> {
    value
        .ok_or_else(|| RegionCoordParseError(format!("region is missing its {name} coordinate")))?
        .trim()
        .parse()
        .map_err(|_| RegionCoordParseError(format!("region {name} coordinate is invalid")))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionCoordParseError(String);

impl fmt::Display for RegionCoordParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for RegionCoordParseError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GalaxyRegion {
    galaxy_id: GalaxyId,
    coordinate: RegionCoord,
}

impl GalaxyRegion {
    pub const fn new(galaxy_id: GalaxyId, coordinate: RegionCoord) -> Self {
        Self {
            galaxy_id,
            coordinate,
        }
    }

    pub const fn galaxy_id(self) -> GalaxyId {
        self.galaxy_id
    }

    pub const fn coordinate(self) -> RegionCoord {
        self.coordinate
    }
}

impl fmt::Display for GalaxyRegion {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}@{}", self.galaxy_id, self.coordinate)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MaterializedRegion {
    coordinate: RegionCoord,
    last_updated_tick: u64,
    objects: ObjectRegistry,
}

impl MaterializedRegion {
    pub(crate) const fn new(
        coordinate: RegionCoord,
        last_updated_tick: u64,
        objects: ObjectRegistry,
    ) -> Self {
        Self {
            coordinate,
            last_updated_tick,
            objects,
        }
    }

    pub const fn coordinate(&self) -> RegionCoord {
        self.coordinate
    }

    pub const fn last_updated_tick(&self) -> u64 {
        self.last_updated_tick
    }

    pub const fn objects(&self) -> &ObjectRegistry {
        &self.objects
    }

    pub(crate) fn objects_mut(&mut self) -> &mut ObjectRegistry {
        &mut self.objects
    }

    pub(crate) fn set_last_updated_tick(&mut self, tick: u64) {
        self.last_updated_tick = tick;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coordinates_round_trip_through_text() {
        let coordinate = RegionCoord::new(-18, 4, 91);
        assert_eq!(coordinate.to_string(), "-18,4,91");
        assert_eq!(coordinate.to_string().parse(), Ok(coordinate));
    }
}
