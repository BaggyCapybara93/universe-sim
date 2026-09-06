use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use std::error::Error;
use std::fmt;
use std::str::FromStr;

use rayon::prelude::*;

use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, CometCloudId, DebrisDiskId, GalaxyId, IdParseError, MoonId, PlanetId, StarId,
    StarSystemId,
};
use crate::orbital::{OrbitalElements, OrbitalError, propagate_kepler};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration, TimeError};
use crate::units::{
    Angle, Length, Luminosity, METERS_PER_PARSEC, Mass, Position, Temperature, UnitError, Velocity,
};

mod dynamics;
mod persistence;
mod planetary;
mod stellar;

const MAX_SYSTEM_SPEED_METERS_PER_SECOND: f64 = 250_000.0;
const PARALLEL_SYSTEM_THRESHOLD: u32 = 128;
const PARALLEL_BODY_THRESHOLD: usize = 256;
const SOLAR_RADIUS_METERS: f64 = 695_700_000.0;

#[derive(Debug, Clone, PartialEq)]
pub struct Galaxy {
    id: GalaxyId,
    position: Position,
    star_system_ids: Vec<StarSystemId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GalaxyMorphology {
    Spiral,
    Elliptical,
    Irregular,
    Dwarf,
}

impl fmt::Display for GalaxyMorphology {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Spiral => "Spiral",
            Self::Elliptical => "Elliptical",
            Self::Irregular => "Irregular",
            Self::Dwarf => "Dwarf",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GalaxyCatalogEntry {
    id: GalaxyId,
    position: Position,
    morphology: GalaxyMorphology,
    radius_parsecs: f64,
}

impl GalaxyCatalogEntry {
    pub const fn id(self) -> GalaxyId {
        self.id
    }

    pub const fn position(self) -> Position {
        self.position
    }

    pub const fn morphology(self) -> GalaxyMorphology {
        self.morphology
    }

    pub const fn radius_parsecs(self) -> f64 {
        self.radius_parsecs
    }

    pub(crate) fn primary(seed: u64, radius_parsecs: f64) -> Self {
        Self {
            id: GalaxyId::new(derive_object_seed(seed, 0, "galaxy-id-v1")),
            position: Position::default(),
            morphology: GalaxyMorphology::Spiral,
            radius_parsecs,
        }
    }

    pub(crate) fn generate(
        seed: u64,
        index: u16,
        base_radius_parsecs: f64,
    ) -> Result<Self, ObjectError> {
        if index == 0 {
            return Ok(Self::primary(seed, base_radius_parsecs));
        }
        let id = GalaxyId::new(derive_object_seed(seed, u64::from(index), "galaxy-id-v2"));
        let mut rng = StableRng::new(derive_object_seed(
            seed,
            u64::from(index),
            "galaxy-catalog-v1",
        ));
        let morphology = match rng.next_f64() {
            value if value < 0.55 => GalaxyMorphology::Spiral,
            value if value < 0.75 => GalaxyMorphology::Elliptical,
            value if value < 0.90 => GalaxyMorphology::Irregular,
            _ => GalaxyMorphology::Dwarf,
        };
        let radius_factor = match morphology {
            GalaxyMorphology::Spiral => random_range(&mut rng, 0.7, 1.5),
            GalaxyMorphology::Elliptical => random_range(&mut rng, 0.5, 1.2),
            GalaxyMorphology::Irregular => random_range(&mut rng, 0.25, 0.7),
            GalaxyMorphology::Dwarf => random_range(&mut rng, 0.05, 0.25),
        };
        let distance = random_range(&mut rng, 500_000.0, 10_000_000.0) * METERS_PER_PARSEC;
        let direction = random_unit_vector(&mut rng);
        Ok(Self {
            id,
            position: Position::from_meters(
                distance * direction[0],
                distance * direction[1],
                distance * direction[2],
            )?,
            morphology,
            radius_parsecs: base_radius_parsecs * radius_factor,
        })
    }
}

impl Galaxy {
    pub const fn id(&self) -> GalaxyId {
        self.id
    }

    pub const fn position(&self) -> Position {
        self.position
    }

    pub fn star_system_ids(&self) -> &[StarSystemId] {
        &self.star_system_ids
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StarSystem {
    id: StarSystemId,
    galaxy_id: GalaxyId,
    position: Position,
    velocity: Velocity,
    star_id: StarId,
    companion_star_ids: Vec<StarId>,
    planet_ids: Vec<PlanetId>,
    moon_ids: Vec<MoonId>,
    asteroid_belt_ids: Vec<AsteroidBeltId>,
    comet_cloud_ids: Vec<CometCloudId>,
    debris_disk_ids: Vec<DebrisDiskId>,
    component: GalacticComponent,
}

impl StarSystem {
    pub const fn id(&self) -> StarSystemId {
        self.id
    }

    pub const fn galaxy_id(&self) -> GalaxyId {
        self.galaxy_id
    }

    pub const fn position(&self) -> Position {
        self.position
    }

    pub const fn velocity(&self) -> Velocity {
        self.velocity
    }

    pub const fn star_id(&self) -> StarId {
        self.star_id
    }

    pub fn companion_star_ids(&self) -> &[StarId] {
        &self.companion_star_ids
    }

    pub fn planet_ids(&self) -> &[PlanetId] {
        &self.planet_ids
    }

    pub fn moon_ids(&self) -> &[MoonId] {
        &self.moon_ids
    }

    pub fn asteroid_belt_ids(&self) -> &[AsteroidBeltId] {
        &self.asteroid_belt_ids
    }

    pub fn comet_cloud_ids(&self) -> &[CometCloudId] {
        &self.comet_cloud_ids
    }

    pub fn debris_disk_ids(&self) -> &[DebrisDiskId] {
        &self.debris_disk_ids
    }

    pub const fn component(&self) -> GalacticComponent {
        self.component
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GalacticComponent {
    Disk,
    Bulge,
    Halo,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DynamicsMode {
    Keplerian,
    Interacting,
}

impl fmt::Display for DynamicsMode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Keplerian => "Keplerian",
            Self::Interacting => "Interacting N-body",
        })
    }
}

impl fmt::Display for GalacticComponent {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Disk => "Disk",
            Self::Bulge => "Bulge",
            Self::Halo => "Halo",
        })
    }
}

impl GalacticComponent {
    fn code(self) -> &'static str {
        match self {
            Self::Disk => "disk",
            Self::Bulge => "bulge",
            Self::Halo => "halo",
        }
    }

    fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "disk" => Ok(Self::Disk),
            "bulge" => Ok(Self::Bulge),
            "halo" => Ok(Self::Halo),
            _ => Err(ObjectError(format!("unknown galactic component '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Star {
    id: StarId,
    system_id: StarSystemId,
    initial_mass: Mass,
    current_mass: Mass,
    epoch_age: SimDuration,
    formed_at: SimDuration,
    age: SimDuration,
    metallicity: f64,
    luminosity: Luminosity,
    temperature: Temperature,
    stage: StellarStage,
    position: Position,
    velocity: Velocity,
    orbital_elements: Option<OrbitalElements>,
}

impl Star {
    pub const fn id(&self) -> StarId {
        self.id
    }

    pub const fn system_id(&self) -> StarSystemId {
        self.system_id
    }

    pub const fn mass(&self) -> Mass {
        self.current_mass
    }

    pub const fn initial_mass(&self) -> Mass {
        self.initial_mass
    }

    pub const fn current_mass(&self) -> Mass {
        self.current_mass
    }

    pub const fn epoch_age(&self) -> SimDuration {
        self.epoch_age
    }

    pub const fn formed_at(&self) -> SimDuration {
        self.formed_at
    }

    pub const fn age(&self) -> SimDuration {
        self.age
    }

    pub const fn metallicity(&self) -> f64 {
        self.metallicity
    }

    pub const fn luminosity(&self) -> Luminosity {
        self.luminosity
    }

    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }

    pub const fn stage(&self) -> StellarStage {
        self.stage
    }

    pub const fn position(&self) -> Position {
        self.position
    }

    pub const fn velocity(&self) -> Velocity {
        self.velocity
    }

    pub const fn orbital_elements(&self) -> Option<OrbitalElements> {
        self.orbital_elements
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum StellarStage {
    Protostar,
    PreMainSequence,
    MainSequence,
    Subgiant,
    Giant,
    CoreHeliumBurning,
    AsymptoticGiantBranch,
    PlanetaryNebula,
    Supergiant,
    WolfRayet,
    PulsationalPairInstability,
    CoreCollapseSupernova,
    WhiteDwarf,
    NeutronStar,
    BlackHole,
    Disrupted,
}

impl fmt::Display for StellarStage {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Protostar => "Protostar",
            Self::PreMainSequence => "Pre-Main Sequence",
            Self::MainSequence => "Main Sequence",
            Self::Subgiant => "Subgiant",
            Self::Giant => "Giant",
            Self::CoreHeliumBurning => "Core Helium Burning",
            Self::AsymptoticGiantBranch => "Asymptotic Giant Branch",
            Self::PlanetaryNebula => "Planetary Nebula",
            Self::Supergiant => "Supergiant",
            Self::WolfRayet => "Wolf-Rayet",
            Self::PulsationalPairInstability => "Pulsational Pair Instability",
            Self::CoreCollapseSupernova => "Core-Collapse Supernova",
            Self::WhiteDwarf => "White Dwarf",
            Self::NeutronStar => "Neutron Star",
            Self::BlackHole => "Black Hole",
            Self::Disrupted => "Disrupted",
        })
    }
}

impl StellarStage {
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::Protostar => "protostar",
            Self::PreMainSequence => "pre_main_sequence",
            Self::MainSequence => "main_sequence",
            Self::Subgiant => "subgiant",
            Self::Giant => "giant",
            Self::CoreHeliumBurning => "core_helium_burning",
            Self::AsymptoticGiantBranch => "asymptotic_giant_branch",
            Self::PlanetaryNebula => "planetary_nebula",
            Self::Supergiant => "supergiant",
            Self::WolfRayet => "wolf_rayet",
            Self::PulsationalPairInstability => "pulsational_pair_instability",
            Self::CoreCollapseSupernova => "core_collapse_supernova",
            Self::WhiteDwarf => "white_dwarf",
            Self::NeutronStar => "neutron_star",
            Self::BlackHole => "black_hole",
            Self::Disrupted => "disrupted",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "protostar" => Ok(Self::Protostar),
            "pre_main_sequence" => Ok(Self::PreMainSequence),
            "main_sequence" => Ok(Self::MainSequence),
            "subgiant" => Ok(Self::Subgiant),
            "giant" => Ok(Self::Giant),
            "core_helium_burning" => Ok(Self::CoreHeliumBurning),
            "asymptotic_giant_branch" => Ok(Self::AsymptoticGiantBranch),
            "planetary_nebula" => Ok(Self::PlanetaryNebula),
            "supergiant" => Ok(Self::Supergiant),
            "wolf_rayet" => Ok(Self::WolfRayet),
            "pulsational_pair_instability" => Ok(Self::PulsationalPairInstability),
            "core_collapse_supernova" => Ok(Self::CoreCollapseSupernova),
            "white_dwarf" => Ok(Self::WhiteDwarf),
            "neutron_star" => Ok(Self::NeutronStar),
            "black_hole" => Ok(Self::BlackHole),
            "disrupted" => Ok(Self::Disrupted),
            _ => Err(ObjectError(format!("unknown stellar stage '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StellarTransition {
    pub tick: u64,
    pub star_id: StarId,
    pub from: StellarStage,
    pub to: StellarStage,
    pub impacts: Vec<PlanetImpact>,
    pub natal_kick_meters_per_second: u32,
    pub interaction: Option<StellarInteractionKind>,
    pub related_star_id: Option<StarId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StellarInteractionKind {
    StarFormation,
    MassTransfer,
    StellarMerger,
    PlanetaryCollision,
}

impl StellarInteractionKind {
    pub(crate) const fn code(self) -> &'static str {
        match self {
            Self::StarFormation => "star_formation",
            Self::MassTransfer => "mass_transfer",
            Self::StellarMerger => "stellar_merger",
            Self::PlanetaryCollision => "planetary_collision",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "star_formation" => Ok(Self::StarFormation),
            "mass_transfer" => Ok(Self::MassTransfer),
            "stellar_merger" => Ok(Self::StellarMerger),
            "planetary_collision" => Ok(Self::PlanetaryCollision),
            _ => Err(ObjectError(format!(
                "unknown stellar interaction '{value}'"
            ))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlanetImpact {
    pub planet_id: PlanetId,
    pub kind: PlanetImpactKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanetImpactKind {
    Irradiated,
    Unbound,
    Destroyed,
}

impl fmt::Display for PlanetImpactKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Irradiated => "irradiated",
            Self::Unbound => "unbound",
            Self::Destroyed => "destroyed",
        })
    }
}

impl PlanetImpactKind {
    pub(crate) const fn code(self) -> &'static str {
        match self {
            Self::Irradiated => "irradiated",
            Self::Unbound => "unbound",
            Self::Destroyed => "destroyed",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "irradiated" => Ok(Self::Irradiated),
            "unbound" => Ok(Self::Unbound),
            "destroyed" => Ok(Self::Destroyed),
            _ => Err(ObjectError(format!("unknown planetary impact '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanetStatus {
    Stable,
    Irradiated,
    Unbound,
    Destroyed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanetAtmosphere {
    None,
    Thin,
    NitrogenRich,
    CarbonDioxideRich,
    HydrogenHelium,
}

impl fmt::Display for PlanetAtmosphere {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::None => "None",
            Self::Thin => "Thin",
            Self::NitrogenRich => "Nitrogen-rich",
            Self::CarbonDioxideRich => "Carbon-dioxide-rich",
            Self::HydrogenHelium => "Hydrogen-helium",
        })
    }
}

impl PlanetAtmosphere {
    fn code(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Thin => "thin",
            Self::NitrogenRich => "nitrogen_rich",
            Self::CarbonDioxideRich => "carbon_dioxide_rich",
            Self::HydrogenHelium => "hydrogen_helium",
        }
    }

    fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "none" => Ok(Self::None),
            "thin" => Ok(Self::Thin),
            "nitrogen_rich" => Ok(Self::NitrogenRich),
            "carbon_dioxide_rich" => Ok(Self::CarbonDioxideRich),
            "hydrogen_helium" => Ok(Self::HydrogenHelium),
            _ => Err(ObjectError(format!(
                "unknown planetary atmosphere '{value}'"
            ))),
        }
    }
}

impl fmt::Display for PlanetStatus {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Stable => "Stable",
            Self::Irradiated => "Irradiated",
            Self::Unbound => "Unbound",
            Self::Destroyed => "Destroyed",
        })
    }
}

impl PlanetStatus {
    fn code(self) -> &'static str {
        match self {
            Self::Stable => "stable",
            Self::Irradiated => "irradiated",
            Self::Unbound => "unbound",
            Self::Destroyed => "destroyed",
        }
    }

    fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "stable" => Ok(Self::Stable),
            "irradiated" => Ok(Self::Irradiated),
            "unbound" => Ok(Self::Unbound),
            "destroyed" => Ok(Self::Destroyed),
            _ => Err(ObjectError(format!("unknown planet status '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Planet {
    id: PlanetId,
    system_id: StarSystemId,
    star_id: StarId,
    mass: Mass,
    radius: Length,
    atmosphere: PlanetAtmosphere,
    equilibrium_temperature: Temperature,
    orbital_radius: Length,
    position: Position,
    velocity: Velocity,
    orbital_elements: OrbitalElements,
    status: PlanetStatus,
}

impl Planet {
    pub const fn id(&self) -> PlanetId {
        self.id
    }

    pub const fn system_id(&self) -> StarSystemId {
        self.system_id
    }

    pub const fn star_id(&self) -> StarId {
        self.star_id
    }

    pub const fn mass(&self) -> Mass {
        self.mass
    }

    pub const fn radius(&self) -> Length {
        self.radius
    }

    pub const fn atmosphere(&self) -> PlanetAtmosphere {
        self.atmosphere
    }

    pub const fn equilibrium_temperature(&self) -> Temperature {
        self.equilibrium_temperature
    }

    pub fn is_potentially_habitable(&self) -> bool {
        matches!(
            self.atmosphere,
            PlanetAtmosphere::NitrogenRich | PlanetAtmosphere::CarbonDioxideRich
        ) && (240.0..=330.0).contains(&self.equilibrium_temperature.kelvin())
            && matches!(self.status, PlanetStatus::Stable | PlanetStatus::Irradiated)
    }

    pub const fn orbital_radius(&self) -> Length {
        self.orbital_radius
    }

    pub const fn position(&self) -> Position {
        self.position
    }

    pub const fn velocity(&self) -> Velocity {
        self.velocity
    }

    pub const fn orbital_elements(&self) -> OrbitalElements {
        self.orbital_elements
    }

    pub const fn status(&self) -> PlanetStatus {
        self.status
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Moon {
    id: MoonId,
    planet_id: PlanetId,
    system_id: StarSystemId,
    mass: Mass,
    radius: Length,
    orbital_radius: Length,
    position: Position,
    velocity: Velocity,
    orbital_elements: OrbitalElements,
}

impl Moon {
    pub const fn id(&self) -> MoonId {
        self.id
    }
    pub const fn planet_id(&self) -> PlanetId {
        self.planet_id
    }
    pub const fn system_id(&self) -> StarSystemId {
        self.system_id
    }
    pub const fn mass(&self) -> Mass {
        self.mass
    }

    pub const fn radius(&self) -> Length {
        self.radius
    }
    pub const fn orbital_radius(&self) -> Length {
        self.orbital_radius
    }
    pub const fn position(&self) -> Position {
        self.position
    }
    pub const fn velocity(&self) -> Velocity {
        self.velocity
    }
    pub const fn orbital_elements(&self) -> OrbitalElements {
        self.orbital_elements
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AsteroidBelt {
    id: AsteroidBeltId,
    system_id: StarSystemId,
    inner_radius: Length,
    outer_radius: Length,
    estimated_body_count: u64,
    total_mass: Mass,
}

impl AsteroidBelt {
    pub const fn id(&self) -> AsteroidBeltId {
        self.id
    }
    pub const fn system_id(&self) -> StarSystemId {
        self.system_id
    }
    pub const fn inner_radius(&self) -> Length {
        self.inner_radius
    }
    pub const fn outer_radius(&self) -> Length {
        self.outer_radius
    }
    pub const fn estimated_body_count(&self) -> u64 {
        self.estimated_body_count
    }
    pub const fn total_mass(&self) -> Mass {
        self.total_mass
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CometCloud {
    id: CometCloudId,
    system_id: StarSystemId,
    inner_radius: Length,
    outer_radius: Length,
    estimated_body_count: u64,
    total_mass: Mass,
}

impl CometCloud {
    pub const fn id(&self) -> CometCloudId {
        self.id
    }
    pub const fn system_id(&self) -> StarSystemId {
        self.system_id
    }
    pub const fn inner_radius(&self) -> Length {
        self.inner_radius
    }
    pub const fn outer_radius(&self) -> Length {
        self.outer_radius
    }
    pub const fn estimated_body_count(&self) -> u64 {
        self.estimated_body_count
    }
    pub const fn total_mass(&self) -> Mass {
        self.total_mass
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DebrisDisk {
    id: DebrisDiskId,
    system_id: StarSystemId,
    inner_radius: Length,
    outer_radius: Length,
    estimated_body_count: u64,
    total_mass: Mass,
}

impl DebrisDisk {
    pub const fn id(&self) -> DebrisDiskId {
        self.id
    }
    pub const fn system_id(&self) -> StarSystemId {
        self.system_id
    }
    pub const fn inner_radius(&self) -> Length {
        self.inner_radius
    }
    pub const fn outer_radius(&self) -> Length {
        self.outer_radius
    }
    pub const fn estimated_body_count(&self) -> u64 {
        self.estimated_body_count
    }
    pub const fn total_mass(&self) -> Mass {
        self.total_mass
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObjectCounts {
    pub galaxies: usize,
    pub star_systems: usize,
    pub stars: usize,
    pub planets: usize,
    pub moons: usize,
    pub asteroid_belts: usize,
    pub comet_clouds: usize,
    pub debris_disks: usize,
    pub gas_reservoirs: usize,
}

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

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ObjectRegistry {
    galaxies: BTreeMap<GalaxyId, Galaxy>,
    star_systems: BTreeMap<StarSystemId, StarSystem>,
    stars: BTreeMap<StarId, Star>,
    planets: BTreeMap<PlanetId, Planet>,
    moons: BTreeMap<MoonId, Moon>,
    asteroid_belts: BTreeMap<AsteroidBeltId, AsteroidBelt>,
    comet_clouds: BTreeMap<CometCloudId, CometCloud>,
    debris_disks: BTreeMap<DebrisDiskId, DebrisDisk>,
    gas_reservoir_solar_masses: f64,
    gas_composition: GasComposition,
    formation_sequence: u64,
}

impl ObjectRegistry {
    pub fn generate(seed: u64, config: GenerationConfig) -> Result<Self, ObjectError> {
        let galaxy = GalaxyCatalogEntry::primary(seed, config.galaxy_radius_parsecs);
        Self::generate_in_galaxy(seed, config, galaxy)
    }

    pub(crate) fn generate_in_galaxy(
        seed: u64,
        config: GenerationConfig,
        galaxy: GalaxyCatalogEntry,
    ) -> Result<Self, ObjectError> {
        config
            .validate()
            .map_err(|error| ObjectError(error.to_string()))?;

        let galaxy_id = galaxy.id;
        let mut registry = Self::default();
        let mut system_ids = Vec::with_capacity(config.star_systems as usize);

        let generate = |system_index| generate_system(seed, config, galaxy, system_index);
        let generated: Result<Vec<_>, ObjectError> =
            if config.star_systems >= PARALLEL_SYSTEM_THRESHOLD {
                (0..config.star_systems)
                    .into_par_iter()
                    .map(generate)
                    .collect()
            } else {
                (0..config.star_systems).map(generate).collect()
            };
        for generated in generated? {
            let system_id = generated.system.id;
            for star in generated.stars {
                insert_unique(&mut registry.stars, star.id, star, "star")?;
            }
            for planet in generated.planets {
                insert_unique(&mut registry.planets, planet.id, planet, "planet")?;
            }
            for moon in generated.moons {
                insert_unique(&mut registry.moons, moon.id, moon, "moon")?;
            }
            for belt in generated.asteroid_belts {
                insert_unique(&mut registry.asteroid_belts, belt.id, belt, "asteroid belt")?;
            }
            for cloud in generated.comet_clouds {
                insert_unique(&mut registry.comet_clouds, cloud.id, cloud, "comet cloud")?;
            }
            for disk in generated.debris_disks {
                insert_unique(&mut registry.debris_disks, disk.id, disk, "debris disk")?;
            }
            insert_unique(
                &mut registry.star_systems,
                system_id,
                generated.system,
                "star system",
            )?;
            system_ids.push(system_id);
        }

        let gas_factor = match galaxy.morphology {
            GalaxyMorphology::Spiral => 2.0,
            GalaxyMorphology::Elliptical => 0.25,
            GalaxyMorphology::Irregular => 3.0,
            GalaxyMorphology::Dwarf => 1.0,
        };
        let gas_metallicity = match galaxy.morphology {
            GalaxyMorphology::Spiral => 0.014,
            GalaxyMorphology::Elliptical => 0.020,
            GalaxyMorphology::Irregular => 0.006,
            GalaxyMorphology::Dwarf => 0.002,
        };
        let galaxy = Galaxy {
            id: galaxy_id,
            position: galaxy.position,
            star_system_ids: system_ids,
        };
        insert_unique(&mut registry.galaxies, galaxy_id, galaxy, "galaxy")?;
        registry.gas_reservoir_solar_masses = f64::from(config.star_systems) * gas_factor;
        registry.gas_composition = GasComposition::from_metallicity(gas_metallicity)
            .map_err(|error| ObjectError(error.to_string()))?;
        registry.validate()?;
        Ok(registry)
    }

    pub fn counts(&self) -> ObjectCounts {
        ObjectCounts {
            galaxies: self.galaxies.len(),
            star_systems: self.star_systems.len(),
            stars: self.stars.len(),
            planets: self.planets.len(),
            moons: self.moons.len(),
            asteroid_belts: self.asteroid_belts.len(),
            comet_clouds: self.comet_clouds.len(),
            debris_disks: self.debris_disks.len(),
            gas_reservoirs: usize::from(self.gas_reservoir_solar_masses > 0.0),
        }
    }

    pub const fn gas_reservoir_solar_masses(&self) -> f64 {
        self.gas_reservoir_solar_masses
    }

    pub const fn gas_composition(&self) -> GasComposition {
        self.gas_composition
    }

    pub(crate) fn form_stars_between(
        &mut self,
        seed: u64,
        config: GenerationConfig,
        from: SimDuration,
        to: SimDuration,
        tick_duration: SimDuration,
    ) -> Result<Vec<StellarTransition>, ObjectError> {
        const FORMATION_INTERVAL_YEARS: u128 = 10_000_000;
        const STAR_FORMATION_EFFICIENCY: f64 = 0.30;
        let interval = SimDuration::from_years(FORMATION_INTERVAL_YEARS)
            .map_err(|error| ObjectError(error.to_string()))?;
        let galaxy_id = self
            .galaxies
            .keys()
            .next()
            .copied()
            .ok_or_else(|| ObjectError("star formation requires a galaxy".to_owned()))?;
        let mut events = Vec::new();
        loop {
            let sequence = self.formation_sequence.saturating_add(1);
            let formation_time = interval
                .checked_mul(sequence)
                .ok_or_else(|| ObjectError("star-formation time overflow".to_owned()))?;
            if formation_time > to {
                break;
            }
            self.formation_sequence = sequence;
            if formation_time <= from {
                continue;
            }
            let mut rng = StableRng::new(derive_object_seed(
                seed,
                sequence,
                "dynamic-star-formation-v1",
            ));
            let available_stellar_mass =
                self.gas_reservoir_solar_masses * STAR_FORMATION_EFFICIENCY;
            let maximum_mass = config.star_mass_solar_max.min(available_stellar_mass);
            if maximum_mass < config.star_mass_solar_min {
                break;
            }
            let mass_solar =
                sample_stellar_mass(&mut rng, config.star_mass_solar_min, maximum_mass);
            let gas_consumed = mass_solar / STAR_FORMATION_EFFICIENCY;
            let system_id = StarSystemId::new(derive_object_seed(
                seed,
                sequence,
                "formed-star-system-id-v1",
            ));
            let star_id = StarId::new(derive_object_seed(seed, sequence, "formed-star-id-v1"));
            if self.star_systems.contains_key(&system_id) || self.stars.contains_key(&star_id) {
                return Err(ObjectError(
                    "dynamic star-formation ID collision".to_owned(),
                ));
            }
            let radius = config.galaxy_radius_parsecs * METERS_PER_PARSEC * rng.next_f64().sqrt();
            let angle = std::f64::consts::TAU * rng.next_f64();
            let (sin_angle, cos_angle) = angle.sin_cos();
            let [gx, gy, gz] = self
                .galaxies
                .get(&galaxy_id)
                .expect("formation galaxy exists")
                .position
                .meters();
            let position =
                Position::from_meters(gx + radius * cos_angle, gy + radius * sin_angle, gz)?;
            let velocity = Velocity::from_meters_per_second(
                -220_000.0 * sin_angle,
                220_000.0 * cos_angle,
                0.0,
            )?;
            let metallicity = (0.004 + 0.020 * rng.next_f64()).clamp(0.0001, 0.04);
            let mut star = generate_star_with_birth(
                star_id,
                system_id,
                formation_time,
                metallicity,
                mass_solar,
            )?;
            star.position = position;
            star.velocity = velocity;
            self.stars.insert(star_id, star);
            self.star_systems.insert(
                system_id,
                StarSystem {
                    id: system_id,
                    galaxy_id,
                    position,
                    velocity,
                    star_id,
                    companion_star_ids: Vec::new(),
                    planet_ids: Vec::new(),
                    moon_ids: Vec::new(),
                    asteroid_belt_ids: Vec::new(),
                    comet_cloud_ids: Vec::new(),
                    debris_disk_ids: Vec::new(),
                    component: GalacticComponent::Disk,
                },
            );
            self.galaxies
                .get_mut(&galaxy_id)
                .expect("formation galaxy exists")
                .star_system_ids
                .push(system_id);
            self.gas_reservoir_solar_masses =
                (self.gas_reservoir_solar_masses - gas_consumed).max(0.0);
            events.push(StellarTransition {
                tick: tick_at_or_after(formation_time, tick_duration)
                    .map_err(|error| ObjectError(error.to_string()))?,
                star_id,
                from: StellarStage::Protostar,
                to: StellarStage::Protostar,
                impacts: Vec::new(),
                natal_kick_meters_per_second: 0,
                interaction: Some(StellarInteractionKind::StarFormation),
                related_star_id: None,
            });
        }
        self.validate()?;
        Ok(events)
    }

    pub fn galaxy(&self, id: GalaxyId) -> Option<&Galaxy> {
        self.galaxies.get(&id)
    }

    pub fn star_system(&self, id: StarSystemId) -> Option<&StarSystem> {
        self.star_systems.get(&id)
    }

    pub fn star(&self, id: StarId) -> Option<&Star> {
        self.stars.get(&id)
    }

    pub fn planet(&self, id: PlanetId) -> Option<&Planet> {
        self.planets.get(&id)
    }

    pub fn moon(&self, id: MoonId) -> Option<&Moon> {
        self.moons.get(&id)
    }

    pub fn asteroid_belt(&self, id: AsteroidBeltId) -> Option<&AsteroidBelt> {
        self.asteroid_belts.get(&id)
    }

    pub fn comet_cloud(&self, id: CometCloudId) -> Option<&CometCloud> {
        self.comet_clouds.get(&id)
    }

    pub fn debris_disk(&self, id: DebrisDiskId) -> Option<&DebrisDisk> {
        self.debris_disks.get(&id)
    }

    pub fn find(&self, id: ObjectId) -> Option<ObjectRef<'_>> {
        match id {
            ObjectId::Galaxy(id) => self.galaxy(id).map(ObjectRef::Galaxy),
            ObjectId::StarSystem(id) => self.star_system(id).map(ObjectRef::StarSystem),
            ObjectId::Star(id) => self.star(id).map(ObjectRef::Star),
            ObjectId::Planet(id) => self.planet(id).map(ObjectRef::Planet),
            ObjectId::Moon(id) => self.moon(id).map(ObjectRef::Moon),
            ObjectId::AsteroidBelt(id) => self.asteroid_belt(id).map(ObjectRef::AsteroidBelt),
            ObjectId::CometCloud(id) => self.comet_cloud(id).map(ObjectRef::CometCloud),
            ObjectId::DebrisDisk(id) => self.debris_disk(id).map(ObjectRef::DebrisDisk),
        }
    }

    pub fn galaxies(&self) -> impl Iterator<Item = &Galaxy> {
        self.galaxies.values()
    }

    pub fn star_systems(&self) -> impl Iterator<Item = &StarSystem> {
        self.star_systems.values()
    }

    pub fn stars(&self) -> impl Iterator<Item = &Star> {
        self.stars.values()
    }

    pub fn planets(&self) -> impl Iterator<Item = &Planet> {
        self.planets.values()
    }

    pub fn moons(&self) -> impl Iterator<Item = &Moon> {
        self.moons.values()
    }

    pub fn asteroid_belts(&self) -> impl Iterator<Item = &AsteroidBelt> {
        self.asteroid_belts.values()
    }

    pub fn comet_clouds(&self) -> impl Iterator<Item = &CometCloud> {
        self.comet_clouds.values()
    }

    pub fn debris_disks(&self) -> impl Iterator<Item = &DebrisDisk> {
        self.debris_disks.values()
    }

    pub fn dynamics_mode(&self, system_id: StarSystemId) -> Option<DynamicsMode> {
        self.star_systems.get(&system_id).map(|system| {
            if self.system_requires_nbody(system) {
                DynamicsMode::Interacting
            } else {
                DynamicsMode::Keplerian
            }
        })
    }

    pub fn validate(&self) -> Result<(), ObjectError> {
        self.gas_composition
            .validate()
            .map_err(|error| ObjectError(error.to_string()))?;
        for (galaxy_id, galaxy) in &self.galaxies {
            if *galaxy_id != galaxy.id {
                return Err(ObjectError(
                    "galaxy registry key does not match its ID".to_owned(),
                ));
            }
            ensure_unique(&galaxy.star_system_ids, "galaxy star-system relationship")?;
            for system_id in &galaxy.star_system_ids {
                let system = self.star_systems.get(system_id).ok_or_else(|| {
                    ObjectError(format!("galaxy {galaxy_id} references missing {system_id}"))
                })?;
                if system.galaxy_id != *galaxy_id {
                    return Err(ObjectError(format!(
                        "{system_id} does not point back to {galaxy_id}"
                    )));
                }
            }
        }

        for (system_id, system) in &self.star_systems {
            if *system_id != system.id {
                return Err(ObjectError(
                    "star-system registry key does not match its ID".to_owned(),
                ));
            }
            ensure_unique(&system.planet_ids, "star-system planet relationship")?;
            let galaxy = self.galaxies.get(&system.galaxy_id).ok_or_else(|| {
                ObjectError(format!(
                    "{system_id} references missing galaxy {}",
                    system.galaxy_id
                ))
            })?;
            if !galaxy.star_system_ids.contains(system_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {system_id}",
                    system.galaxy_id
                )));
            }
            let star = self.stars.get(&system.star_id).ok_or_else(|| {
                ObjectError(format!(
                    "{system_id} references missing star {}",
                    system.star_id
                ))
            })?;
            if star.system_id != *system_id {
                return Err(ObjectError(format!(
                    "{} does not point back to {system_id}",
                    system.star_id
                )));
            }
            ensure_unique(
                &system.companion_star_ids,
                "star-system companion relationship",
            )?;
            for companion_id in &system.companion_star_ids {
                let companion = self.stars.get(companion_id).ok_or_else(|| {
                    ObjectError(format!(
                        "{system_id} references missing companion {companion_id}"
                    ))
                })?;
                if companion.system_id != *system_id || *companion_id == system.star_id {
                    return Err(ObjectError(format!(
                        "{companion_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for planet_id in &system.planet_ids {
                let planet = self.planets.get(planet_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {planet_id}"))
                })?;
                if planet.system_id != *system_id || planet.star_id != system.star_id {
                    return Err(ObjectError(format!(
                        "{planet_id} has inconsistent parent relationships"
                    )));
                }
            }
            ensure_unique(&system.moon_ids, "star-system moon relationship")?;
            for moon_id in &system.moon_ids {
                let moon = self.moons.get(moon_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {moon_id}"))
                })?;
                if moon.system_id != *system_id || !system.planet_ids.contains(&moon.planet_id) {
                    return Err(ObjectError(format!(
                        "{moon_id} has inconsistent parent relationships"
                    )));
                }
            }
            ensure_unique(
                &system.asteroid_belt_ids,
                "star-system asteroid-belt relationship",
            )?;
            ensure_unique(
                &system.comet_cloud_ids,
                "star-system comet-cloud relationship",
            )?;
            ensure_unique(
                &system.debris_disk_ids,
                "star-system debris-disk relationship",
            )?;
            for belt_id in &system.asteroid_belt_ids {
                let belt = self.asteroid_belts.get(belt_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {belt_id}"))
                })?;
                if belt.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{belt_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for cloud_id in &system.comet_cloud_ids {
                let cloud = self.comet_clouds.get(cloud_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {cloud_id}"))
                })?;
                if cloud.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{cloud_id} has an inconsistent parent relationship"
                    )));
                }
            }
            for disk_id in &system.debris_disk_ids {
                let disk = self.debris_disks.get(disk_id).ok_or_else(|| {
                    ObjectError(format!("{system_id} references missing {disk_id}"))
                })?;
                if disk.system_id != *system_id {
                    return Err(ObjectError(format!(
                        "{disk_id} has an inconsistent parent relationship"
                    )));
                }
            }
        }

        for (star_id, star) in &self.stars {
            if *star_id != star.id {
                return Err(ObjectError(
                    "star registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&star.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{star_id} references missing system {}",
                    star.system_id
                ))
            })?;
            if system.star_id != *star_id && !system.companion_star_ids.contains(star_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {star_id}",
                    star.system_id
                )));
            }
        }

        for (planet_id, planet) in &self.planets {
            if *planet_id != planet.id {
                return Err(ObjectError(
                    "planet registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&planet.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{planet_id} references missing system {}",
                    planet.system_id
                ))
            })?;
            if !system.planet_ids.contains(planet_id) || system.star_id != planet.star_id {
                return Err(ObjectError(format!(
                    "{} has inconsistent child relationship for {planet_id}",
                    planet.system_id
                )));
            }
        }
        for (moon_id, moon) in &self.moons {
            if *moon_id != moon.id {
                return Err(ObjectError(
                    "moon registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&moon.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{moon_id} references missing system {}",
                    moon.system_id
                ))
            })?;
            let planet = self.planets.get(&moon.planet_id).ok_or_else(|| {
                ObjectError(format!(
                    "{moon_id} references missing planet {}",
                    moon.planet_id
                ))
            })?;
            if !system.moon_ids.contains(moon_id)
                || planet.system_id != moon.system_id
                || planet.star_id != system.star_id
            {
                return Err(ObjectError(format!(
                    "{moon_id} has inconsistent child relationship"
                )));
            }
        }
        for (belt_id, belt) in &self.asteroid_belts {
            if *belt_id != belt.id {
                return Err(ObjectError(
                    "asteroid-belt registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&belt.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{belt_id} references missing system {}",
                    belt.system_id
                ))
            })?;
            if !system.asteroid_belt_ids.contains(belt_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {belt_id}",
                    belt.system_id
                )));
            }
            planetary::validate_small_body_population(
                belt.inner_radius,
                belt.outer_radius,
                belt.estimated_body_count,
                "asteroid belt",
            )?;
        }
        for (cloud_id, cloud) in &self.comet_clouds {
            if *cloud_id != cloud.id {
                return Err(ObjectError(
                    "comet-cloud registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&cloud.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{cloud_id} references missing system {}",
                    cloud.system_id
                ))
            })?;
            if !system.comet_cloud_ids.contains(cloud_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {cloud_id}",
                    cloud.system_id
                )));
            }
            planetary::validate_small_body_population(
                cloud.inner_radius,
                cloud.outer_radius,
                cloud.estimated_body_count,
                "comet cloud",
            )?;
        }
        for (disk_id, disk) in &self.debris_disks {
            if *disk_id != disk.id {
                return Err(ObjectError(
                    "debris-disk registry key does not match its ID".to_owned(),
                ));
            }
            let system = self.star_systems.get(&disk.system_id).ok_or_else(|| {
                ObjectError(format!(
                    "{disk_id} references missing system {}",
                    disk.system_id
                ))
            })?;
            if !system.debris_disk_ids.contains(disk_id) {
                return Err(ObjectError(format!(
                    "{} does not list child {disk_id}",
                    disk.system_id
                )));
            }
            planetary::validate_small_body_population(
                disk.inner_radius,
                disk.outer_radius,
                disk.estimated_body_count,
                "debris disk",
            )?;
        }
        Ok(())
    }

    pub(crate) fn advance_stars_to(
        &mut self,
        elapsed: SimDuration,
        tick_duration: SimDuration,
    ) -> Result<Vec<StellarTransition>, TimeError> {
        let mut scheduled = BinaryHeap::new();
        let stars: Vec<_> = self.stars.values().collect();
        let schedule_star = |star: &&Star| -> Result<Vec<_>, TimeError> {
            let target_age = star
                .epoch_age
                .checked_add(SimDuration::from_seconds(
                    elapsed.seconds().saturating_sub(star.formed_at.seconds()),
                ))
                .ok_or(TimeError::Overflow)?;
            let mut items = Vec::new();
            for (threshold, stage) in
                stellar::stellar_schedule(star.initial_mass, star.metallicity)?
            {
                if star.age < threshold && threshold <= target_age && stage > star.stage {
                    let event_elapsed =
                        SimDuration::from_seconds(star.formed_at.seconds().saturating_add(
                            threshold.seconds().saturating_sub(star.epoch_age.seconds()),
                        ));
                    let tick = tick_at_or_after(event_elapsed, tick_duration)?;
                    items.push((tick, threshold.seconds(), star.id, stage));
                }
            }
            Ok(items)
        };
        let schedules: Result<Vec<_>, TimeError> = if stars.len() >= PARALLEL_BODY_THRESHOLD {
            stars.par_iter().map(schedule_star).collect()
        } else {
            stars.iter().map(schedule_star).collect()
        };
        for item in schedules?.into_iter().flatten() {
            scheduled.push(Reverse(item));
        }

        let mut transitions = Vec::new();
        while let Some(Reverse((tick, _threshold, star_id, stage))) = scheduled.pop() {
            let star = self.stars.get_mut(&star_id).expect("scheduled star exists");
            let from = star.stage;
            let previous_mass = star.current_mass;
            star.stage = stage;
            let (mass, luminosity, temperature) =
                stellar::stellar_properties(star.initial_mass, star.metallicity, stage)
                    .map_err(|_| TimeError::Overflow)?;
            star.current_mass = mass;
            star.luminosity = luminosity;
            star.temperature = temperature;
            let natal_kick_meters_per_second = stellar::natal_kick_speed(star, stage, tick);
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
        }
        let trigger_ticks: BTreeMap<_, _> = transitions
            .iter()
            .filter(|transition| transition.interaction.is_none())
            .map(|transition| (transition.star_id, transition.tick))
            .collect();
        transitions.extend(self.resolve_binary_interactions(&trigger_ticks, elapsed)?);
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

struct GeneratedSystem {
    system: StarSystem,
    stars: Vec<Star>,
    planets: Vec<Planet>,
    moons: Vec<Moon>,
    asteroid_belts: Vec<AsteroidBelt>,
    comet_clouds: Vec<CometCloud>,
    debris_disks: Vec<DebrisDisk>,
}

fn generate_system(
    seed: u64,
    config: GenerationConfig,
    galaxy: GalaxyCatalogEntry,
    system_index: u32,
) -> Result<GeneratedSystem, ObjectError> {
    let galaxy_id = galaxy.id;
    let system_key = u64::from(system_index);
    let system_id = StarSystemId::new(derive_object_seed(seed, system_key, "star-system-id-v1"));
    let star_id = StarId::new(derive_object_seed(seed, system_key, "star-id-v1"));
    let mut system_rng = StableRng::new(derive_object_seed(
        seed,
        system_id.value(),
        "star-system-properties-v1",
    ));
    let population = generate_system_population(&mut system_rng, galaxy)?;
    let [gx, gy, gz] = galaxy.position.meters();
    let [rx, ry, rz] = population.position.meters();
    let position = Position::from_meters(gx + rx, gy + ry, gz + rz)?;
    let velocity = population.velocity;
    let mut primary = generate_star(
        seed,
        star_id,
        system_id,
        population.epoch_age,
        population.metallicity,
        config,
    )?;
    primary.position = position;
    primary.velocity = velocity;
    let primary_mass = primary.initial_mass;
    let current_primary_mass = primary.current_mass;
    let current_primary_luminosity = primary.luminosity;
    let primary_stage = primary.stage;
    let mut stars = vec![primary];

    let binary_probability =
        (0.20 + 0.15 * primary_mass.solar_masses().ln().max(0.0)).clamp(0.15, 0.70);
    let mut companion_star_ids = Vec::new();
    let mut inner_companion_axis_au = None;
    if system_rng.next_f64() < binary_probability {
        let triple_probability =
            (0.05 + 0.04 * primary_mass.solar_masses().sqrt()).clamp(0.05, 0.25);
        let companion_count = if system_rng.next_f64() < triple_probability {
            2
        } else {
            1
        };
        for companion_index in 0..companion_count {
            let namespace = format!("companion-star-id-v1-{companion_index}");
            let companion_id = StarId::new(derive_object_seed(seed, system_key, &namespace));
            let mass_ratio = random_range(&mut system_rng, 0.10, 1.0);
            let companion_mass_solar = (primary_mass.solar_masses() * mass_ratio)
                .clamp(config.star_mass_solar_min, config.star_mass_solar_max);
            let mut companion = generate_star_with_mass(
                companion_id,
                system_id,
                population.epoch_age,
                population.metallicity,
                companion_mass_solar,
            )?;
            let semi_major_axis_au = if let Some(inner_axis) = inner_companion_axis_au {
                inner_axis * random_range(&mut system_rng, 10.0, 30.0)
            } else {
                (random_range(&mut system_rng, 0.05_f64.ln(), 50.0_f64.ln())).exp()
            };
            inner_companion_axis_au.get_or_insert(semi_major_axis_au);
            let orbit = OrbitalElements::new(
                Length::from_astronomical_units(semi_major_axis_au)?,
                0.6 * system_rng.next_f64().powi(2),
                Angle::from_degrees(random_signed(&mut system_rng) * 30.0)?,
                Angle::from_radians(std::f64::consts::TAU * system_rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * system_rng.next_f64())?,
                Angle::from_radians(std::f64::consts::TAU * system_rng.next_f64())?,
                SimDuration::ZERO,
            )?;
            let interior_mass: f64 = stars.iter().map(|star| star.initial_mass.kilograms()).sum();
            let orbital_mass =
                Mass::from_kilograms(interior_mass + companion.initial_mass.kilograms())?;
            let (companion_position, companion_velocity) =
                propagate_kepler(orbit, orbital_mass, SimDuration::ZERO, position, velocity)?;
            companion.position = companion_position;
            companion.velocity = companion_velocity;
            companion.orbital_elements = Some(orbit);
            stars.push(companion);
            companion_star_ids.push(companion_id);
        }
    }

    let planet_count = if primary_stage == StellarStage::Disrupted {
        0
    } else {
        planetary::generate_planet_count(
            &mut system_rng,
            config,
            primary_mass,
            population.metallicity,
        )
    };
    let orbital_radii = planetary::stable_orbital_radii(&mut system_rng, planet_count, config);
    let max_planet_orbit_au = orbital_radii.iter().copied().fold(0.0_f64, f64::max);
    let mut planet_ids = Vec::with_capacity(usize::from(planet_count));
    let mut planets = Vec::with_capacity(usize::from(planet_count));
    let mut moon_ids = Vec::new();
    let mut moons = Vec::new();
    for (planet_index, orbit_au) in (0..planet_count).zip(orbital_radii) {
        let planet_key = (system_key << 8) | u64::from(planet_index);
        let planet_id = PlanetId::new(derive_object_seed(seed, planet_key, "planet-id-v1"));
        let planet = planetary::generate_planet(
            seed,
            planet_id,
            system_id,
            star_id,
            position,
            velocity,
            current_primary_mass,
            current_primary_luminosity,
            orbit_au,
            population.metallicity,
            config,
        )?;
        planet_ids.push(planet_id);
        let moon_count = planetary::generate_moon_count(&mut system_rng, planet.mass);
        for moon_index in 0..moon_count {
            let moon_key =
                (system_key << 16) | (u64::from(planet_index) << 8) | u64::from(moon_index);
            let moon_id = MoonId::new(derive_object_seed(seed, moon_key, "moon-id-v1"));
            if let Some(moon) = planetary::generate_moon(
                moon_id,
                planet.id,
                system_id,
                planet.position,
                planet.velocity,
                planet.mass,
                current_primary_mass,
                orbit_au,
                &mut system_rng,
            )? {
                moons.push(moon);
                moon_ids.push(moon_id);
            }
        }
        planets.push(planet);
    }

    let mut asteroid_belt_ids = Vec::new();
    let mut asteroid_belts = Vec::new();
    if primary_stage != StellarStage::Disrupted && system_rng.next_f64() < 0.70 {
        let id = AsteroidBeltId::new(derive_object_seed(seed, system_key, "asteroid-belt-id-v1"));
        asteroid_belts.push(planetary::generate_asteroid_belt(
            id,
            system_id,
            &mut system_rng,
            max_planet_orbit_au,
        )?);
        asteroid_belt_ids.push(id);
    }
    let mut comet_cloud_ids = Vec::new();
    let mut comet_clouds = Vec::new();
    if primary_stage != StellarStage::Disrupted && system_rng.next_f64() < 0.55 {
        let id = CometCloudId::new(derive_object_seed(seed, system_key, "comet-cloud-id-v1"));
        comet_clouds.push(planetary::generate_comet_cloud(
            id,
            system_id,
            &mut system_rng,
        )?);
        comet_cloud_ids.push(id);
    }
    let mut debris_disk_ids = Vec::new();
    let mut debris_disks = Vec::new();
    if primary_stage != StellarStage::Disrupted
        && population.epoch_age.seconds() <= 200_000_000 * SECONDS_PER_YEAR
        && system_rng.next_f64() < 0.80
    {
        let id = DebrisDiskId::new(derive_object_seed(seed, system_key, "debris-disk-id-v1"));
        debris_disks.push(planetary::generate_debris_disk(
            id,
            system_id,
            &mut system_rng,
        )?);
        debris_disk_ids.push(id);
    }

    Ok(GeneratedSystem {
        system: StarSystem {
            id: system_id,
            galaxy_id,
            position,
            velocity,
            star_id,
            companion_star_ids,
            planet_ids,
            moon_ids,
            asteroid_belt_ids,
            comet_cloud_ids,
            debris_disk_ids,
            component: population.component,
        },
        stars,
        planets,
        moons,
        asteroid_belts,
        comet_clouds,
        debris_disks,
    })
}

struct SystemPopulation {
    component: GalacticComponent,
    position: Position,
    velocity: Velocity,
    epoch_age: SimDuration,
    metallicity: f64,
}

fn generate_system_population(
    rng: &mut StableRng,
    galaxy: GalaxyCatalogEntry,
) -> Result<SystemPopulation, ObjectError> {
    let component = match (galaxy.morphology, rng.next_f64()) {
        (GalaxyMorphology::Elliptical, value) if value < 0.85 => GalacticComponent::Bulge,
        (GalaxyMorphology::Dwarf | GalaxyMorphology::Irregular, value) if value < 0.70 => {
            GalacticComponent::Disk
        }
        (_, value) if value < 0.80 => GalacticComponent::Disk,
        (_, value) if value < 0.95 => GalacticComponent::Bulge,
        _ => GalacticComponent::Halo,
    };
    let galaxy_radius = galaxy.radius_parsecs * METERS_PER_PARSEC;
    let (position, velocity, age_years, base_metallicity) = match component {
        GalacticComponent::Disk => {
            let radius = galaxy_radius * rng.next_f64().sqrt();
            let angle = std::f64::consts::TAU * rng.next_f64();
            let (sin_angle, cos_angle) = angle.sin_cos();
            let scale_height = galaxy_radius * 0.015;
            let vertical =
                random_signed(rng) * scale_height * (-(1.0 - rng.next_f64()).max(1e-12).ln());
            let circular_speed = 220_000.0;
            (
                Position::from_meters(radius * cos_angle, radius * sin_angle, vertical)?,
                Velocity::from_meters_per_second(
                    -circular_speed * sin_angle + random_signed(rng) * 20_000.0,
                    circular_speed * cos_angle + random_signed(rng) * 20_000.0,
                    random_signed(rng) * 10_000.0,
                )?,
                random_range(rng, 0.0, 10_000_000_000.0),
                0.024 * (-1.5 * radius / galaxy_radius).exp(),
            )
        }
        GalacticComponent::Bulge => {
            let radius = galaxy_radius * 0.15 * rng.next_f64().cbrt();
            let direction = random_unit_vector(rng);
            (
                Position::from_meters(
                    radius * direction[0],
                    radius * direction[1],
                    radius * direction[2],
                )?,
                Velocity::from_meters_per_second(
                    random_signed(rng) * 140_000.0,
                    random_signed(rng) * 140_000.0,
                    random_signed(rng) * 140_000.0,
                )?,
                random_range(rng, 8_000_000_000.0, 12_000_000_000.0),
                0.022,
            )
        }
        GalacticComponent::Halo => {
            let radius = galaxy_radius * rng.next_f64().cbrt();
            let direction = random_unit_vector(rng);
            (
                Position::from_meters(
                    radius * direction[0],
                    radius * direction[1],
                    radius * direction[2],
                )?,
                Velocity::from_meters_per_second(
                    random_signed(rng) * MAX_SYSTEM_SPEED_METERS_PER_SECOND,
                    random_signed(rng) * MAX_SYSTEM_SPEED_METERS_PER_SECOND,
                    random_signed(rng) * MAX_SYSTEM_SPEED_METERS_PER_SECOND,
                )?,
                random_range(rng, 10_000_000_000.0, 13_000_000_000.0),
                0.0015,
            )
        }
    };
    let (age_years, base_metallicity) = match galaxy.morphology {
        GalaxyMorphology::Spiral => (age_years, base_metallicity),
        GalaxyMorphology::Elliptical => (random_range(rng, 8.0e9, 13.0e9), 0.022),
        GalaxyMorphology::Irregular => (random_range(rng, 0.0, 8.0e9), 0.008),
        GalaxyMorphology::Dwarf => (random_range(rng, 0.0, 12.0e9), 0.003),
    };
    let age_metallicity_factor = (1.0 - 0.35 * age_years / 13_000_000_000.0).max(0.2);
    let metallicity = (base_metallicity * age_metallicity_factor).clamp(0.0001, 0.04);
    let epoch_age = SimDuration::from_years(age_years.round() as u128)
        .map_err(|error| ObjectError(error.to_string()))?;
    Ok(SystemPopulation {
        component,
        position,
        velocity,
        epoch_age,
        metallicity,
    })
}

fn generate_star(
    seed: u64,
    id: StarId,
    system_id: StarSystemId,
    epoch_age: SimDuration,
    metallicity: f64,
    config: GenerationConfig,
) -> Result<Star, ObjectError> {
    let mut rng = StableRng::new(derive_object_seed(seed, id.value(), "star-properties-v1"));
    let mass_solar = sample_stellar_mass(
        &mut rng,
        config.star_mass_solar_min,
        config.star_mass_solar_max,
    );
    generate_star_with_mass(id, system_id, epoch_age, metallicity, mass_solar)
}

fn generate_star_with_mass(
    id: StarId,
    system_id: StarSystemId,
    epoch_age: SimDuration,
    metallicity: f64,
    mass_solar: f64,
) -> Result<Star, ObjectError> {
    let initial_mass = Mass::from_solar_masses(mass_solar)?;
    let stage = stellar::stellar_stage_at(initial_mass, metallicity, epoch_age)?;
    let (current_mass, luminosity, temperature) =
        stellar::stellar_properties(initial_mass, metallicity, stage)?;
    Ok(Star {
        id,
        system_id,
        initial_mass,
        current_mass,
        epoch_age,
        formed_at: SimDuration::ZERO,
        age: epoch_age,
        metallicity,
        luminosity,
        temperature,
        stage,
        position: Position::default(),
        velocity: Velocity::default(),
        orbital_elements: None,
    })
}

fn generate_star_with_birth(
    id: StarId,
    system_id: StarSystemId,
    formed_at: SimDuration,
    metallicity: f64,
    mass_solar: f64,
) -> Result<Star, ObjectError> {
    let initial_mass = Mass::from_solar_masses(mass_solar)?;
    let stage = StellarStage::Protostar;
    let (current_mass, luminosity, temperature) =
        stellar::stellar_properties(initial_mass, metallicity, stage)?;
    Ok(Star {
        id,
        system_id,
        initial_mass,
        current_mass,
        epoch_age: SimDuration::ZERO,
        formed_at,
        age: SimDuration::ZERO,
        metallicity,
        luminosity,
        temperature,
        stage,
        position: Position::default(),
        velocity: Velocity::default(),
        orbital_elements: None,
    })
}

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

fn ensure_unique<T>(ids: &[T], relationship: &str) -> Result<(), ObjectError>
where
    T: Ord + Copy + fmt::Display,
{
    let mut seen = BTreeSet::new();
    for id in ids {
        if !seen.insert(*id) {
            return Err(ObjectError(format!("duplicate {id} in {relationship}")));
        }
    }
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generation_is_deterministic_and_relationships_are_valid() {
        let first = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
        let second = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
        assert_eq!(first, second);
        first.validate().unwrap();
        assert_eq!(first.counts().galaxies, 1);
        assert_eq!(first.counts().star_systems, 4);
        assert!(first.counts().stars >= 4);
        assert!(first.counts().stars <= 12);
    }

    #[test]
    fn generation_adds_moons_and_aggregate_small_body_populations() {
        let mut config = GenerationConfig::DEFAULT;
        config.star_systems = 128;
        let objects = ObjectRegistry::generate(42, config).unwrap();
        let counts = objects.counts();
        assert!(counts.moons > 0);
        assert!(counts.asteroid_belts > 0);
        assert!(counts.comet_clouds > 0);
        assert!(counts.debris_disks > 0);
        assert!(objects.planets().all(|planet| {
            planet.radius().meters() > 0.0 && planet.equilibrium_temperature().kelvin() > 0.0
        }));
        for moon in objects.moons() {
            let planet = objects.planet(moon.planet_id()).unwrap();
            assert_eq!(planet.system_id(), moon.system_id());
            assert!(moon.radius().meters() > 0.0);
            assert!(moon.orbital_radius().meters() > 0.0);
            assert!(
                objects
                    .star_system(moon.system_id())
                    .unwrap()
                    .moon_ids()
                    .contains(&moon.id())
            );
        }
        for belt in objects.asteroid_belts() {
            assert!(belt.inner_radius() < belt.outer_radius());
            assert!(belt.estimated_body_count() > 0);
        }
        objects.validate().unwrap();
    }

    #[test]
    fn different_seeds_generate_different_objects() {
        let first = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
        let second = ObjectRegistry::generate(43, GenerationConfig::DEFAULT).unwrap();
        assert_ne!(first, second);
    }

    #[test]
    fn generation_supports_three_hundred_solar_mass_stars() {
        let mut config = GenerationConfig::DEFAULT;
        config.star_systems = 1;
        config.planets_per_system_min = 0;
        config.planets_per_system_max = 0;
        config.star_mass_solar_min = 300.0;
        config.star_mass_solar_max = 300.0;
        let objects = ObjectRegistry::generate(42, config).unwrap();
        assert!(
            objects
                .stars()
                .all(|star| (star.initial_mass().solar_masses() - 300.0).abs() < 1e-9)
        );
    }

    #[test]
    fn very_massive_low_metallicity_stars_cover_pair_instability_fates() {
        let mass = |solar| Mass::from_solar_masses(solar).unwrap();
        assert_eq!(
            stellar::stellar_fate(mass(100.0), 0.0001),
            stellar::StellarFate::PulsationalPairInstability
        );
        assert_eq!(
            stellar::stellar_fate(mass(180.0), 0.0001),
            stellar::StellarFate::PairInstability
        );
        assert_eq!(
            stellar::stellar_fate(mass(300.0), 0.0001),
            stellar::StellarFate::DirectCollapse
        );
        assert!(
            stellar::stellar_schedule(mass(100.0), 0.0001)
                .unwrap()
                .iter()
                .any(|(_, stage)| *stage == StellarStage::PulsationalPairInstability)
        );
        assert!(
            stellar::stellar_schedule(mass(180.0), 0.0001)
                .unwrap()
                .iter()
                .any(|(_, stage)| *stage == StellarStage::Disrupted)
        );
    }

    #[test]
    fn solar_mass_star_has_a_complete_low_mass_lifecycle() {
        let initial_mass = Mass::from_solar_masses(1.0).unwrap();
        let schedule = stellar::stellar_schedule(initial_mass, 0.014).unwrap();
        let stages: Vec<_> = schedule.iter().map(|(_, stage)| *stage).collect();
        assert_eq!(
            stages,
            vec![
                StellarStage::PreMainSequence,
                StellarStage::MainSequence,
                StellarStage::Subgiant,
                StellarStage::Giant,
                StellarStage::CoreHeliumBurning,
                StellarStage::AsymptoticGiantBranch,
                StellarStage::PlanetaryNebula,
                StellarStage::WhiteDwarf,
            ]
        );
        assert!(
            schedule
                .windows(2)
                .all(|entries| entries[0].0 < entries[1].0)
        );
        assert_eq!(
            stellar::stellar_stage_at(initial_mass, 0.014, SimDuration::ZERO).unwrap(),
            StellarStage::Protostar
        );
        assert_eq!(
            stellar::stellar_stage_at(initial_mass, 0.014, schedule.last().unwrap().0).unwrap(),
            StellarStage::WhiteDwarf
        );
    }

    #[test]
    fn massive_solar_metallicity_star_explodes_after_a_wolf_rayet_phase() {
        let schedule =
            stellar::stellar_schedule(Mass::from_solar_masses(30.0).unwrap(), 0.014).unwrap();
        let stages: Vec<_> = schedule.iter().map(|(_, stage)| *stage).collect();
        assert!(stages.contains(&StellarStage::Supergiant));
        assert!(stages.contains(&StellarStage::WolfRayet));
        assert!(stages.contains(&StellarStage::CoreCollapseSupernova));
        assert_eq!(stages.last(), Some(&StellarStage::BlackHole));
        assert!(
            schedule
                .windows(2)
                .all(|entries| entries[0].0 < entries[1].0)
        );
    }

    #[test]
    fn gas_reservoir_forms_a_persisted_protostar_system() {
        let mut config = GenerationConfig::DEFAULT;
        config.star_systems = 1;
        config.star_mass_solar_min = 0.1;
        config.star_mass_solar_max = 0.2;
        let mut objects = ObjectRegistry::generate(42, config).unwrap();
        let initial_systems = objects.counts().star_systems;
        let initial_gas = objects.gas_reservoir_solar_masses();
        let initial_composition = objects.gas_composition();
        let ten_myr = SimDuration::from_years(10_000_000).unwrap();
        let events = objects
            .form_stars_between(
                42,
                config,
                SimDuration::ZERO,
                ten_myr,
                SimDuration::from_years(1_000).unwrap(),
            )
            .unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0].interaction,
            Some(StellarInteractionKind::StarFormation)
        );
        assert_eq!(objects.counts().star_systems, initial_systems + 1);
        assert!(objects.gas_reservoir_solar_masses() < initial_gas);
        assert_eq!(objects.gas_composition(), initial_composition);
        assert_eq!(
            objects.star(events[0].star_id).unwrap().formed_at(),
            ten_myr
        );
        let mut encoded = String::new();
        objects.encode(&mut encoded);
        assert_eq!(
            ObjectRegistry::decode(encoded.lines().skip(1)).unwrap(),
            objects
        );
    }

    #[test]
    fn galaxy_morphology_sets_the_initial_gas_metallicity() {
        let expected = [
            (GalaxyMorphology::Spiral, 0.014),
            (GalaxyMorphology::Elliptical, 0.020),
            (GalaxyMorphology::Irregular, 0.006),
            (GalaxyMorphology::Dwarf, 0.002),
        ];
        let mut config = GenerationConfig::DEFAULT;
        config.star_systems = 1;
        for (index, (morphology, metallicity)) in expected.into_iter().enumerate() {
            let galaxy = GalaxyCatalogEntry {
                id: GalaxyId::new(index as u64 + 1),
                position: Position::default(),
                morphology,
                radius_parsecs: config.galaxy_radius_parsecs,
            };
            let objects = ObjectRegistry::generate_in_galaxy(index as u64, config, galaxy).unwrap();
            assert!((objects.gas_composition().metallicity() - metallicity).abs() < 1.0e-12);
        }
    }

    #[test]
    fn roche_overflow_transfers_mass_and_contact_merges() {
        let mut selected = (0..100)
            .find_map(|seed| {
                let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
                objects
                    .star_systems()
                    .find(|system| !system.companion_star_ids().is_empty())
                    .map(|system| (objects.clone(), system.id()))
            })
            .expect("a deterministic binary system");
        let system = selected.0.star_system(selected.1).unwrap().clone();
        let primary_id = system.star_id;
        let companion_id = system.companion_star_ids[0];
        for id in [primary_id, companion_id] {
            let star = selected.0.stars.get_mut(&id).unwrap();
            star.current_mass = Mass::from_solar_masses(1.0).unwrap();
            star.initial_mass = star.current_mass;
            star.stage = StellarStage::MainSequence;
        }
        selected.0.stars.get_mut(&primary_id).unwrap().stage = StellarStage::Giant;
        selected
            .0
            .stars
            .get_mut(&companion_id)
            .unwrap()
            .orbital_elements = Some(test_orbit(1.0));
        let before = selected.0.star(primary_id).unwrap().current_mass();
        let outcomes = selected
            .0
            .resolve_binary_interactions(&BTreeMap::from([(primary_id, 1)]), SimDuration::ZERO)
            .unwrap();
        assert_eq!(
            outcomes[0].interaction,
            Some(StellarInteractionKind::MassTransfer)
        );
        assert!(selected.0.star(primary_id).unwrap().current_mass() < before);

        selected
            .0
            .stars
            .get_mut(&companion_id)
            .unwrap()
            .orbital_elements = Some(test_orbit(0.1));
        let outcomes = selected
            .0
            .resolve_binary_interactions(&BTreeMap::from([(primary_id, 2)]), SimDuration::ZERO)
            .unwrap();
        assert_eq!(
            outcomes[0].interaction,
            Some(StellarInteractionKind::StellarMerger)
        );
        assert!(selected.0.star(companion_id).is_none());
        selected.0.validate().unwrap();
    }

    #[test]
    fn expanded_star_physically_consumes_intersecting_planet() {
        let mut config = GenerationConfig::DEFAULT;
        config.star_systems = 1;
        config.planets_per_system_min = 1;
        config.planets_per_system_max = 1;
        config.orbit_au_min = 0.1;
        config.orbit_au_max = 0.1;
        let mut objects = ObjectRegistry::generate(42, config).unwrap();
        let star_id = objects.star_systems().next().unwrap().star_id();
        {
            let star = objects.stars.get_mut(&star_id).unwrap();
            star.current_mass = Mass::from_solar_masses(1.0).unwrap();
            star.stage = StellarStage::Giant;
        }
        let impacts = objects.apply_direct_planet_collisions(star_id, StellarStage::Giant);
        assert_eq!(impacts.len(), 1);
        assert_eq!(impacts[0].kind, PlanetImpactKind::Destroyed);
        assert_eq!(
            objects.planets().next().unwrap().status(),
            PlanetStatus::Destroyed
        );
    }

    fn test_orbit(semi_major_axis_au: f64) -> OrbitalElements {
        OrbitalElements::new(
            Length::from_astronomical_units(semi_major_axis_au).unwrap(),
            0.0,
            Angle::from_radians(0.0).unwrap(),
            Angle::from_radians(0.0).unwrap(),
            Angle::from_radians(0.0).unwrap(),
            Angle::from_radians(0.0).unwrap(),
            SimDuration::ZERO,
        )
        .unwrap()
    }

    #[test]
    fn generated_companions_have_hierarchical_orbits() {
        let objects = (0..100)
            .find_map(|seed| {
                let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
                let has_companion = objects
                    .star_systems()
                    .any(|system| !system.companion_star_ids().is_empty());
                has_companion.then_some(objects)
            })
            .expect("at least one deterministic seed generates a companion");
        for system in objects
            .star_systems()
            .filter(|system| !system.companion_star_ids().is_empty())
        {
            let mut axes: Vec<_> = system
                .companion_star_ids()
                .iter()
                .map(|id| {
                    objects
                        .star(*id)
                        .unwrap()
                        .orbital_elements()
                        .unwrap()
                        .semi_major_axis()
                        .astronomical_units()
                })
                .collect();
            axes.sort_by(f64::total_cmp);
            if axes.len() == 2 {
                assert!(axes[1] / axes[0] >= 10.0);
            }
        }
    }

    #[test]
    fn interacting_systems_use_the_selective_nbody_path_for_short_windows() {
        let elapsed = SimDuration::from_years(1).unwrap();
        let mut selected = None;
        for seed in 0..200 {
            let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
            let interacting_system = objects.star_systems().find_map(|system| {
                (objects.dynamics_mode(system.id()) == Some(DynamicsMode::Interacting))
                    .then_some(system.id())
            });
            if let Some(system_id) = interacting_system {
                selected = Some((objects, system_id));
                break;
            }
        }
        let (mut objects, system_id) = selected.expect("an interacting generated system");
        assert!(
            objects
                .integrate_system_from_epoch(system_id, elapsed)
                .unwrap()
        );
    }

    #[test]
    fn object_records_round_trip_exactly() {
        let objects = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
        let mut encoded = String::new();
        objects.encode(&mut encoded);
        let decoded = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
        assert_eq!(decoded, objects);
    }

    #[test]
    fn orbital_catch_up_changes_positions_without_changing_elements() {
        let mut objects = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
        let before: Vec<_> = objects
            .planets()
            .map(|planet| (planet.id(), planet.position(), planet.orbital_elements()))
            .collect();
        objects
            .advance_orbits_to(SimDuration::from_years(1_000_000).unwrap())
            .unwrap();
        assert!(
            before
                .iter()
                .any(|(id, position, _)| { objects.planet(*id).unwrap().position() != *position })
        );
        assert!(before.iter().all(|(id, _, elements)| {
            objects.planet(*id).unwrap().orbital_elements() == *elements
        }));
    }

    #[test]
    fn parallel_generation_and_advancement_match_single_threaded_results() {
        let mut config = GenerationConfig::DEFAULT;
        config.star_systems = 512;
        let run = || {
            let mut objects = ObjectRegistry::generate(42, config).unwrap();
            let elapsed = SimDuration::from_years(1_000_000_000).unwrap();
            let tick = SimDuration::from_years(1).unwrap();
            let transitions = objects.advance_stars_to(elapsed, tick).unwrap();
            objects.advance_orbits_to(elapsed).unwrap();
            (objects, transitions)
        };
        let single_thread = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(run);
        let four_threads = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(run);
        assert_eq!(single_thread, four_threads);
    }
}
