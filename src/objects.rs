use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use std::error::Error;
use std::fmt;
use std::fmt::Write as _;
use std::str::FromStr;

use rayon::prelude::*;

use crate::chemistry::GasComposition;
use crate::config::GenerationConfig;
use crate::id::{
    AsteroidBeltId, CometCloudId, DebrisDiskId, GalaxyId, IdParseError, MoonId, PlanetId, StarId,
    StarSystemId,
};
use crate::orbital::{
    GRAVITATIONAL_CONSTANT, GravityBody, OrbitalElements, OrbitalError, leapfrog_step,
    propagate_kepler,
};
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration, TimeError};
use crate::units::{
    Angle, Length, Luminosity, METERS_PER_ASTRONOMICAL_UNIT, METERS_PER_PARSEC, Mass, Position,
    Temperature, UnitError, Velocity,
};

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
            validate_small_body_population(
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
            validate_small_body_population(
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
            validate_small_body_population(
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
            for (threshold, stage) in stellar_schedule(star.initial_mass, star.metallicity)? {
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
                stellar_properties(star.initial_mass, star.metallicity, stage)
                    .map_err(|_| TimeError::Overflow)?;
            star.current_mass = mass;
            star.luminosity = luminosity;
            star.temperature = temperature;
            let natal_kick_meters_per_second = natal_kick_speed(star, stage, tick);
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

    fn refresh_planet_temperatures(&mut self) -> Result<(), TimeError> {
        let mut updates = Vec::new();
        for planet in self.planets.values() {
            let Some(star) = self.stars.get(&planet.star_id) else {
                continue;
            };
            let temperature = Temperature::from_kelvin(equilibrium_temperature_kelvin(
                star.luminosity,
                planet.orbital_radius.astronomical_units(),
            ))
            .map_err(|_| TimeError::Overflow)?;
            updates.push((planet.id, temperature));
        }
        for (planet_id, temperature) in updates {
            self.planets
                .get_mut(&planet_id)
                .expect("refreshed planet exists")
                .equilibrium_temperature = temperature;
        }
        Ok(())
    }

    fn apply_direct_planet_collisions(
        &mut self,
        star_id: StarId,
        stage: StellarStage,
    ) -> Vec<PlanetImpact> {
        let Some(star) = self.stars.get(&star_id) else {
            return Vec::new();
        };
        let collision_radius_au = stellar_radius_solar(star.current_mass, stage)
            * SOLAR_RADIUS_METERS
            / METERS_PER_ASTRONOMICAL_UNIT;
        let mut impacts = Vec::new();
        for planet in self
            .planets
            .values_mut()
            .filter(|planet| planet.star_id == star_id && planet.status != PlanetStatus::Destroyed)
        {
            let periapsis_au = planet
                .orbital_elements
                .semi_major_axis()
                .astronomical_units()
                * (1.0 - planet.orbital_elements.eccentricity());
            if periapsis_au <= collision_radius_au {
                planet.status = PlanetStatus::Destroyed;
                impacts.push(PlanetImpact {
                    planet_id: planet.id,
                    kind: PlanetImpactKind::Destroyed,
                });
            }
        }
        impacts
    }

    fn resolve_binary_interactions(
        &mut self,
        trigger_ticks: &BTreeMap<StarId, u64>,
        elapsed: SimDuration,
    ) -> Result<Vec<StellarTransition>, TimeError> {
        let mut candidates = Vec::new();
        for system in self.star_systems.values() {
            for companion_id in &system.companion_star_ids {
                let Some(companion) = self.stars.get(companion_id) else {
                    continue;
                };
                let Some(orbit) = companion.orbital_elements else {
                    continue;
                };
                let tick = trigger_ticks
                    .get(&system.star_id)
                    .or_else(|| trigger_ticks.get(companion_id));
                if let Some(tick) = tick {
                    candidates.push((system.id, system.star_id, *companion_id, orbit, *tick));
                }
            }
        }
        let mut outcomes = Vec::new();
        for (system_id, primary_id, companion_id, orbit, tick) in candidates {
            let (Some(primary), Some(companion)) =
                (self.stars.get(&primary_id), self.stars.get(&companion_id))
            else {
                continue;
            };
            let separation_au =
                orbit.semi_major_axis().astronomical_units() * (1.0 - orbit.eccentricity());
            let primary_radius_au = stellar_radius_solar(primary.current_mass, primary.stage)
                * SOLAR_RADIUS_METERS
                / METERS_PER_ASTRONOMICAL_UNIT;
            let companion_radius_au = stellar_radius_solar(companion.current_mass, companion.stage)
                * SOLAR_RADIUS_METERS
                / METERS_PER_ASTRONOMICAL_UNIT;
            let primary_lobe = roche_lobe_radius(
                separation_au,
                primary.current_mass.solar_masses() / companion.current_mass.solar_masses(),
            );
            let companion_lobe = roche_lobe_radius(
                separation_au,
                companion.current_mass.solar_masses() / primary.current_mass.solar_masses(),
            );
            let primary_overflow = primary_radius_au > primary_lobe;
            let companion_overflow = companion_radius_au > companion_lobe;
            if !primary_overflow
                && !companion_overflow
                && primary_radius_au + companion_radius_au < separation_au
            {
                continue;
            }
            let mass_ratio = primary
                .current_mass
                .solar_masses()
                .max(companion.current_mass.solar_masses())
                / primary
                    .current_mass
                    .solar_masses()
                    .min(companion.current_mass.solar_masses());
            if primary_radius_au + companion_radius_au >= separation_au
                || (primary_overflow || companion_overflow) && mass_ratio > 3.0
            {
                let previous_stage = primary.stage;
                let merged_mass = 0.9
                    * (primary.current_mass.solar_masses() + companion.current_mass.solar_masses());
                let merged_mass =
                    Mass::from_solar_masses(merged_mass).map_err(|_| TimeError::Overflow)?;
                let primary = self
                    .stars
                    .get_mut(&primary_id)
                    .expect("merger primary exists");
                primary.initial_mass = merged_mass;
                primary.current_mass = merged_mass;
                primary.epoch_age = SimDuration::ZERO;
                primary.formed_at = elapsed;
                primary.age = SimDuration::ZERO;
                primary.stage = StellarStage::MainSequence;
                let (_, luminosity, temperature) = stellar_properties(
                    merged_mass,
                    primary.metallicity,
                    StellarStage::MainSequence,
                )
                .map_err(|_| TimeError::Overflow)?;
                primary.luminosity = luminosity;
                primary.temperature = temperature;
                self.stars.remove(&companion_id);
                if let Some(system) = self.star_systems.get_mut(&system_id) {
                    system.companion_star_ids.retain(|id| *id != companion_id);
                }
                outcomes.push(StellarTransition {
                    tick,
                    star_id: primary_id,
                    from: previous_stage,
                    to: StellarStage::MainSequence,
                    impacts: Vec::new(),
                    natal_kick_meters_per_second: 0,
                    interaction: Some(StellarInteractionKind::StellarMerger),
                    related_star_id: Some(companion_id),
                });
                continue;
            }
            let (donor_id, accretor_id) = if primary_overflow {
                (primary_id, companion_id)
            } else {
                (companion_id, primary_id)
            };
            let donor = self.stars.get(&donor_id).expect("binary donor exists");
            let donor_stage = donor.stage;
            let donor_old_solar = donor.current_mass.solar_masses();
            let accretor_old_solar = self
                .stars
                .get(&accretor_id)
                .expect("binary accretor exists")
                .current_mass
                .solar_masses();
            let transferred_solar = 0.05 * donor_old_solar;
            let donor_mass = Mass::from_solar_masses(donor_old_solar - transferred_solar)
                .map_err(|_| TimeError::Overflow)?;
            let accretor_mass =
                Mass::from_solar_masses(accretor_old_solar + 0.70 * transferred_solar)
                    .map_err(|_| TimeError::Overflow)?;
            self.stars
                .get_mut(&donor_id)
                .expect("binary donor exists")
                .current_mass = donor_mass;
            self.stars
                .get_mut(&accretor_id)
                .expect("binary accretor exists")
                .current_mass = accretor_mass;
            let old_product = donor_old_solar * accretor_old_solar;
            let new_product = donor_mass.solar_masses() * accretor_mass.solar_masses();
            let axis_scale = (old_product / new_product).powi(2).clamp(0.25, 4.0);
            if let Some(companion) = self.stars.get_mut(&companion_id)
                && let Some(existing) = companion.orbital_elements
            {
                companion.orbital_elements = Some(
                    OrbitalElements::new(
                        Length::from_meters(existing.semi_major_axis().meters() * axis_scale)
                            .map_err(|_| TimeError::Overflow)?,
                        existing.eccentricity() * 0.5,
                        existing.inclination(),
                        existing.longitude_of_ascending_node(),
                        existing.argument_of_periapsis(),
                        existing.mean_anomaly_at_epoch(),
                        existing.epoch(),
                    )
                    .map_err(|_| TimeError::Overflow)?,
                );
            }
            outcomes.push(StellarTransition {
                tick,
                star_id: donor_id,
                from: donor_stage,
                to: donor_stage,
                impacts: Vec::new(),
                natal_kick_meters_per_second: 0,
                interaction: Some(StellarInteractionKind::MassTransfer),
                related_star_id: Some(accretor_id),
            });
        }
        Ok(outcomes)
    }

    fn apply_stellar_impacts(
        &mut self,
        star_id: StarId,
        previous_mass: Mass,
        current_mass: Mass,
        stage: StellarStage,
    ) -> Vec<PlanetImpact> {
        let retained_mass = current_mass.kilograms() / previous_mass.kilograms();
        let mut impacts = Vec::new();
        for planet in self
            .planets
            .values_mut()
            .filter(|planet| planet.star_id == star_id)
        {
            if matches!(
                planet.status,
                PlanetStatus::Unbound | PlanetStatus::Destroyed
            ) {
                continue;
            }
            let periapsis_au = planet
                .orbital_elements
                .semi_major_axis()
                .astronomical_units()
                * (1.0 - planet.orbital_elements.eccentricity());
            let outcome = match stage {
                StellarStage::AsymptoticGiantBranch if periapsis_au <= 2.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::AsymptoticGiantBranch if periapsis_au <= 5.0 => {
                    Some(PlanetImpactKind::Irradiated)
                }
                StellarStage::Giant if periapsis_au <= 0.5 => Some(PlanetImpactKind::Destroyed),
                StellarStage::Giant if periapsis_au <= 2.0 => Some(PlanetImpactKind::Irradiated),
                StellarStage::Supergiant if periapsis_au <= 3.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::Supergiant if periapsis_au <= 10.0 => {
                    Some(PlanetImpactKind::Irradiated)
                }
                StellarStage::PulsationalPairInstability if periapsis_au <= 20.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::PulsationalPairInstability => Some(PlanetImpactKind::Irradiated),
                StellarStage::CoreCollapseSupernova if periapsis_au <= 20.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::CoreCollapseSupernova => Some(PlanetImpactKind::Irradiated),
                StellarStage::Disrupted if periapsis_au <= 20.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::Disrupted => Some(PlanetImpactKind::Unbound),
                StellarStage::NeutronStar | StellarStage::BlackHole if periapsis_au <= 5.0 => {
                    Some(PlanetImpactKind::Destroyed)
                }
                StellarStage::PlanetaryNebula
                | StellarStage::WhiteDwarf
                | StellarStage::NeutronStar
                | StellarStage::BlackHole
                    if retained_mass < 0.5 =>
                {
                    Some(PlanetImpactKind::Unbound)
                }
                StellarStage::NeutronStar | StellarStage::BlackHole => {
                    Some(PlanetImpactKind::Irradiated)
                }
                _ => None,
            };
            let Some(kind) = outcome else {
                continue;
            };
            let new_status = match kind {
                PlanetImpactKind::Irradiated => PlanetStatus::Irradiated,
                PlanetImpactKind::Unbound => PlanetStatus::Unbound,
                PlanetImpactKind::Destroyed => PlanetStatus::Destroyed,
            };
            if planet.status == new_status {
                continue;
            }
            planet.status = new_status;
            impacts.push(PlanetImpact {
                planet_id: planet.id,
                kind,
            });
        }
        impacts
    }

    pub(crate) fn advance_orbits_to(&mut self, elapsed: SimDuration) -> Result<(), ObjectError> {
        let system_ids: Vec<_> = self.star_systems.keys().copied().collect();
        let mut integrated_systems = BTreeSet::new();
        for system_id in system_ids {
            if self.integrate_system_from_epoch(system_id, elapsed)? {
                integrated_systems.insert(system_id);
            }
        }
        let companion_inputs: Vec<_> = self
            .stars
            .values()
            .filter_map(|star| {
                (!integrated_systems.contains(&star.system_id))
                    .then_some(star.orbital_elements)
                    .flatten()
                    .map(|orbit| (star.id, star.system_id, orbit))
            })
            .collect();
        let propagate_companion =
            |(star_id, system_id, orbit): &(StarId, StarSystemId, OrbitalElements)| {
                let system = &self.star_systems[system_id];
                let primary = &self.stars[&system.star_id];
                let companion = &self.stars[star_id];
                let interior_mass: f64 = system
                    .companion_star_ids
                    .iter()
                    .filter_map(|id| self.stars.get(id))
                    .filter_map(|star| star.orbital_elements.map(|candidate| (star, candidate)))
                    .filter(|(_, candidate)| {
                        candidate.semi_major_axis().meters() < orbit.semi_major_axis().meters()
                    })
                    .map(|(star, _)| star.current_mass.kilograms())
                    .sum();
                let orbital_mass = Mass::from_kilograms(
                    primary.current_mass.kilograms()
                        + interior_mass
                        + companion.current_mass.kilograms(),
                )?;
                if orbital_mass.kilograms() == 0.0 {
                    return Ok((*star_id, companion.position, companion.velocity));
                }
                let (position, velocity) = propagate_kepler(
                    *orbit,
                    orbital_mass,
                    elapsed,
                    system.position,
                    system.velocity,
                )?;
                Ok((*star_id, position, velocity))
            };
        let companion_updates: Result<Vec<_>, ObjectError> =
            if companion_inputs.len() >= PARALLEL_BODY_THRESHOLD {
                companion_inputs
                    .par_iter()
                    .map(propagate_companion)
                    .collect()
            } else {
                companion_inputs.iter().map(propagate_companion).collect()
            };
        for (star_id, position, velocity) in companion_updates? {
            let star = self
                .stars
                .get_mut(&star_id)
                .expect("propagated star exists");
            star.position = position;
            star.velocity = velocity;
        }

        let planets: Vec<_> = self
            .planets
            .values()
            .filter(|planet| {
                !integrated_systems.contains(&planet.system_id)
                    && !matches!(
                        planet.status,
                        PlanetStatus::Unbound | PlanetStatus::Destroyed
                    )
            })
            .collect();
        let propagate = |planet: &&Planet| -> Result<_, ObjectError> {
            let system = &self.star_systems[&planet.system_id];
            let star = &self.stars[&planet.star_id];
            let (position, velocity) = propagate_kepler(
                planet.orbital_elements,
                star.current_mass,
                elapsed,
                system.position,
                system.velocity,
            )?;
            Ok((planet.id, position, velocity))
        };
        let updates: Result<Vec<_>, ObjectError> = if planets.len() >= PARALLEL_BODY_THRESHOLD {
            planets.par_iter().map(propagate).collect()
        } else {
            planets.iter().map(propagate).collect()
        };
        for (planet_id, position, velocity) in updates? {
            let planet = self
                .planets
                .get_mut(&planet_id)
                .expect("propagated planet exists");
            planet.position = position;
            planet.velocity = velocity;
        }
        let moons: Vec<_> = self
            .moons
            .values()
            .filter_map(|moon| {
                let planet = self.planets.get(&moon.planet_id)?;
                (!matches!(
                    planet.status,
                    PlanetStatus::Unbound | PlanetStatus::Destroyed
                ))
                .then_some((moon, planet.position, planet.velocity, planet.mass))
            })
            .collect();
        let propagate_moon = |(moon, planet_position, planet_velocity, planet_mass): &(
            &Moon,
            Position,
            Velocity,
            Mass,
        )| {
            let (position, velocity) = propagate_kepler(
                moon.orbital_elements,
                *planet_mass,
                elapsed,
                *planet_position,
                *planet_velocity,
            )?;
            Ok((moon.id, position, velocity))
        };
        let moon_updates: Result<Vec<_>, ObjectError> = if moons.len() >= PARALLEL_BODY_THRESHOLD {
            moons.par_iter().map(propagate_moon).collect()
        } else {
            moons.iter().map(propagate_moon).collect()
        };
        for (moon_id, position, velocity) in moon_updates? {
            let moon = self
                .moons
                .get_mut(&moon_id)
                .expect("propagated moon exists");
            moon.position = position;
            moon.velocity = velocity;
        }
        Ok(())
    }

    fn system_requires_nbody(&self, system: &StarSystem) -> bool {
        let mut companion_orbits: Vec<_> = system
            .companion_star_ids
            .iter()
            .filter_map(|id| self.stars[id].orbital_elements)
            .collect();
        companion_orbits.sort_by(|left, right| {
            left.semi_major_axis()
                .meters()
                .total_cmp(&right.semi_major_axis().meters())
        });
        let max_planet_apoapsis = system
            .planet_ids
            .iter()
            .filter_map(|id| self.planets.get(id))
            .filter(|planet| {
                !matches!(
                    planet.status,
                    PlanetStatus::Destroyed | PlanetStatus::Unbound
                )
            })
            .map(|planet| {
                planet.orbital_elements.semi_major_axis().meters()
                    * (1.0 + planet.orbital_elements.eccentricity())
            })
            .fold(0.0_f64, f64::max);
        if max_planet_apoapsis > 0.0
            && companion_orbits.iter().any(|orbit| {
                orbit.semi_major_axis().meters() * (1.0 - orbit.eccentricity())
                    < 5.0 * max_planet_apoapsis
            })
        {
            return true;
        }
        companion_orbits.windows(2).any(|orbits| {
            let inner_apoapsis =
                orbits[0].semi_major_axis().meters() * (1.0 + orbits[0].eccentricity());
            let outer_periapsis =
                orbits[1].semi_major_axis().meters() * (1.0 - orbits[1].eccentricity());
            outer_periapsis < 5.0 * inner_apoapsis
        })
    }

    fn integrate_system_from_epoch(
        &mut self,
        system_id: StarSystemId,
        elapsed: SimDuration,
    ) -> Result<bool, ObjectError> {
        const STEPS_PER_SHORTEST_ORBIT: f64 = 32.0;
        const MAX_NBODY_STEPS: usize = 4_096;
        let system = &self.star_systems[&system_id];
        if !self.system_requires_nbody(system) || elapsed == SimDuration::ZERO {
            return Ok(false);
        }
        let primary = &self.stars[&system.star_id];
        let active_planets: Vec<_> = system
            .planet_ids
            .iter()
            .filter_map(|id| self.planets.get(id))
            .filter(|planet| {
                !matches!(
                    planet.status,
                    PlanetStatus::Destroyed | PlanetStatus::Unbound
                )
            })
            .collect();
        let companions: Vec<_> = system
            .companion_star_ids
            .iter()
            .filter_map(|id| self.stars.get(id))
            .filter(|star| star.orbital_elements.is_some())
            .collect();
        let primary_mass = primary.current_mass.kilograms();
        if primary_mass <= 0.0 {
            return Ok(false);
        }
        let shortest_axis = companions
            .iter()
            .filter_map(|star| star.orbital_elements)
            .map(|orbit| orbit.semi_major_axis().meters())
            .chain(
                active_planets
                    .iter()
                    .map(|planet| planet.orbital_elements.semi_major_axis().meters()),
            )
            .fold(f64::INFINITY, f64::min);
        if !shortest_axis.is_finite() {
            return Ok(false);
        }
        let shortest_period = std::f64::consts::TAU
            * (shortest_axis.powi(3) / (GRAVITATIONAL_CONSTANT * primary_mass)).sqrt();
        let desired_step = shortest_period / STEPS_PER_SHORTEST_ORBIT;
        let elapsed_seconds = elapsed.seconds() as f64;
        let step_count = (elapsed_seconds / desired_step).ceil() as usize;
        if step_count == 0 || step_count > MAX_NBODY_STEPS {
            return Ok(false);
        }
        let step_seconds = elapsed_seconds / step_count as f64;
        let mut keys = Vec::with_capacity(1 + companions.len() + active_planets.len());
        let mut bodies = Vec::with_capacity(keys.capacity());
        keys.push(NBodyObject::Star(primary.id));
        bodies.push(GravityBody {
            mass: primary.current_mass,
            position: Position::default(),
            velocity: Velocity::default(),
        });
        for companion in companions {
            let orbit = companion
                .orbital_elements
                .expect("filtered companion orbit");
            let interior_mass: f64 = self.star_systems[&system_id]
                .companion_star_ids
                .iter()
                .filter_map(|id| self.stars.get(id))
                .filter_map(|star| star.orbital_elements.map(|candidate| (star, candidate)))
                .filter(|(_, candidate)| {
                    candidate.semi_major_axis().meters() < orbit.semi_major_axis().meters()
                })
                .map(|(star, _)| star.current_mass.kilograms())
                .sum();
            let total_mass = Mass::from_kilograms(
                primary.current_mass.kilograms()
                    + interior_mass
                    + companion.current_mass.kilograms(),
            )?;
            let (position, velocity) = propagate_kepler(
                orbit,
                total_mass,
                SimDuration::ZERO,
                Position::default(),
                Velocity::default(),
            )?;
            keys.push(NBodyObject::Star(companion.id));
            bodies.push(GravityBody {
                mass: companion.current_mass,
                position,
                velocity,
            });
        }
        for planet in active_planets {
            let (position, velocity) = propagate_kepler(
                planet.orbital_elements,
                primary.current_mass,
                SimDuration::ZERO,
                Position::default(),
                Velocity::default(),
            )?;
            keys.push(NBodyObject::Planet(planet.id));
            bodies.push(GravityBody {
                mass: planet.mass,
                position,
                velocity,
            });
        }
        for _ in 0..step_count {
            leapfrog_step(&mut bodies, step_seconds)?;
        }
        let [sx, sy, sz] = system.position.meters();
        let [svx, svy, svz] = system.velocity.meters_per_second();
        for (key, body) in keys.into_iter().zip(bodies) {
            let [x, y, z] = body.position.meters();
            let [vx, vy, vz] = body.velocity.meters_per_second();
            let position = Position::from_meters(sx + x, sy + y, sz + z)?;
            let velocity = Velocity::from_meters_per_second(svx + vx, svy + vy, svz + vz)?;
            match key {
                NBodyObject::Star(id) => {
                    let star = self.stars.get_mut(&id).expect("integrated star exists");
                    star.position = position;
                    star.velocity = velocity;
                }
                NBodyObject::Planet(id) => {
                    let planet = self.planets.get_mut(&id).expect("integrated planet exists");
                    planet.position = position;
                    planet.velocity = velocity;
                }
            }
        }
        Ok(true)
    }

    pub(crate) fn encode(&self, output: &mut String) {
        output.push_str("[objects]\n");
        writeln!(
            output,
            "M|{}|{}|{}|{}|{}|{}|{}|{}",
            float_bits(self.gas_reservoir_solar_masses),
            self.formation_sequence,
            float_bits(self.gas_composition.hydrogen_fraction()),
            float_bits(self.gas_composition.helium_fraction()),
            float_bits(self.gas_composition.carbon_fraction()),
            float_bits(self.gas_composition.oxygen_fraction()),
            float_bits(self.gas_composition.iron_fraction()),
            float_bits(self.gas_composition.other_metals_fraction()),
        )
        .expect("writing to a String cannot fail");
        for galaxy in self.galaxies.values() {
            let [x, y, z] = galaxy.position.meters();
            writeln!(
                output,
                "G|{}|{}|{}|{}|{}",
                galaxy.id.value(),
                float_bits(x),
                float_bits(y),
                float_bits(z),
                encode_ids(&galaxy.star_system_ids)
            )
            .expect("writing to a String cannot fail");
        }
        for system in self.star_systems.values() {
            let [px, py, pz] = system.position.meters();
            let [vx, vy, vz] = system.velocity.meters_per_second();
            writeln!(
                output,
                "Y|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                system.id.value(),
                system.galaxy_id.value(),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                system.star_id.value(),
                encode_ids(&system.planet_ids),
                system.component.code(),
                encode_ids(&system.companion_star_ids),
                encode_ids(&system.moon_ids),
                encode_ids(&system.asteroid_belt_ids),
                encode_ids(&system.comet_cloud_ids),
                encode_ids(&system.debris_disk_ids),
            )
            .expect("writing to a String cannot fail");
        }
        for star in self.stars.values() {
            let [px, py, pz] = star.position.meters();
            let [vx, vy, vz] = star.velocity.meters_per_second();
            writeln!(
                output,
                "S|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                star.id.value(),
                star.system_id.value(),
                float_bits(star.initial_mass.kilograms()),
                float_bits(star.current_mass.kilograms()),
                star.epoch_age.seconds(),
                star.age.seconds(),
                float_bits(star.metallicity),
                float_bits(star.luminosity.watts()),
                float_bits(star.temperature.kelvin()),
                star.stage.code(),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                encode_optional_orbit(star.orbital_elements),
                star.formed_at.seconds(),
            )
            .expect("writing to a String cannot fail");
        }
        for planet in self.planets.values() {
            let [px, py, pz] = planet.position.meters();
            let [vx, vy, vz] = planet.velocity.meters_per_second();
            writeln!(
                output,
                "P|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                planet.id.value(),
                planet.system_id.value(),
                planet.star_id.value(),
                float_bits(planet.mass.kilograms()),
                float_bits(planet.orbital_radius.meters()),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                float_bits(planet.orbital_elements.semi_major_axis().meters()),
                float_bits(planet.orbital_elements.eccentricity()),
                float_bits(planet.orbital_elements.inclination().radians()),
                float_bits(
                    planet
                        .orbital_elements
                        .longitude_of_ascending_node()
                        .radians()
                ),
                float_bits(planet.orbital_elements.argument_of_periapsis().radians()),
                float_bits(planet.orbital_elements.mean_anomaly_at_epoch().radians()),
                planet.orbital_elements.epoch().seconds(),
                planet.status.code(),
                float_bits(planet.radius.meters()),
                float_bits(planet.equilibrium_temperature.kelvin()),
                planet.atmosphere.code(),
            )
            .expect("writing to a String cannot fail");
        }
        for moon in self.moons.values() {
            let [px, py, pz] = moon.position.meters();
            let [vx, vy, vz] = moon.velocity.meters_per_second();
            writeln!(
                output,
                "N|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                moon.id.value(),
                moon.planet_id.value(),
                moon.system_id.value(),
                float_bits(moon.mass.kilograms()),
                float_bits(moon.radius.meters()),
                float_bits(moon.orbital_radius.meters()),
                float_bits(px),
                float_bits(py),
                float_bits(pz),
                float_bits(vx),
                float_bits(vy),
                float_bits(vz),
                float_bits(moon.orbital_elements.semi_major_axis().meters()),
                float_bits(moon.orbital_elements.eccentricity()),
                float_bits(moon.orbital_elements.inclination().radians()),
                float_bits(
                    moon.orbital_elements
                        .longitude_of_ascending_node()
                        .radians()
                ),
                float_bits(moon.orbital_elements.argument_of_periapsis().radians()),
                float_bits(moon.orbital_elements.mean_anomaly_at_epoch().radians()),
                moon.orbital_elements.epoch().seconds(),
            )
            .expect("writing to a String cannot fail");
        }
        for belt in self.asteroid_belts.values() {
            writeln!(
                output,
                "A|{}|{}|{}|{}|{}|{}",
                belt.id.value(),
                belt.system_id.value(),
                float_bits(belt.inner_radius.meters()),
                float_bits(belt.outer_radius.meters()),
                belt.estimated_body_count,
                float_bits(belt.total_mass.kilograms()),
            )
            .expect("writing to a String cannot fail");
        }
        for cloud in self.comet_clouds.values() {
            writeln!(
                output,
                "C|{}|{}|{}|{}|{}|{}",
                cloud.id.value(),
                cloud.system_id.value(),
                float_bits(cloud.inner_radius.meters()),
                float_bits(cloud.outer_radius.meters()),
                cloud.estimated_body_count,
                float_bits(cloud.total_mass.kilograms()),
            )
            .expect("writing to a String cannot fail");
        }
        for disk in self.debris_disks.values() {
            writeln!(
                output,
                "D|{}|{}|{}|{}|{}|{}",
                disk.id.value(),
                disk.system_id.value(),
                float_bits(disk.inner_radius.meters()),
                float_bits(disk.outer_radius.meters()),
                disk.estimated_body_count,
                float_bits(disk.total_mass.kilograms()),
            )
            .expect("writing to a String cannot fail");
        }
    }

    pub(crate) fn decode<'a>(lines: impl Iterator<Item = &'a str>) -> Result<Self, ObjectError> {
        let mut registry = Self::default();
        for line in lines {
            if line.is_empty() {
                continue;
            }
            let fields: Vec<_> = line.split('|').collect();
            match fields.first().copied() {
                Some("M") => {
                    if !matches!(fields.len(), 3 | 9) {
                        return Err(ObjectError(format!(
                            "matter-reservoir record has {} fields; expected 3 or 9",
                            fields.len()
                        )));
                    }
                    registry.gas_reservoir_solar_masses =
                        parse_float_bits(fields[1], "gas reservoir mass")?;
                    registry.formation_sequence = parse_number(fields[2], "formation sequence")?;
                    if fields.len() == 9 {
                        registry.gas_composition = GasComposition::new(
                            parse_float_bits(fields[3], "hydrogen fraction")?,
                            parse_float_bits(fields[4], "helium fraction")?,
                            parse_float_bits(fields[5], "carbon fraction")?,
                            parse_float_bits(fields[6], "oxygen fraction")?,
                            parse_float_bits(fields[7], "iron fraction")?,
                            parse_float_bits(fields[8], "other-metals fraction")?,
                        )
                        .map_err(|error| ObjectError(error.to_string()))?;
                    }
                    if !registry.gas_reservoir_solar_masses.is_finite()
                        || registry.gas_reservoir_solar_masses < 0.0
                    {
                        return Err(ObjectError(
                            "gas reservoir mass must be finite and nonnegative".to_owned(),
                        ));
                    }
                }
                Some("G") => {
                    require_field_count(&fields, 6, "galaxy")?;
                    let galaxy = Galaxy {
                        id: GalaxyId::new(parse_number(fields[1], "galaxy ID")?),
                        position: parse_position(&fields[2..5])?,
                        star_system_ids: parse_ids(fields[5], StarSystemId::new)?,
                    };
                    insert_unique(&mut registry.galaxies, galaxy.id, galaxy, "galaxy")?;
                }
                Some("Y") => {
                    if !matches!(fields.len(), 11..=13 | 17) {
                        return Err(ObjectError(format!(
                            "star-system record has {} fields; expected 11, 12, 13, or 17",
                            fields.len()
                        )));
                    }
                    let system = StarSystem {
                        id: StarSystemId::new(parse_number(fields[1], "star-system ID")?),
                        galaxy_id: GalaxyId::new(parse_number(fields[2], "galaxy ID")?),
                        position: parse_position(&fields[3..6])?,
                        velocity: parse_velocity(&fields[6..9])?,
                        star_id: StarId::new(parse_number(fields[9], "star ID")?),
                        planet_ids: parse_ids(fields[10], PlanetId::new)?,
                        component: if fields.len() >= 12 {
                            GalacticComponent::parse_code(fields[11])?
                        } else {
                            GalacticComponent::Disk
                        },
                        companion_star_ids: if fields.len() == 13 || fields.len() == 17 {
                            parse_ids(fields[12], StarId::new)?
                        } else {
                            Vec::new()
                        },
                        moon_ids: if fields.len() == 17 {
                            parse_ids(fields[13], MoonId::new)?
                        } else {
                            Vec::new()
                        },
                        asteroid_belt_ids: if fields.len() == 17 {
                            parse_ids(fields[14], AsteroidBeltId::new)?
                        } else {
                            Vec::new()
                        },
                        comet_cloud_ids: if fields.len() == 17 {
                            parse_ids(fields[15], CometCloudId::new)?
                        } else {
                            Vec::new()
                        },
                        debris_disk_ids: if fields.len() == 17 {
                            parse_ids(fields[16], DebrisDiskId::new)?
                        } else {
                            Vec::new()
                        },
                    };
                    insert_unique(&mut registry.star_systems, system.id, system, "star system")?;
                }
                Some("S") => {
                    if !matches!(fields.len(), 8 | 9 | 11 | 18 | 19) {
                        return Err(ObjectError(format!(
                            "star record has {} fields; expected 8, 9, 11, 18, or 19",
                            fields.len()
                        )));
                    }
                    let is_current = matches!(fields.len(), 11 | 18 | 19);
                    let metallicity_index = if is_current { 7 } else { 5 };
                    let metallicity = parse_float_bits(fields[metallicity_index], "metallicity")?;
                    if !metallicity.is_finite() || metallicity < 0.0 {
                        return Err(ObjectError(
                            "star metallicity must be finite and nonnegative".to_owned(),
                        ));
                    }
                    let initial_mass =
                        Mass::from_kilograms(parse_float_bits(fields[3], "initial star mass")?)?;
                    let current_mass = if is_current {
                        Mass::from_kilograms(parse_float_bits(fields[4], "current star mass")?)?
                    } else {
                        initial_mass
                    };
                    let epoch_age = if is_current {
                        SimDuration::from_seconds(parse_number(fields[5], "epoch star age")?)
                    } else {
                        SimDuration::ZERO
                    };
                    let age_index = if is_current { 6 } else { 4 };
                    let age =
                        SimDuration::from_seconds(parse_number(fields[age_index], "star age")?);
                    let stage = if is_current {
                        StellarStage::parse_code(fields[10])?
                    } else if fields.len() == 9 {
                        StellarStage::parse_code(fields[8])?
                    } else {
                        stellar_stage_at(initial_mass, metallicity, age)?
                    };
                    let luminosity_index = if is_current { 8 } else { 6 };
                    let temperature_index = if is_current { 9 } else { 7 };
                    let star = Star {
                        id: StarId::new(parse_number(fields[1], "star ID")?),
                        system_id: StarSystemId::new(parse_number(fields[2], "system ID")?),
                        initial_mass,
                        current_mass,
                        epoch_age,
                        formed_at: if fields.len() == 19 {
                            SimDuration::from_seconds(parse_number(fields[18], "formation time")?)
                        } else {
                            SimDuration::ZERO
                        },
                        age,
                        metallicity,
                        luminosity: Luminosity::from_watts(parse_float_bits(
                            fields[luminosity_index],
                            "star luminosity",
                        )?)?,
                        temperature: Temperature::from_kelvin(parse_float_bits(
                            fields[temperature_index],
                            "star temperature",
                        )?)?,
                        stage,
                        position: if fields.len() >= 18 {
                            parse_position(&fields[11..14])?
                        } else {
                            Position::default()
                        },
                        velocity: if fields.len() >= 18 {
                            parse_velocity(&fields[14..17])?
                        } else {
                            Velocity::default()
                        },
                        orbital_elements: if fields.len() >= 18 {
                            parse_optional_orbit(fields[17])?
                        } else {
                            None
                        },
                    };
                    insert_unique(&mut registry.stars, star.id, star, "star")?;
                }
                Some("P") => {
                    if !matches!(fields.len(), 12 | 19 | 20 | 23) {
                        return Err(ObjectError(format!(
                            "planet record has {} fields; expected 12, 19, 20, or 23",
                            fields.len()
                        )));
                    }
                    let orbital_radius =
                        Length::from_meters(parse_float_bits(fields[5], "orbital radius")?)?;
                    let orbital_elements = if fields.len() >= 19 {
                        OrbitalElements::new(
                            Length::from_meters(parse_float_bits(fields[12], "semi-major axis")?)?,
                            parse_float_bits(fields[13], "eccentricity")?,
                            Angle::from_radians(parse_float_bits(fields[14], "inclination")?)?,
                            Angle::from_radians(parse_float_bits(fields[15], "ascending node")?)?,
                            Angle::from_radians(parse_float_bits(
                                fields[16],
                                "periapsis argument",
                            )?)?,
                            Angle::from_radians(parse_float_bits(fields[17], "mean anomaly")?)?,
                            SimDuration::from_seconds(parse_number(fields[18], "orbit epoch")?),
                        )?
                    } else {
                        OrbitalElements::new(
                            orbital_radius,
                            0.0,
                            Angle::default(),
                            Angle::default(),
                            Angle::default(),
                            Angle::default(),
                            SimDuration::ZERO,
                        )?
                    };
                    let mass = Mass::from_kilograms(parse_float_bits(fields[4], "planet mass")?)?;
                    let atmosphere = if fields.len() == 23 {
                        PlanetAtmosphere::parse_code(fields[22])?
                    } else {
                        PlanetAtmosphere::Thin
                    };
                    let planet = Planet {
                        id: PlanetId::new(parse_number(fields[1], "planet ID")?),
                        system_id: StarSystemId::new(parse_number(fields[2], "system ID")?),
                        star_id: StarId::new(parse_number(fields[3], "star ID")?),
                        mass,
                        radius: if fields.len() == 23 {
                            Length::from_meters(parse_float_bits(fields[20], "planet radius")?)?
                        } else {
                            Length::from_meters(planet_radius_meters(mass, atmosphere))?
                        },
                        atmosphere,
                        equilibrium_temperature: if fields.len() == 23 {
                            Temperature::from_kelvin(parse_float_bits(
                                fields[21],
                                "planet equilibrium temperature",
                            )?)?
                        } else {
                            Temperature::from_kelvin(288.0)?
                        },
                        orbital_radius,
                        position: parse_position(&fields[6..9])?,
                        velocity: parse_velocity(&fields[9..12])?,
                        orbital_elements,
                        status: if fields.len() == 20 || fields.len() == 23 {
                            PlanetStatus::parse_code(fields[19])?
                        } else {
                            PlanetStatus::Stable
                        },
                    };
                    insert_unique(&mut registry.planets, planet.id, planet, "planet")?;
                }
                Some("N") => {
                    if !matches!(fields.len(), 19 | 20) {
                        return Err(ObjectError(format!(
                            "moon record has {} fields; expected 19 or 20",
                            fields.len()
                        )));
                    }
                    let has_radius = fields.len() == 20;
                    let mass = Mass::from_kilograms(parse_float_bits(fields[4], "moon mass")?)?;
                    let radius = if has_radius {
                        Length::from_meters(parse_float_bits(fields[5], "moon radius")?)?
                    } else {
                        Length::from_meters(planet_radius_meters(mass, PlanetAtmosphere::None))?
                    };
                    let orbital_radius_index = if has_radius { 6 } else { 5 };
                    let position_start = if has_radius { 7 } else { 6 };
                    let orbit_start = if has_radius { 13 } else { 12 };
                    let orbital_radius = Length::from_meters(parse_float_bits(
                        fields[orbital_radius_index],
                        "moon orbital radius",
                    )?)?;
                    let orbital_elements = OrbitalElements::new(
                        Length::from_meters(parse_float_bits(
                            fields[orbit_start],
                            "moon semi-major axis",
                        )?)?,
                        parse_float_bits(fields[orbit_start + 1], "moon eccentricity")?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 2],
                            "moon inclination",
                        )?)?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 3],
                            "moon ascending node",
                        )?)?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 4],
                            "moon periapsis argument",
                        )?)?,
                        Angle::from_radians(parse_float_bits(
                            fields[orbit_start + 5],
                            "moon mean anomaly",
                        )?)?,
                        SimDuration::from_seconds(parse_number(
                            fields[orbit_start + 6],
                            "moon orbit epoch",
                        )?),
                    )?;
                    let moon = Moon {
                        id: MoonId::new(parse_number(fields[1], "moon ID")?),
                        planet_id: PlanetId::new(parse_number(fields[2], "moon planet ID")?),
                        system_id: StarSystemId::new(parse_number(fields[3], "moon system ID")?),
                        mass,
                        radius,
                        orbital_radius,
                        position: parse_position(&fields[position_start..position_start + 3])?,
                        velocity: parse_velocity(&fields[position_start + 3..position_start + 6])?,
                        orbital_elements,
                    };
                    insert_unique(&mut registry.moons, moon.id, moon, "moon")?;
                }
                Some("A") => {
                    require_field_count(&fields, 7, "asteroid belt")?;
                    let belt = AsteroidBelt {
                        id: AsteroidBeltId::new(parse_number(fields[1], "asteroid-belt ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "asteroid-belt system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "asteroid-belt inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "asteroid-belt outer radius",
                        )?)?,
                        estimated_body_count: parse_number(fields[5], "asteroid-belt body count")?,
                        total_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "asteroid-belt mass",
                        )?)?,
                    };
                    insert_unique(&mut registry.asteroid_belts, belt.id, belt, "asteroid belt")?;
                }
                Some("C") => {
                    require_field_count(&fields, 7, "comet cloud")?;
                    let cloud = CometCloud {
                        id: CometCloudId::new(parse_number(fields[1], "comet-cloud ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "comet-cloud system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "comet-cloud inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "comet-cloud outer radius",
                        )?)?,
                        estimated_body_count: parse_number(fields[5], "comet-cloud body count")?,
                        total_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "comet-cloud mass",
                        )?)?,
                    };
                    insert_unique(&mut registry.comet_clouds, cloud.id, cloud, "comet cloud")?;
                }
                Some("D") => {
                    require_field_count(&fields, 7, "debris disk")?;
                    let disk = DebrisDisk {
                        id: DebrisDiskId::new(parse_number(fields[1], "debris-disk ID")?),
                        system_id: StarSystemId::new(parse_number(
                            fields[2],
                            "debris-disk system ID",
                        )?),
                        inner_radius: Length::from_meters(parse_float_bits(
                            fields[3],
                            "debris-disk inner radius",
                        )?)?,
                        outer_radius: Length::from_meters(parse_float_bits(
                            fields[4],
                            "debris-disk outer radius",
                        )?)?,
                        estimated_body_count: parse_number(fields[5], "debris-disk body count")?,
                        total_mass: Mass::from_kilograms(parse_float_bits(
                            fields[6],
                            "debris-disk mass",
                        )?)?,
                    };
                    insert_unique(&mut registry.debris_disks, disk.id, disk, "debris disk")?;
                }
                Some(tag) => return Err(ObjectError(format!("unknown object record '{tag}'"))),
                None => return Err(ObjectError("empty object record".to_owned())),
            }
        }
        registry.validate()?;
        Ok(registry)
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

#[derive(Debug, Clone, Copy)]
enum NBodyObject {
    Star(StarId),
    Planet(PlanetId),
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
        generate_planet_count(
            &mut system_rng,
            config,
            primary_mass,
            population.metallicity,
        )
    };
    let orbital_radii = stable_orbital_radii(&mut system_rng, planet_count, config);
    let max_planet_orbit_au = orbital_radii.iter().copied().fold(0.0_f64, f64::max);
    let mut planet_ids = Vec::with_capacity(usize::from(planet_count));
    let mut planets = Vec::with_capacity(usize::from(planet_count));
    let mut moon_ids = Vec::new();
    let mut moons = Vec::new();
    for (planet_index, orbit_au) in (0..planet_count).zip(orbital_radii) {
        let planet_key = (system_key << 8) | u64::from(planet_index);
        let planet_id = PlanetId::new(derive_object_seed(seed, planet_key, "planet-id-v1"));
        let planet = generate_planet(
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
        let moon_count = generate_moon_count(&mut system_rng, planet.mass);
        for moon_index in 0..moon_count {
            let moon_key =
                (system_key << 16) | (u64::from(planet_index) << 8) | u64::from(moon_index);
            let moon_id = MoonId::new(derive_object_seed(seed, moon_key, "moon-id-v1"));
            if let Some(moon) = generate_moon(
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
        asteroid_belts.push(generate_asteroid_belt(
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
        comet_clouds.push(generate_comet_cloud(id, system_id, &mut system_rng)?);
        comet_cloud_ids.push(id);
    }
    let mut debris_disk_ids = Vec::new();
    let mut debris_disks = Vec::new();
    if primary_stage != StellarStage::Disrupted
        && population.epoch_age.seconds() <= 200_000_000 * SECONDS_PER_YEAR
        && system_rng.next_f64() < 0.80
    {
        let id = DebrisDiskId::new(derive_object_seed(seed, system_key, "debris-disk-id-v1"));
        debris_disks.push(generate_debris_disk(id, system_id, &mut system_rng)?);
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
    let stage = stellar_stage_at(initial_mass, metallicity, epoch_age)?;
    let (current_mass, luminosity, temperature) =
        stellar_properties(initial_mass, metallicity, stage)?;
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
        stellar_properties(initial_mass, metallicity, stage)?;
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

fn stellar_schedule(
    initial_mass: Mass,
    metallicity: f64,
) -> Result<Vec<(SimDuration, StellarStage)>, TimeError> {
    let solar_masses = initial_mass.solar_masses();
    let metallicity_factor = (metallicity / 0.014).clamp(0.05, 3.0).powf(0.12);
    let main_sequence_years = (10_000_000_000.0 * solar_masses.powf(-2.5) * metallicity_factor)
        .clamp(2_500_000.0, 200_000_000_000.0);
    let pre_main_sequence_years =
        (30_000_000.0 * solar_masses.powf(-2.5)).clamp(50_000.0, 2_000_000_000.0);
    let protostar_years = (pre_main_sequence_years * 0.15).clamp(10_000.0, 500_000.0);
    let main_sequence_end = pre_main_sequence_years + main_sequence_years;
    let mut schedule = vec![
        (
            years_to_duration(protostar_years)?,
            StellarStage::PreMainSequence,
        ),
        (
            years_to_duration(pre_main_sequence_years)?,
            StellarStage::MainSequence,
        ),
    ];
    match stellar_fate(initial_mass, metallicity) {
        StellarFate::WhiteDwarf => {
            schedule.push((
                years_to_duration(main_sequence_end)?,
                StellarStage::Subgiant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.05 * main_sequence_years)?,
                StellarStage::Giant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.10 * main_sequence_years)?,
                StellarStage::CoreHeliumBurning,
            ));
            let agb_start = main_sequence_end + 0.11 * main_sequence_years;
            let agb_years = (1_000_000.0 * solar_masses.powf(-1.5)).clamp(30_000.0, 2_000_000.0);
            schedule.push((
                years_to_duration(agb_start)?,
                StellarStage::AsymptoticGiantBranch,
            ));
            schedule.push((
                years_to_duration(agb_start + agb_years)?,
                StellarStage::PlanetaryNebula,
            ));
            schedule.push((
                years_to_duration(agb_start + agb_years + 30_000.0)?,
                StellarStage::WhiteDwarf,
            ));
        }
        StellarFate::NeutronStar | StellarFate::BlackHole => {
            schedule_massive_post_main_sequence(
                &mut schedule,
                main_sequence_end,
                main_sequence_years,
                solar_masses,
                metallicity,
            )?;
            let collapse = main_sequence_end + 0.05 * main_sequence_years;
            schedule.push((
                years_to_duration(collapse)?,
                StellarStage::CoreCollapseSupernova,
            ));
            let remnant = if matches!(
                stellar_fate(initial_mass, metallicity),
                StellarFate::NeutronStar
            ) {
                StellarStage::NeutronStar
            } else {
                StellarStage::BlackHole
            };
            schedule.push((years_to_duration(collapse + 1.0 / 365.25)?, remnant));
        }
        StellarFate::DirectCollapse => {
            schedule_massive_post_main_sequence(
                &mut schedule,
                main_sequence_end,
                main_sequence_years,
                solar_masses,
                metallicity,
            )?;
            schedule.push((
                years_to_duration(main_sequence_end + 0.05 * main_sequence_years)?,
                StellarStage::BlackHole,
            ));
        }
        StellarFate::PulsationalPairInstability => {
            schedule.push((
                years_to_duration(main_sequence_end)?,
                StellarStage::Supergiant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.03 * main_sequence_years)?,
                StellarStage::PulsationalPairInstability,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.08 * main_sequence_years)?,
                StellarStage::BlackHole,
            ));
        }
        StellarFate::PairInstability => {
            schedule.push((
                years_to_duration(main_sequence_end)?,
                StellarStage::Supergiant,
            ));
            schedule.push((
                years_to_duration(main_sequence_end + 0.05 * main_sequence_years)?,
                StellarStage::Disrupted,
            ));
        }
    }
    Ok(schedule)
}

fn schedule_massive_post_main_sequence(
    schedule: &mut Vec<(SimDuration, StellarStage)>,
    main_sequence_end: f64,
    main_sequence_years: f64,
    solar_masses: f64,
    metallicity: f64,
) -> Result<(), TimeError> {
    schedule.push((
        years_to_duration(main_sequence_end)?,
        StellarStage::Supergiant,
    ));
    if solar_masses >= 25.0 && metallicity >= 0.008 {
        schedule.push((
            years_to_duration(main_sequence_end + 0.025 * main_sequence_years)?,
            StellarStage::WolfRayet,
        ));
    }
    Ok(())
}

fn years_to_duration(years: f64) -> Result<SimDuration, TimeError> {
    let seconds = years * SECONDS_PER_YEAR as f64;
    if !seconds.is_finite() || seconds < 0.0 || seconds > u128::MAX as f64 {
        return Err(TimeError::Overflow);
    }
    Ok(SimDuration::from_seconds(seconds.round() as u128))
}

fn stellar_stage_at(
    initial_mass: Mass,
    metallicity: f64,
    age: SimDuration,
) -> Result<StellarStage, ObjectError> {
    let schedule = stellar_schedule(initial_mass, metallicity)
        .map_err(|error| ObjectError(error.to_string()))?;
    let mut stage = StellarStage::Protostar;
    for (threshold, scheduled_stage) in schedule {
        if age < threshold {
            break;
        }
        stage = scheduled_stage;
    }
    Ok(stage)
}

fn stellar_properties(
    initial_mass: Mass,
    metallicity: f64,
    stage: StellarStage,
) -> Result<(Mass, Luminosity, Temperature), ObjectError> {
    let solar_masses = initial_mass.solar_masses();
    let retained = wind_retained_fraction(solar_masses, metallicity);
    let pre_collapse_mass = solar_masses * retained;
    let main_sequence_luminosity = if solar_masses <= 20.0 {
        solar_masses.powf(3.5)
    } else {
        20.0_f64.powf(3.5) * (solar_masses / 20.0).powf(1.5)
    };
    let (current_mass_solar, luminosity_solar, temperature_kelvin) = match stage {
        StellarStage::Protostar => (
            solar_masses,
            (2.0 * solar_masses.powf(2.0)).max(0.001),
            (3_000.0 * solar_masses.powf(0.08)).clamp(2_500.0, 5_000.0),
        ),
        StellarStage::PreMainSequence => (
            solar_masses * 0.995,
            (1.5 * solar_masses.powf(2.5)).max(0.0005),
            (4_200.0 * solar_masses.powf(0.15)).clamp(2_700.0, 15_000.0),
        ),
        StellarStage::MainSequence => (
            solar_masses * (1.0 - 0.5 * (1.0 - retained)),
            main_sequence_luminosity.max(0.0001),
            (5_772.0 * solar_masses.powf(0.505)).min(60_000.0),
        ),
        StellarStage::Subgiant => (
            solar_masses * 0.97,
            (3.0 * main_sequence_luminosity).max(0.001),
            (5_200.0 * solar_masses.powf(0.08)).clamp(3_800.0, 8_000.0),
        ),
        StellarStage::Giant => (
            solar_masses * 0.85,
            1_000.0 * solar_masses.powf(2.0),
            3_800.0,
        ),
        StellarStage::CoreHeliumBurning => {
            (solar_masses * 0.80, 50.0 * solar_masses.powf(2.0), 5_000.0)
        }
        StellarStage::AsymptoticGiantBranch => (
            solar_masses * 0.65,
            10_000.0 * solar_masses.powf(1.2),
            3_200.0,
        ),
        StellarStage::PlanetaryNebula => (
            (0.109 * solar_masses + 0.394).clamp(0.45, 1.38),
            1_000.0,
            100_000.0,
        ),
        StellarStage::Supergiant => (
            pre_collapse_mass,
            30_000.0 * solar_masses.powf(1.5),
            4_500.0,
        ),
        StellarStage::WolfRayet => (
            pre_collapse_mass,
            50_000.0 * solar_masses.powf(1.3),
            80_000.0,
        ),
        StellarStage::PulsationalPairInstability => {
            (pre_collapse_mass.min(60.0), 500_000.0, 12_000.0)
        }
        StellarStage::CoreCollapseSupernova => (pre_collapse_mass, 1_000_000_000.0, 10_000.0),
        StellarStage::WhiteDwarf => (
            (0.109 * solar_masses + 0.394).clamp(0.45, 1.38),
            0.01,
            12_000.0,
        ),
        StellarStage::NeutronStar => (
            (1.25 + 0.05 * (solar_masses - 8.0)).clamp(1.2, 2.2),
            0.0001,
            600_000.0,
        ),
        StellarStage::BlackHole => {
            let remnant = match stellar_fate(initial_mass, metallicity) {
                StellarFate::PulsationalPairInstability => pre_collapse_mass.min(50.0) * 0.9,
                StellarFate::DirectCollapse => pre_collapse_mass * 0.9,
                _ => (3.0 + 0.25 * solar_masses).min(pre_collapse_mass).max(3.0),
            };
            (remnant, 0.0, 0.0)
        }
        StellarStage::Disrupted => (0.0, 0.0, 0.0),
    };
    Ok((
        Mass::from_solar_masses(current_mass_solar)?,
        Luminosity::from_solar_luminosities(luminosity_solar)?,
        Temperature::from_kelvin(temperature_kelvin)?,
    ))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StellarFate {
    WhiteDwarf,
    NeutronStar,
    BlackHole,
    PulsationalPairInstability,
    PairInstability,
    DirectCollapse,
}

fn wind_retained_fraction(initial_mass_solar: f64, metallicity: f64) -> f64 {
    if initial_mass_solar < 25.0 {
        return 1.0;
    }
    let metallicity_ratio = (metallicity / 0.014).clamp(0.001, 3.0);
    let wind_loss = (0.03
        + 0.45 * metallicity_ratio.powf(0.7) * (initial_mass_solar / 100.0).sqrt())
    .clamp(0.03, 0.80);
    1.0 - wind_loss
}

fn stellar_fate(initial_mass: Mass, metallicity: f64) -> StellarFate {
    let solar_masses = initial_mass.solar_masses();
    if solar_masses < 8.0 {
        return StellarFate::WhiteDwarf;
    }
    if solar_masses < 25.0 {
        return StellarFate::NeutronStar;
    }
    let helium_core = 0.55 * solar_masses * wind_retained_fraction(solar_masses, metallicity);
    if helium_core < 32.0 {
        StellarFate::BlackHole
    } else if helium_core < 64.0 {
        StellarFate::PulsationalPairInstability
    } else if helium_core < 135.0 {
        StellarFate::PairInstability
    } else {
        StellarFate::DirectCollapse
    }
}

fn natal_kick_speed(star: &Star, stage: StellarStage, tick: u64) -> u32 {
    let mut rng = StableRng::new(derive_object_seed(
        star.id.value(),
        tick,
        "natal-kick-speed-v1",
    ));
    match stage {
        StellarStage::NeutronStar => 100_000 + rng.next_u64_below(400_001).unwrap_or(0) as u32,
        StellarStage::BlackHole => match stellar_fate(star.initial_mass, star.metallicity) {
            StellarFate::DirectCollapse => 0,
            StellarFate::PulsationalPairInstability => {
                10_000 + rng.next_u64_below(40_001).unwrap_or(0) as u32
            }
            _ => 20_000 + rng.next_u64_below(80_001).unwrap_or(0) as u32,
        },
        _ => 0,
    }
}

fn stellar_radius_solar(mass: Mass, stage: StellarStage) -> f64 {
    let solar_masses = mass.solar_masses().max(0.01);
    match stage {
        StellarStage::Protostar => 3.0 * solar_masses.powf(0.5),
        StellarStage::PreMainSequence => 2.0 * solar_masses.powf(0.6),
        StellarStage::MainSequence => {
            solar_masses.powf(if solar_masses < 1.0 { 0.8 } else { 0.57 })
        }
        StellarStage::Subgiant => 3.0 * solar_masses.powf(0.6),
        StellarStage::Giant => 100.0 * solar_masses.powf(0.5),
        StellarStage::CoreHeliumBurning => 10.0 * solar_masses.powf(0.5),
        StellarStage::AsymptoticGiantBranch => 200.0 * solar_masses.powf(0.5),
        StellarStage::PlanetaryNebula => 0.02,
        StellarStage::Supergiant => 500.0 * (solar_masses / 10.0).sqrt(),
        StellarStage::WolfRayet => 3.0 * solar_masses.powf(0.4),
        StellarStage::PulsationalPairInstability => 100.0,
        StellarStage::CoreCollapseSupernova => 500.0,
        StellarStage::WhiteDwarf => 0.012,
        StellarStage::NeutronStar => 0.000_02,
        StellarStage::BlackHole | StellarStage::Disrupted => 0.0,
    }
}

fn roche_lobe_radius(separation: f64, mass_ratio: f64) -> f64 {
    let q_one_third = mass_ratio.max(1e-12).cbrt();
    let q_two_thirds = q_one_third * q_one_third;
    separation * 0.49 * q_two_thirds / (0.6 * q_two_thirds + (1.0 + q_one_third).ln())
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

#[allow(clippy::too_many_arguments)]
fn generate_planet(
    seed: u64,
    id: PlanetId,
    system_id: StarSystemId,
    star_id: StarId,
    planet_position: Position,
    planet_velocity: Velocity,
    star_mass: Mass,
    star_luminosity: Luminosity,
    orbit_au: f64,
    metallicity: f64,
    config: GenerationConfig,
) -> Result<Planet, ObjectError> {
    let mut rng = StableRng::new(derive_object_seed(seed, id.value(), "planet-properties-v1"));
    let metal_mass_bias = (metallicity / 0.014).clamp(0.25, 2.0);
    let mass_min = config.planet_mass_earth_min.ln();
    let mass_max = config.planet_mass_earth_max.ln();
    let mass_earth = (mass_min + rng.next_f64().powf(1.5) * (mass_max - mass_min))
        .exp()
        .mul_add(metal_mass_bias.sqrt(), 0.0)
        .clamp(config.planet_mass_earth_min, config.planet_mass_earth_max);
    let orbital_radius = Length::from_astronomical_units(orbit_au)?;
    let atmosphere = generate_planet_atmosphere(&mut rng, mass_earth, metallicity);
    let radius = Length::from_meters(planet_radius_meters(
        Mass::from_earth_masses(mass_earth)?,
        atmosphere,
    ))?;
    let equilibrium_temperature =
        Temperature::from_kelvin(equilibrium_temperature_kelvin(star_luminosity, orbit_au))?;
    let orbital_elements = OrbitalElements::new(
        orbital_radius,
        0.25 * rng.next_f64().powi(2),
        Angle::from_degrees(random_signed(&mut rng) * 5.0)?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        SimDuration::ZERO,
    )?;
    let (position, velocity) = propagate_kepler(
        orbital_elements,
        star_mass,
        SimDuration::ZERO,
        planet_position,
        planet_velocity,
    )?;
    Ok(Planet {
        id,
        system_id,
        star_id,
        mass: Mass::from_earth_masses(mass_earth)?,
        radius,
        atmosphere,
        equilibrium_temperature,
        orbital_radius,
        position,
        velocity,
        orbital_elements,
        status: PlanetStatus::Stable,
    })
}

fn generate_moon_count(rng: &mut StableRng, planet_mass: Mass) -> u8 {
    let mass_earth = planet_mass.earth_masses();
    if mass_earth < 0.2 {
        return 0;
    }
    let occurrence = (0.10 + 0.25 * mass_earth.ln().max(0.0)).clamp(0.10, 0.85);
    if rng.next_f64() >= occurrence {
        return 0;
    }
    let maximum: u64 = if mass_earth >= 20.0 {
        4
    } else if mass_earth >= 2.0 {
        2
    } else {
        1
    };
    1 + u8::try_from(rng.next_u64_below(maximum).unwrap_or(0)).unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn generate_moon(
    id: MoonId,
    planet_id: PlanetId,
    system_id: StarSystemId,
    planet_position: Position,
    planet_velocity: Velocity,
    planet_mass: Mass,
    star_mass: Mass,
    planet_orbit_au: f64,
    rng: &mut StableRng,
) -> Result<Option<Moon>, ObjectError> {
    let planet_radius = planet_radius_meters(planet_mass, PlanetAtmosphere::None);
    let hill_radius = planet_orbit_au
        * METERS_PER_ASTRONOMICAL_UNIT
        * (planet_mass.kilograms() / (3.0 * star_mass.kilograms())).cbrt();
    let minimum_axis = planet_radius * 6.0;
    let maximum_axis = hill_radius * 0.25;
    if !maximum_axis.is_finite() || maximum_axis <= minimum_axis {
        return Ok(None);
    }
    let axis = random_range(rng, minimum_axis, maximum_axis);
    let mass_earth = (random_range(rng, 0.00001, 0.02) * planet_mass.earth_masses().sqrt())
        .clamp(0.000001, 0.08);
    let orbital_elements = OrbitalElements::new(
        Length::from_meters(axis)?,
        random_range(rng, 0.0, 0.15),
        Angle::from_degrees(random_signed(rng) * 8.0)?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        Angle::from_radians(std::f64::consts::TAU * rng.next_f64())?,
        SimDuration::ZERO,
    )?;
    let moon_mass = Mass::from_earth_masses(mass_earth)?;
    let (position, velocity) = propagate_kepler(
        orbital_elements,
        moon_mass,
        SimDuration::ZERO,
        planet_position,
        planet_velocity,
    )?;
    Ok(Some(Moon {
        id,
        planet_id,
        system_id,
        mass: moon_mass,
        radius: Length::from_meters(planet_radius_meters(moon_mass, PlanetAtmosphere::None))?,
        orbital_radius: Length::from_meters(axis)?,
        position,
        velocity,
        orbital_elements,
    }))
}

fn generate_asteroid_belt(
    id: AsteroidBeltId,
    system_id: StarSystemId,
    rng: &mut StableRng,
    max_planet_orbit_au: f64,
) -> Result<AsteroidBelt, ObjectError> {
    let inner_au = if max_planet_orbit_au > 0.0 {
        max_planet_orbit_au * random_range(rng, 1.25, 2.0)
    } else {
        random_range(rng, 1.0, 3.0)
    };
    let outer_au = inner_au * random_range(rng, 1.4, 3.0);
    Ok(AsteroidBelt {
        id,
        system_id,
        inner_radius: Length::from_astronomical_units(inner_au)?,
        outer_radius: Length::from_astronomical_units(outer_au)?,
        estimated_body_count: 1_000_000 + rng.next_u64_below(9_000_000_000).unwrap_or(0),
        total_mass: Mass::from_earth_masses(random_range(rng, 0.00001, 0.25))?,
    })
}

fn generate_comet_cloud(
    id: CometCloudId,
    system_id: StarSystemId,
    rng: &mut StableRng,
) -> Result<CometCloud, ObjectError> {
    let inner_au = random_range(rng, 30.0, 80.0);
    let outer_au = inner_au * random_range(rng, 4.0, 20.0);
    Ok(CometCloud {
        id,
        system_id,
        inner_radius: Length::from_astronomical_units(inner_au)?,
        outer_radius: Length::from_astronomical_units(outer_au)?,
        estimated_body_count: 100_000 + rng.next_u64_below(10_000_000_000).unwrap_or(0),
        total_mass: Mass::from_earth_masses(random_range(rng, 0.001, 2.0))?,
    })
}

fn generate_debris_disk(
    id: DebrisDiskId,
    system_id: StarSystemId,
    rng: &mut StableRng,
) -> Result<DebrisDisk, ObjectError> {
    let inner_au = random_range(rng, 0.05, 0.5);
    let outer_au = inner_au * random_range(rng, 5.0, 30.0);
    Ok(DebrisDisk {
        id,
        system_id,
        inner_radius: Length::from_astronomical_units(inner_au)?,
        outer_radius: Length::from_astronomical_units(outer_au)?,
        estimated_body_count: 10_000 + rng.next_u64_below(100_000_000).unwrap_or(0),
        total_mass: Mass::from_earth_masses(random_range(rng, 0.000001, 0.05))?,
    })
}

fn generate_planet_atmosphere(
    rng: &mut StableRng,
    mass_earth: f64,
    metallicity: f64,
) -> PlanetAtmosphere {
    if mass_earth < 0.01 {
        PlanetAtmosphere::None
    } else if mass_earth >= 15.0 {
        PlanetAtmosphere::HydrogenHelium
    } else if mass_earth < 0.1 || rng.next_f64() > (0.45 + metallicity * 4.0).clamp(0.45, 0.90) {
        PlanetAtmosphere::Thin
    } else if rng.next_f64() < 0.70 {
        PlanetAtmosphere::NitrogenRich
    } else {
        PlanetAtmosphere::CarbonDioxideRich
    }
}

fn planet_radius_meters(mass: Mass, atmosphere: PlanetAtmosphere) -> f64 {
    let mass_earth = mass.earth_masses().max(0.000001);
    if atmosphere == PlanetAtmosphere::HydrogenHelium {
        11.2e6 * (mass_earth / 318.0).powf(0.01).clamp(0.75, 1.15)
    } else {
        6.371e6 * mass_earth.powf(0.27)
    }
}

fn equilibrium_temperature_kelvin(luminosity: Luminosity, orbit_au: f64) -> f64 {
    let luminosity_solar = luminosity.solar_luminosities().max(0.0);
    (278.5 * luminosity_solar.powf(0.25) / orbit_au.max(0.001).sqrt()).clamp(2.7, 10_000.0)
}

fn validate_small_body_population(
    inner_radius: Length,
    outer_radius: Length,
    estimated_body_count: u64,
    kind: &str,
) -> Result<(), ObjectError> {
    if inner_radius >= outer_radius || estimated_body_count == 0 {
        return Err(ObjectError(format!(
            "{kind} must have an increasing radial range and a nonzero population"
        )));
    }
    Ok(())
}

fn generate_planet_count(
    rng: &mut StableRng,
    config: GenerationConfig,
    star_mass: Mass,
    metallicity: f64,
) -> u8 {
    let occurrence = (0.25
        + 0.35 * (metallicity / 0.014).clamp(0.0, 2.0)
        + 0.15 * (1.0 / star_mass.solar_masses()).clamp(0.0, 2.0))
    .clamp(0.05, 0.95);
    if config.planets_per_system_min == 0 && rng.next_f64() > occurrence {
        return 0;
    }
    let span = u64::from(config.planets_per_system_max - config.planets_per_system_min + 1);
    config.planets_per_system_min
        + u8::try_from(
            rng.next_u64_below(span)
                .expect("planet count span is nonzero"),
        )
        .expect("planet count span fits in u8")
}

fn stable_orbital_radii(rng: &mut StableRng, count: u8, config: GenerationConfig) -> Vec<f64> {
    if count == 0 {
        return Vec::new();
    }
    if count == 1 {
        return vec![
            (config.orbit_au_min.ln()
                + rng.next_f64() * (config.orbit_au_max.ln() - config.orbit_au_min.ln()))
            .exp(),
        ];
    }
    let ratio = (config.orbit_au_max / config.orbit_au_min).powf(1.0 / f64::from(count - 1));
    (0..count)
        .map(|index| config.orbit_au_min * ratio.powi(i32::from(index)))
        .collect()
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

fn encode_ids<T>(ids: &[T]) -> String
where
    T: Copy,
    u64: From<T>,
{
    ids.iter()
        .map(|id| u64::from(*id).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn float_bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn encode_optional_orbit(orbit: Option<OrbitalElements>) -> String {
    let Some(orbit) = orbit else {
        return "-".to_owned();
    };
    [
        float_bits(orbit.semi_major_axis().meters()),
        float_bits(orbit.eccentricity()),
        float_bits(orbit.inclination().radians()),
        float_bits(orbit.longitude_of_ascending_node().radians()),
        float_bits(orbit.argument_of_periapsis().radians()),
        float_bits(orbit.mean_anomaly_at_epoch().radians()),
        orbit.epoch().seconds().to_string(),
    ]
    .join(",")
}

fn parse_optional_orbit(value: &str) -> Result<Option<OrbitalElements>, ObjectError> {
    if value == "-" {
        return Ok(None);
    }
    let fields: Vec<_> = value.split(',').collect();
    if fields.len() != 7 {
        return Err(ObjectError("malformed companion orbit".to_owned()));
    }
    Ok(Some(OrbitalElements::new(
        Length::from_meters(parse_float_bits(fields[0], "companion semi-major axis")?)?,
        parse_float_bits(fields[1], "companion eccentricity")?,
        Angle::from_radians(parse_float_bits(fields[2], "companion inclination")?)?,
        Angle::from_radians(parse_float_bits(fields[3], "companion ascending node")?)?,
        Angle::from_radians(parse_float_bits(fields[4], "companion periapsis")?)?,
        Angle::from_radians(parse_float_bits(fields[5], "companion mean anomaly")?)?,
        SimDuration::from_seconds(parse_number(fields[6], "companion orbit epoch")?),
    )?))
}

fn parse_float_bits(value: &str, field: &str) -> Result<f64, ObjectError> {
    u64::from_str_radix(value, 16)
        .map(f64::from_bits)
        .map_err(|_| ObjectError(format!("invalid {field}")))
}

fn parse_number<T>(value: &str, field: &str) -> Result<T, ObjectError>
where
    T: FromStr,
{
    value
        .parse()
        .map_err(|_| ObjectError(format!("invalid {field}")))
}

fn parse_position(fields: &[&str]) -> Result<Position, ObjectError> {
    Position::from_meters(
        parse_float_bits(fields[0], "position x")?,
        parse_float_bits(fields[1], "position y")?,
        parse_float_bits(fields[2], "position z")?,
    )
    .map_err(ObjectError::from)
}

fn parse_velocity(fields: &[&str]) -> Result<Velocity, ObjectError> {
    Velocity::from_meters_per_second(
        parse_float_bits(fields[0], "velocity x")?,
        parse_float_bits(fields[1], "velocity y")?,
        parse_float_bits(fields[2], "velocity z")?,
    )
    .map_err(ObjectError::from)
}

fn parse_ids<T>(value: &str, constructor: fn(u64) -> T) -> Result<Vec<T>, ObjectError> {
    if value.is_empty() {
        return Ok(Vec::new());
    }
    value
        .split(',')
        .map(|value| parse_number(value, "relationship ID").map(constructor))
        .collect()
}

fn require_field_count(fields: &[&str], expected: usize, kind: &str) -> Result<(), ObjectError> {
    if fields.len() == expected {
        Ok(())
    } else {
        Err(ObjectError(format!(
            "{kind} record has {} fields; expected {expected}",
            fields.len()
        )))
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
            stellar_fate(mass(100.0), 0.0001),
            StellarFate::PulsationalPairInstability
        );
        assert_eq!(
            stellar_fate(mass(180.0), 0.0001),
            StellarFate::PairInstability
        );
        assert_eq!(
            stellar_fate(mass(300.0), 0.0001),
            StellarFate::DirectCollapse
        );
        assert!(
            stellar_schedule(mass(100.0), 0.0001)
                .unwrap()
                .iter()
                .any(|(_, stage)| *stage == StellarStage::PulsationalPairInstability)
        );
        assert!(
            stellar_schedule(mass(180.0), 0.0001)
                .unwrap()
                .iter()
                .any(|(_, stage)| *stage == StellarStage::Disrupted)
        );
    }

    #[test]
    fn solar_mass_star_has_a_complete_low_mass_lifecycle() {
        let initial_mass = Mass::from_solar_masses(1.0).unwrap();
        let schedule = stellar_schedule(initial_mass, 0.014).unwrap();
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
            stellar_stage_at(initial_mass, 0.014, SimDuration::ZERO).unwrap(),
            StellarStage::Protostar
        );
        assert_eq!(
            stellar_stage_at(initial_mass, 0.014, schedule.last().unwrap().0).unwrap(),
            StellarStage::WhiteDwarf
        );
    }

    #[test]
    fn massive_solar_metallicity_star_explodes_after_a_wolf_rayet_phase() {
        let schedule = stellar_schedule(Mass::from_solar_masses(30.0).unwrap(), 0.014).unwrap();
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
