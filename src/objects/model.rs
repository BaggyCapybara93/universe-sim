use std::fmt;

use crate::id::{
    AsteroidBeltId, CometCloudId, DebrisDiskId, GalaxyId, MoonId, PlanetId, StarId, StarSystemId,
};
use crate::orbital::OrbitalElements;
use crate::rng::{StableRng, derive_object_seed};
use crate::time::SimDuration;
use crate::units::{Length, Luminosity, METERS_PER_PARSEC, Mass, Position, Temperature, Velocity};

use super::ObjectError;

#[derive(Debug, Clone, PartialEq)]
pub struct Galaxy {
    pub(super) id: GalaxyId,
    pub(super) position: Position,
    pub(super) star_system_ids: Vec<StarSystemId>,
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
    pub(super) id: GalaxyId,
    pub(super) position: Position,
    pub(super) morphology: GalaxyMorphology,
    pub(super) radius_parsecs: f64,
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
            GalaxyMorphology::Spiral => super::random_range(&mut rng, 0.7, 1.5),
            GalaxyMorphology::Elliptical => super::random_range(&mut rng, 0.5, 1.2),
            GalaxyMorphology::Irregular => super::random_range(&mut rng, 0.25, 0.7),
            GalaxyMorphology::Dwarf => super::random_range(&mut rng, 0.05, 0.25),
        };
        let distance = super::random_range(&mut rng, 500_000.0, 10_000_000.0) * METERS_PER_PARSEC;
        let direction = super::random_unit_vector(&mut rng);
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
    pub(super) id: StarSystemId,
    pub(super) galaxy_id: GalaxyId,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) star_id: StarId,
    pub(super) companion_star_ids: Vec<StarId>,
    pub(super) planet_ids: Vec<PlanetId>,
    pub(super) moon_ids: Vec<MoonId>,
    pub(super) asteroid_belt_ids: Vec<AsteroidBeltId>,
    pub(super) comet_cloud_ids: Vec<CometCloudId>,
    pub(super) debris_disk_ids: Vec<DebrisDiskId>,
    pub(super) component: GalacticComponent,
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
    pub(super) fn code(self) -> &'static str {
        match self {
            Self::Disk => "disk",
            Self::Bulge => "bulge",
            Self::Halo => "halo",
        }
    }

    pub(super) fn parse_code(value: &str) -> Result<Self, ObjectError> {
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
    pub(super) id: StarId,
    pub(super) system_id: StarSystemId,
    pub(super) initial_mass: Mass,
    pub(super) current_mass: Mass,
    pub(super) epoch_age: SimDuration,
    pub(super) formed_at: SimDuration,
    pub(super) age: SimDuration,
    pub(super) metallicity: f64,
    pub(super) luminosity: Luminosity,
    pub(super) temperature: Temperature,
    pub(super) stage: StellarStage,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) orbital_elements: Option<OrbitalElements>,
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
    pub(super) fn code(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Thin => "thin",
            Self::NitrogenRich => "nitrogen_rich",
            Self::CarbonDioxideRich => "carbon_dioxide_rich",
            Self::HydrogenHelium => "hydrogen_helium",
        }
    }

    pub(super) fn parse_code(value: &str) -> Result<Self, ObjectError> {
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
    pub(super) fn code(self) -> &'static str {
        match self {
            Self::Stable => "stable",
            Self::Irradiated => "irradiated",
            Self::Unbound => "unbound",
            Self::Destroyed => "destroyed",
        }
    }

    pub(super) fn parse_code(value: &str) -> Result<Self, ObjectError> {
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
    pub(super) id: PlanetId,
    pub(super) system_id: StarSystemId,
    pub(super) star_id: StarId,
    pub(super) mass: Mass,
    pub(super) radius: Length,
    pub(super) atmosphere: PlanetAtmosphere,
    pub(super) equilibrium_temperature: Temperature,
    pub(super) orbital_radius: Length,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) orbital_elements: OrbitalElements,
    pub(super) status: PlanetStatus,
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
    pub(super) id: MoonId,
    pub(super) planet_id: PlanetId,
    pub(super) system_id: StarSystemId,
    pub(super) mass: Mass,
    pub(super) radius: Length,
    pub(super) orbital_radius: Length,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) orbital_elements: OrbitalElements,
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
    pub(super) id: AsteroidBeltId,
    pub(super) system_id: StarSystemId,
    pub(super) inner_radius: Length,
    pub(super) outer_radius: Length,
    pub(super) estimated_body_count: u64,
    pub(super) total_mass: Mass,
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
    pub(super) id: CometCloudId,
    pub(super) system_id: StarSystemId,
    pub(super) inner_radius: Length,
    pub(super) outer_radius: Length,
    pub(super) estimated_body_count: u64,
    pub(super) total_mass: Mass,
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
    pub(super) id: DebrisDiskId,
    pub(super) system_id: StarSystemId,
    pub(super) inner_radius: Length,
    pub(super) outer_radius: Length,
    pub(super) estimated_body_count: u64,
    pub(super) total_mass: Mass,
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
