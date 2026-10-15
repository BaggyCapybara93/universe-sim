use std::fmt;

use crate::chemistry::GasComposition;
use crate::id::{
    AsteroidBeltId, AsteroidId, CometCloudId, CometId, DebrisDiskId, GalaxyId, GasCloudId, MoonId,
    PlanetId, ProtoplanetaryDiskId, StarClusterId, StarId, StarSystemId,
};
use crate::orbital::OrbitalElements;
use crate::rng::{StableRng, derive_object_seed};
use crate::time::{SECONDS_PER_YEAR, SimDuration};
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

impl GalaxyMorphology {
    pub(crate) const fn code(self) -> &'static str {
        match self {
            Self::Spiral => "spiral",
            Self::Elliptical => "elliptical",
            Self::Irregular => "irregular",
            Self::Dwarf => "dwarf",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "spiral" => Ok(Self::Spiral),
            "elliptical" => Ok(Self::Elliptical),
            "irregular" => Ok(Self::Irregular),
            "dwarf" => Ok(Self::Dwarf),
            _ => Err(ObjectError(format!("unknown galaxy morphology '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GalaxyCatalogEntry {
    pub(super) id: GalaxyId,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) morphology: GalaxyMorphology,
    pub(super) radius_parsecs: f64,
    pub(super) merged_into: Option<GalaxyId>,
}

impl GalaxyCatalogEntry {
    pub const fn id(self) -> GalaxyId {
        self.id
    }

    pub const fn position(self) -> Position {
        self.position
    }

    pub const fn velocity(self) -> Velocity {
        self.velocity
    }

    pub const fn morphology(self) -> GalaxyMorphology {
        self.morphology
    }

    pub const fn radius_parsecs(self) -> f64 {
        self.radius_parsecs
    }

    pub const fn merged_into(self) -> Option<GalaxyId> {
        self.merged_into
    }

    pub(crate) fn with_dynamics(
        mut self,
        position: Position,
        velocity: Velocity,
        morphology: GalaxyMorphology,
        radius_parsecs: f64,
        merged_into: Option<GalaxyId>,
    ) -> Self {
        self.position = position;
        self.velocity = velocity;
        self.morphology = morphology;
        self.radius_parsecs = radius_parsecs;
        self.merged_into = merged_into;
        self
    }

    pub(crate) fn primary(seed: u64, radius_parsecs: f64) -> Self {
        Self {
            id: GalaxyId::new(derive_object_seed(seed, 0, "galaxy-id-v1")),
            position: Position::default(),
            velocity: Velocity::default(),
            morphology: GalaxyMorphology::Spiral,
            radius_parsecs,
            merged_into: None,
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
        let tangent_seed = super::random_unit_vector(&mut rng);
        let projection: f64 = direction.iter().zip(tangent_seed).map(|(a, b)| a * b).sum();
        let mut tangent = [0.0; 3];
        for axis in 0..3 {
            tangent[axis] = tangent_seed[axis] - projection * direction[axis];
        }
        let norm = tangent
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            .sqrt();
        if norm > 0.0 {
            for component in &mut tangent {
                *component /= norm;
            }
        }
        let radial_speed = super::random_range(&mut rng, 150_000.0, 350_000.0);
        let transverse_speed = super::random_range(&mut rng, 0.0, 1_500.0);
        Ok(Self {
            id,
            position: Position::from_meters(
                distance * direction[0],
                distance * direction[1],
                distance * direction[2],
            )?,
            velocity: Velocity::from_meters_per_second(
                -radial_speed * direction[0] + transverse_speed * tangent[0],
                -radial_speed * direction[1] + transverse_speed * tangent[1],
                -radial_speed * direction[2] + transverse_speed * tangent[2],
            )?,
            morphology,
            radius_parsecs: base_radius_parsecs * radius_factor,
            merged_into: None,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GasCloudKind {
    Diffuse,
    Molecular,
}

impl GasCloudKind {
    pub(crate) const fn code(self) -> &'static str {
        match self {
            Self::Diffuse => "D",
            Self::Molecular => "M",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "D" => Ok(Self::Diffuse),
            "M" => Ok(Self::Molecular),
            _ => Err(ObjectError(format!("unknown gas-cloud kind '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct GasCloud {
    pub(super) id: GasCloudId,
    pub(super) galaxy_id: GalaxyId,
    pub(super) kind: GasCloudKind,
    pub(super) position: Position,
    pub(super) mass_solar: f64,
    pub(super) composition: GasComposition,
    pub(super) cluster_ids: Vec<StarClusterId>,
    pub(super) pending_shock: Option<CloudShock>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CloudShock {
    pub(super) source_star_id: StarId,
    pub(super) arrived_at: SimDuration,
    pub(super) strength: f64,
}

impl CloudShock {
    pub const fn source_star_id(self) -> StarId {
        self.source_star_id
    }
    pub const fn arrived_at(self) -> SimDuration {
        self.arrived_at
    }
    pub const fn strength(self) -> f64 {
        self.strength
    }
}

impl GasCloud {
    pub const fn id(&self) -> GasCloudId {
        self.id
    }
    pub const fn galaxy_id(&self) -> GalaxyId {
        self.galaxy_id
    }
    pub const fn kind(&self) -> GasCloudKind {
        self.kind
    }
    pub const fn position(&self) -> Position {
        self.position
    }
    pub const fn mass_solar(&self) -> f64 {
        self.mass_solar
    }
    pub const fn composition(&self) -> GasComposition {
        self.composition
    }
    pub fn cluster_ids(&self) -> &[StarClusterId] {
        &self.cluster_ids
    }
    pub const fn pending_shock(&self) -> Option<CloudShock> {
        self.pending_shock
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StarCluster {
    pub(super) id: StarClusterId,
    pub(super) cloud_id: GasCloudId,
    pub(super) formed_at: SimDuration,
    pub(super) position: Position,
    pub(super) system_ids: Vec<StarSystemId>,
    pub(super) triggering_star_id: Option<StarId>,
}

/// A coarse structural state for a cluster at a universal simulation time.
/// Membership in `StarCluster::system_ids` remains the permanent birth record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClusterPhase {
    Embedded,
    Bound,
    Dispersing,
    Dispersed,
}

impl fmt::Display for ClusterPhase {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Embedded => "embedded",
            Self::Bound => "bound",
            Self::Dispersing => "dispersing",
            Self::Dispersed => "dispersed",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterEvolution {
    phase: ClusterPhase,
    characteristic_radius_parsecs: f64,
    bound_system_count: usize,
}

impl ClusterEvolution {
    pub const fn phase(self) -> ClusterPhase {
        self.phase
    }
    pub const fn characteristic_radius_parsecs(self) -> f64 {
        self.characteristic_radius_parsecs
    }
    pub const fn bound_system_count(self) -> usize {
        self.bound_system_count
    }
}

impl StarCluster {
    pub const fn id(&self) -> StarClusterId {
        self.id
    }
    pub const fn cloud_id(&self) -> GasCloudId {
        self.cloud_id
    }
    pub const fn formed_at(&self) -> SimDuration {
        self.formed_at
    }
    pub const fn position(&self) -> Position {
        self.position
    }
    pub fn system_ids(&self) -> &[StarSystemId] {
        &self.system_ids
    }
    pub const fn triggering_star_id(&self) -> Option<StarId> {
        self.triggering_star_id
    }

    /// Deterministic association expansion and evaporation prescription.
    /// This is a population estimate, not an N-body orbit integration.
    pub fn evolution_at(&self, elapsed: SimDuration) -> ClusterEvolution {
        let age_myr = self.age_myr_at(elapsed);
        let bound_system_count = self
            .system_ids
            .iter()
            .enumerate()
            .filter(|(rank, _)| age_myr < self.escape_age_myr(*rank))
            .count();
        let phase = if age_myr < 3.0 {
            ClusterPhase::Embedded
        } else if bound_system_count == 0 {
            ClusterPhase::Dispersed
        } else if bound_system_count == self.system_ids.len() {
            ClusterPhase::Bound
        } else {
            ClusterPhase::Dispersing
        };
        ClusterEvolution {
            phase,
            characteristic_radius_parsecs: 0.1
                + 0.1 * age_myr.min(3.0)
                + 0.8 * (age_myr - 3.0).max(0.0),
            bound_system_count,
        }
    }

    fn age_myr_at(&self, elapsed: SimDuration) -> f64 {
        elapsed.seconds().saturating_sub(self.formed_at.seconds()) as f64
            / (1_000_000 * SECONDS_PER_YEAR) as f64
    }

    fn escape_age_myr(&self, rank: usize) -> f64 {
        3.0 + 50.0 * (self.system_ids.len() as f64 / (rank as f64 + 0.5)).ln()
    }

    /// Deterministic birth offset used by newly formed multi-system clusters.
    pub(super) fn member_birth_offset_meters(&self, system_id: StarSystemId) -> [f64; 3] {
        if self.formed_at == SimDuration::ZERO {
            return [0.0; 3];
        }
        let mut rng = StableRng::new(derive_object_seed(
            self.id.value(),
            system_id.value(),
            "cluster-member-motion-v1",
        ));
        let radius = 0.1 * rng.next_f64().sqrt() * METERS_PER_PARSEC;
        let angle = std::f64::consts::TAU * rng.next_f64();
        [radius * angle.cos(), radius * angle.sin(), 0.0]
    }

    /// Offset from the cluster birth center, including post-escape drift.
    pub(super) fn member_offset_at(
        &self,
        system_id: StarSystemId,
        elapsed: SimDuration,
    ) -> Option<[f64; 3]> {
        let rank = self.system_ids.iter().position(|id| *id == system_id)?;
        let age_myr = self.age_myr_at(elapsed);
        let escape_age = self.escape_age_myr(rank);
        let birth = self.member_birth_offset_meters(system_id);
        let mut rng = StableRng::new(derive_object_seed(
            self.id.value(),
            system_id.value(),
            "cluster-member-motion-v1",
        ));
        let _ = rng.next_f64();
        let angle = std::f64::consts::TAU * rng.next_f64();
        let travel_parsecs = 0.1 * age_myr.min(escape_age) + 0.8 * (age_myr - escape_age).max(0.0);
        let travel_meters = travel_parsecs * METERS_PER_PARSEC;
        Some([
            birth[0] + travel_meters * angle.cos(),
            birth[1] + travel_meters * angle.sin(),
            birth[2],
        ])
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StarSystem {
    pub(super) id: StarSystemId,
    pub(super) galaxy_id: GalaxyId,
    pub(super) cluster_id: Option<StarClusterId>,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) star_id: StarId,
    pub(super) companion_star_ids: Vec<StarId>,
    pub(super) planet_ids: Vec<PlanetId>,
    pub(super) moon_ids: Vec<MoonId>,
    pub(super) asteroid_ids: Vec<AsteroidId>,
    pub(super) comet_ids: Vec<CometId>,
    pub(super) asteroid_belt_ids: Vec<AsteroidBeltId>,
    pub(super) comet_cloud_ids: Vec<CometCloudId>,
    pub(super) debris_disk_ids: Vec<DebrisDiskId>,
    pub(super) protoplanetary_disk_ids: Vec<ProtoplanetaryDiskId>,
    pub(super) component: GalacticComponent,
}

impl StarSystem {
    pub const fn id(&self) -> StarSystemId {
        self.id
    }

    pub const fn galaxy_id(&self) -> GalaxyId {
        self.galaxy_id
    }

    pub const fn cluster_id(&self) -> Option<StarClusterId> {
        self.cluster_id
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

    pub fn asteroid_ids(&self) -> &[AsteroidId] {
        &self.asteroid_ids
    }

    pub fn comet_ids(&self) -> &[CometId] {
        &self.comet_ids
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

    pub fn protoplanetary_disk_ids(&self) -> &[ProtoplanetaryDiskId] {
        &self.protoplanetary_disk_ids
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
    pub(super) birth_composition: GasComposition,
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

    pub const fn birth_composition(&self) -> GasComposition {
        self.birth_composition
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
    TypeIaSupernova,
    Kilonova,
    BlackHoleMerger,
    CompactObjectMerger,
}

impl StellarInteractionKind {
    pub(crate) const fn code(self) -> &'static str {
        match self {
            Self::StarFormation => "star_formation",
            Self::MassTransfer => "mass_transfer",
            Self::StellarMerger => "stellar_merger",
            Self::PlanetaryCollision => "planetary_collision",
            Self::TypeIaSupernova => "type_ia_supernova",
            Self::Kilonova => "kilonova",
            Self::BlackHoleMerger => "black_hole_merger",
            Self::CompactObjectMerger => "compact_object_merger",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "star_formation" => Ok(Self::StarFormation),
            "mass_transfer" => Ok(Self::MassTransfer),
            "stellar_merger" => Ok(Self::StellarMerger),
            "planetary_collision" => Ok(Self::PlanetaryCollision),
            "type_ia_supernova" => Ok(Self::TypeIaSupernova),
            "kilonova" => Ok(Self::Kilonova),
            "black_hole_merger" => Ok(Self::BlackHoleMerger),
            "compact_object_merger" => Ok(Self::CompactObjectMerger),
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanetaryComposition {
    rocky: f64,
    icy: f64,
    iron: f64,
    hydrogen_helium: f64,
}

impl PlanetaryComposition {
    pub(crate) fn new(rocky: f64, icy: f64, iron: f64, hydrogen_helium: f64) -> Self {
        Self {
            rocky,
            icy,
            iron,
            hydrogen_helium,
        }
    }

    pub const fn rocky_fraction(self) -> f64 {
        self.rocky
    }
    pub const fn icy_fraction(self) -> f64 {
        self.icy
    }
    pub const fn iron_fraction(self) -> f64 {
        self.iron
    }
    pub const fn hydrogen_helium_fraction(self) -> f64 {
        self.hydrogen_helium
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoonOrigin {
    InSitu,
    Impact,
    Captured,
}

impl MoonOrigin {
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::InSitu => "in_situ",
            Self::Impact => "impact",
            Self::Captured => "captured",
        }
    }

    pub(crate) fn parse_code(value: &str) -> Result<Self, ObjectError> {
        match value {
            "in_situ" => Ok(Self::InSitu),
            "impact" => Ok(Self::Impact),
            "captured" => Ok(Self::Captured),
            _ => Err(ObjectError(format!("unknown moon origin '{value}'"))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MoonFormationHistory {
    origin: MoonOrigin,
    formed_at: SimDuration,
    captured_at: Option<SimDuration>,
}

impl MoonFormationHistory {
    pub(crate) const fn new(
        origin: MoonOrigin,
        formed_at: SimDuration,
        captured_at: Option<SimDuration>,
    ) -> Self {
        Self {
            origin,
            formed_at,
            captured_at,
        }
    }

    pub const fn origin(self) -> MoonOrigin {
        self.origin
    }
    pub const fn formed_at(self) -> SimDuration {
        self.formed_at
    }
    pub const fn captured_at(self) -> Option<SimDuration> {
        self.captured_at
    }
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
    pub(super) atmosphere_loss_fraction: f64,
    pub(super) equilibrium_temperature: Temperature,
    pub(super) in_habitable_zone: bool,
    pub(super) bulk_composition: PlanetaryComposition,
    pub(super) formed_at: SimDuration,
    pub(super) initial_orbital_radius: Length,
    pub(super) migration_rate_au_per_year: f64,
    pub(super) scattering_count: u32,
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

    pub const fn atmosphere_loss_fraction(&self) -> f64 {
        self.atmosphere_loss_fraction
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

    pub const fn is_in_habitable_zone(&self) -> bool {
        self.in_habitable_zone
    }

    pub const fn bulk_composition(&self) -> PlanetaryComposition {
        self.bulk_composition
    }

    pub const fn formed_at(&self) -> SimDuration {
        self.formed_at
    }
    pub const fn migration_rate_au_per_year(&self) -> f64 {
        self.migration_rate_au_per_year
    }
    pub const fn scattering_count(&self) -> u32 {
        self.scattering_count
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
    pub(super) atmosphere: PlanetAtmosphere,
    pub(super) equilibrium_temperature: Temperature,
    pub(super) formation_history: MoonFormationHistory,
    pub(super) tidal_migration_rate_m_per_year: f64,
    pub(super) initial_orbital_radius: Length,
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
    pub const fn atmosphere(&self) -> PlanetAtmosphere {
        self.atmosphere
    }
    pub const fn equilibrium_temperature(&self) -> Temperature {
        self.equilibrium_temperature
    }
    pub const fn formation_history(&self) -> MoonFormationHistory {
        self.formation_history
    }
    pub const fn tidal_migration_rate_m_per_year(&self) -> f64 {
        self.tidal_migration_rate_m_per_year
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
pub struct Asteroid {
    pub(super) id: AsteroidId,
    pub(super) belt_id: AsteroidBeltId,
    pub(super) system_id: StarSystemId,
    pub(super) mass: Mass,
    pub(super) radius: Length,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) orbital_elements: OrbitalElements,
}

impl Asteroid {
    pub const fn id(&self) -> AsteroidId {
        self.id
    }
    pub const fn belt_id(&self) -> AsteroidBeltId {
        self.belt_id
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
pub struct Comet {
    pub(super) id: CometId,
    pub(super) cloud_id: CometCloudId,
    pub(super) system_id: StarSystemId,
    pub(super) mass: Mass,
    pub(super) radius: Length,
    pub(super) volatile_fraction: f64,
    pub(super) position: Position,
    pub(super) velocity: Velocity,
    pub(super) orbital_elements: OrbitalElements,
}

impl Comet {
    pub const fn id(&self) -> CometId {
        self.id
    }
    pub const fn cloud_id(&self) -> CometCloudId {
        self.cloud_id
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
    pub const fn volatile_fraction(&self) -> f64 {
        self.volatile_fraction
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

#[derive(Debug, Clone, PartialEq)]
pub struct ProtoplanetaryDisk {
    pub(super) id: ProtoplanetaryDiskId,
    pub(super) system_id: StarSystemId,
    pub(super) inner_radius: Length,
    pub(super) outer_radius: Length,
    pub(super) initial_mass: Mass,
    pub(super) remaining_mass: Mass,
    pub(super) formed_at: SimDuration,
    pub(super) dissipates_at: SimDuration,
}

impl ProtoplanetaryDisk {
    pub const fn id(&self) -> ProtoplanetaryDiskId {
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
    pub const fn initial_mass(&self) -> Mass {
        self.initial_mass
    }
    pub const fn remaining_mass(&self) -> Mass {
        self.remaining_mass
    }
    pub const fn formed_at(&self) -> SimDuration {
        self.formed_at
    }
    pub const fn dissipates_at(&self) -> SimDuration {
        self.dissipates_at
    }
    pub fn is_active_at(&self, elapsed: SimDuration) -> bool {
        elapsed < self.dissipates_at
    }
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
    pub gas_clouds: usize,
    pub star_clusters: usize,
    pub star_systems: usize,
    pub stars: usize,
    pub planets: usize,
    pub moons: usize,
    pub asteroids: usize,
    pub comets: usize,
    pub asteroid_belts: usize,
    pub comet_clouds: usize,
    pub debris_disks: usize,
    pub protoplanetary_disks: usize,
    pub gas_reservoirs: usize,
}
