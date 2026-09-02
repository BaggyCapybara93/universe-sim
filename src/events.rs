use crate::id::StarId;
use crate::objects::{PlanetImpact, StellarInteractionKind, StellarStage};
use crate::region::RegionCoord;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AstroEventKind {
    StellarStageChanged,
    StellarIgnition,
    PlanetaryNebulaFormation,
    CoreCollapseSupernova,
    NeutronStarFormation,
    PulsationalPairInstability,
    PairInstabilitySupernova,
    BlackHoleFormation,
    DirectCollapse,
    StarFormation,
    BinaryMassTransfer,
    StellarMerger,
    PlanetaryCollision,
}

impl std::fmt::Display for AstroEventKind {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::StellarStageChanged => formatter.write_str("stellar stage change"),
            Self::StellarIgnition => formatter.write_str("stellar ignition"),
            Self::PlanetaryNebulaFormation => formatter.write_str("planetary-nebula formation"),
            Self::CoreCollapseSupernova => formatter.write_str("core-collapse supernova"),
            Self::NeutronStarFormation => formatter.write_str("neutron-star formation"),
            Self::PulsationalPairInstability => formatter.write_str("pulsational pair-instability"),
            Self::PairInstabilitySupernova => formatter.write_str("pair-instability supernova"),
            Self::BlackHoleFormation => formatter.write_str("black-hole formation"),
            Self::DirectCollapse => formatter.write_str("direct collapse"),
            Self::StarFormation => formatter.write_str("star formation"),
            Self::BinaryMassTransfer => formatter.write_str("binary mass transfer"),
            Self::StellarMerger => formatter.write_str("stellar merger"),
            Self::PlanetaryCollision => formatter.write_str("planetary collision"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AstroEvent {
    id: u64,
    tick: u64,
    galaxy_id: crate::id::GalaxyId,
    region: RegionCoord,
    star_id: StarId,
    from: StellarStage,
    to: StellarStage,
    impacts: Vec<PlanetImpact>,
    natal_kick_meters_per_second: u32,
    interaction: Option<StellarInteractionKind>,
    related_star_id: Option<StarId>,
}

pub(crate) struct AstroEventDetails {
    pub impacts: Vec<PlanetImpact>,
    pub natal_kick_meters_per_second: u32,
    pub interaction: Option<StellarInteractionKind>,
    pub related_star_id: Option<StarId>,
}

pub(crate) struct AstroEventLocation {
    pub galaxy_id: crate::id::GalaxyId,
    pub region: RegionCoord,
    pub star_id: StarId,
}

impl AstroEvent {
    pub(crate) fn new(
        id: u64,
        tick: u64,
        location: AstroEventLocation,
        from: StellarStage,
        to: StellarStage,
        details: AstroEventDetails,
    ) -> Self {
        Self {
            id,
            tick,
            galaxy_id: location.galaxy_id,
            region: location.region,
            star_id: location.star_id,
            from,
            to,
            impacts: details.impacts,
            natal_kick_meters_per_second: details.natal_kick_meters_per_second,
            interaction: details.interaction,
            related_star_id: details.related_star_id,
        }
    }

    pub const fn id(&self) -> u64 {
        self.id
    }

    pub const fn tick(&self) -> u64 {
        self.tick
    }

    pub const fn region(&self) -> RegionCoord {
        self.region
    }

    pub const fn galaxy_id(&self) -> crate::id::GalaxyId {
        self.galaxy_id
    }

    pub const fn star_id(&self) -> StarId {
        self.star_id
    }

    pub const fn from(&self) -> StellarStage {
        self.from
    }

    pub const fn to(&self) -> StellarStage {
        self.to
    }

    pub fn impacts(&self) -> &[PlanetImpact] {
        &self.impacts
    }

    pub const fn natal_kick_meters_per_second(&self) -> u32 {
        self.natal_kick_meters_per_second
    }

    pub const fn related_star_id(&self) -> Option<StarId> {
        self.related_star_id
    }

    pub(crate) const fn interaction(&self) -> Option<StellarInteractionKind> {
        self.interaction
    }

    pub const fn kind(&self) -> AstroEventKind {
        if let Some(interaction) = self.interaction {
            return match interaction {
                StellarInteractionKind::StarFormation => AstroEventKind::StarFormation,
                StellarInteractionKind::MassTransfer => AstroEventKind::BinaryMassTransfer,
                StellarInteractionKind::StellarMerger => AstroEventKind::StellarMerger,
                StellarInteractionKind::PlanetaryCollision => AstroEventKind::PlanetaryCollision,
            };
        }
        match self.to {
            StellarStage::MainSequence => AstroEventKind::StellarIgnition,
            StellarStage::PlanetaryNebula => AstroEventKind::PlanetaryNebulaFormation,
            StellarStage::CoreCollapseSupernova => AstroEventKind::CoreCollapseSupernova,
            StellarStage::NeutronStar => AstroEventKind::NeutronStarFormation,
            StellarStage::PulsationalPairInstability => AstroEventKind::PulsationalPairInstability,
            StellarStage::Disrupted => AstroEventKind::PairInstabilitySupernova,
            StellarStage::BlackHole
                if matches!(
                    self.from,
                    StellarStage::PulsationalPairInstability | StellarStage::CoreCollapseSupernova
                ) =>
            {
                AstroEventKind::BlackHoleFormation
            }
            StellarStage::BlackHole => AstroEventKind::DirectCollapse,
            _ => AstroEventKind::StellarStageChanged,
        }
    }
}
