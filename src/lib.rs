//! Dependency-free foundations for a deterministic astrophysical simulator.

pub mod chemistry;
pub mod config;
pub mod events;
pub mod id;
pub mod matter;
pub mod objects;
pub mod orbital;
pub mod region;
pub mod rng;
pub mod simulation;
pub mod time;
pub mod units;

pub use chemistry::{CompositionError, GasComposition};
pub use config::{GenerationConfig, SimulationConfig};
pub use events::{AstroEvent, AstroEventKind};
pub use id::{
    AsteroidBeltId, AsteroidId, CometCloudId, CometId, DebrisDiskId, GalaxyId, GasCloudId, MoonId,
    PlanetId, ProtoplanetaryDiskId, StarClusterId, StarId, StarSystemId, UniverseId,
};
pub use matter::MatterInventory;
pub use objects::{
    Asteroid, AsteroidBelt, Comet, CometCloud, DebrisDisk, DynamicsMode, Galaxy,
    GalaxyCatalogEntry, GalaxyMorphology, GasCloud, GasCloudKind, Moon, MoonFormationHistory,
    MoonOrigin, ObjectCounts, ObjectId, ObjectRef, ObjectRegistry, Planet, PlanetAtmosphere,
    PlanetImpact, PlanetImpactKind, PlanetStatus, PlanetaryComposition, ProtoplanetaryDisk, Star,
    StarCluster, StarSystem, StellarInteractionKind, StellarStage,
};
pub use orbital::{GravityBody, OrbitalElements, leapfrog_step, propagate_kepler};
pub use region::{GalaxyRegion, MaterializedRegion, RegionCoord};
pub use simulation::{MAX_CONCURRENT_GALAXIES, Universe};
pub use time::{SimDuration, SimulationClock, TickDuration};
