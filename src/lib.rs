//! Dependency-free foundations for a deterministic astrophysical simulator.

pub mod config;
pub mod events;
pub mod id;
pub mod objects;
pub mod orbital;
pub mod region;
pub mod rng;
pub mod simulation;
pub mod time;
pub mod units;

pub use config::{GenerationConfig, SimulationConfig};
pub use events::{AstroEvent, AstroEventKind};
pub use id::{GalaxyId, PlanetId, StarId, StarSystemId, UniverseId};
pub use objects::{
    DynamicsMode, Galaxy, GalaxyCatalogEntry, GalaxyMorphology, ObjectCounts, ObjectId, ObjectRef,
    ObjectRegistry, Planet, PlanetImpact, PlanetImpactKind, PlanetStatus, Star, StarSystem,
    StellarInteractionKind, StellarStage,
};
pub use orbital::{GravityBody, OrbitalElements, leapfrog_step, propagate_kepler};
pub use region::{GalaxyRegion, MaterializedRegion, RegionCoord};
pub use simulation::{MAX_CONCURRENT_GALAXIES, Universe};
pub use time::{SimDuration, SimulationClock, TickDuration};
