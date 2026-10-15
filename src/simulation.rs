use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::ffi::OsString;
use std::fmt;
use std::fmt::Write as _;
use std::fs::{self, File, OpenOptions};
use std::io::Write as _;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use crate::config::{ConfigValidationError, GenerationConfig, SimulationConfig};
use crate::events::{
    AstroEvent, AstroEventDetails, AstroEventLocation, GalaxyEncounter, GalaxyEncounterKind,
};
use crate::galaxy::{GalaxyMotion, advance_encounters, encounter_id};
use crate::id::{GalaxyId, UniverseId};
use crate::matter::MatterInventory;
use crate::objects::{
    CentralBlackHole, GalaxyCatalogEntry, ObjectCounts, ObjectError, ObjectRegistry, PlanetImpact,
    PlanetImpactKind, StellarTransition,
};
use crate::region::{GalaxyRegion, MaterializedRegion, RegionCoord};
use crate::rng::{derive_object_seed, derive_region_seed};
use crate::time::{SimDuration, SimulationClock, TickDuration, TimeError};
use crate::units::{Position, Velocity};

pub const STATE_FORMAT_VERSION: u32 = 19;
pub const SIMULATION_VERSION: u32 = 8;
pub const GENERATION_VERSION: u32 = 15;
pub const RNG_VERSION: u32 = 1;
pub const MAX_CONCURRENT_GALAXIES: usize = 4;

const STATE_HEADER: &str = "COSMOS_STATE";
static TEMP_FILE_SEQUENCE: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone, PartialEq)]
pub struct Universe {
    id: UniverseId,
    config: SimulationConfig,
    clock: SimulationClock,
    galaxies: Vec<GalaxyCatalogEntry>,
    galaxy_motion: BTreeMap<GalaxyId, GalaxyMotion>,
    galaxy_encounters: Vec<GalaxyEncounter>,
    objects: ObjectRegistry,
    known_regions: BTreeMap<GalaxyRegion, u64>,
    regions: BTreeMap<GalaxyRegion, MaterializedRegion>,
    events: Vec<AstroEvent>,
}

impl Universe {
    pub fn new(config: SimulationConfig) -> Result<Self, UniverseCreateError> {
        config
            .validate()
            .map_err(UniverseCreateError::Configuration)?;
        let id = UniverseId::new(derive_object_seed(config.seed, 0, "universe"));
        let galaxies: Vec<_> = (0..config.generation.galaxies)
            .map(|index| {
                GalaxyCatalogEntry::generate(
                    config.seed,
                    index,
                    config.generation.galaxy_radius_parsecs,
                )
            })
            .collect::<Result<_, _>>()
            .map_err(UniverseCreateError::Objects)?;
        let primary_galaxy = galaxies[0];
        let galaxy_motion = galaxies
            .iter()
            .map(|galaxy| {
                (
                    galaxy.id(),
                    GalaxyMotion::from_catalog(*galaxy, SimDuration::ZERO),
                )
            })
            .collect();
        let objects =
            ObjectRegistry::generate_in_galaxy(config.seed, config.generation, primary_galaxy)
                .map_err(UniverseCreateError::Objects)?;
        let origin = GalaxyRegion::new(primary_galaxy.id(), RegionCoord::ORIGIN);
        Ok(Self {
            id,
            config,
            clock: SimulationClock::new(config.tick_duration),
            galaxies,
            galaxy_motion,
            galaxy_encounters: Vec::new(),
            objects,
            known_regions: BTreeMap::from([(origin, 0)]),
            regions: BTreeMap::new(),
            events: Vec::new(),
        })
    }

    pub const fn id(&self) -> UniverseId {
        self.id
    }

    pub const fn seed(&self) -> u64 {
        self.config.seed
    }

    pub const fn config(&self) -> SimulationConfig {
        self.config
    }

    pub const fn clock(&self) -> SimulationClock {
        self.clock
    }

    pub const fn objects(&self) -> &ObjectRegistry {
        &self.objects
    }

    /// Counts explicit objects in every currently materialized galaxy region.
    pub fn materialized_counts(&self) -> ObjectCounts {
        let mut total = self.objects.counts();
        for region in self.regions.values() {
            let counts = region.objects().counts();
            total.gas_clouds += counts.gas_clouds;
            total.star_clusters += counts.star_clusters;
            total.star_systems += counts.star_systems;
            total.stars += counts.stars;
            total.planets += counts.planets;
            total.moons += counts.moons;
            total.asteroids += counts.asteroids;
            total.comets += counts.comets;
            total.asteroid_belts += counts.asteroid_belts;
            total.comet_clouds += counts.comet_clouds;
            total.debris_disks += counts.debris_disks;
            total.protoplanetary_disks += counts.protoplanetary_disks;
            total.gas_reservoirs += counts.gas_reservoirs;
        }
        total.galaxies = self.galaxies.len();
        total
    }

    /// Inventories for materialized regions only; unloaded regions have no
    /// current detailed object graph to total.
    pub fn matter_inventories(&self) -> impl Iterator<Item = (GalaxyRegion, MatterInventory)> + '_ {
        let origin = GalaxyRegion::new(self.primary_galaxy_id(), RegionCoord::ORIGIN);
        std::iter::once((origin, self.objects.matter_inventory())).chain(
            self.regions
                .iter()
                .map(|(key, region)| (*key, region.objects().matter_inventory())),
        )
    }

    pub fn galaxies(&self) -> impl Iterator<Item = GalaxyCatalogEntry> + '_ {
        self.galaxies.iter().map(|galaxy| {
            self.galaxy_motion[&galaxy.id()]
                .snapshot(*galaxy, self.clock.elapsed())
                .expect("galaxy position was validated during advancement")
        })
    }

    pub fn galaxy_encounters(&self) -> &[GalaxyEncounter] {
        &self.galaxy_encounters
    }

    /// Current central black-hole and AGN state for a galaxy or its merger remnant.
    pub fn central_black_hole(&self, galaxy_id: GalaxyId) -> Option<CentralBlackHole> {
        let remnant = self.surviving_galaxy_id(galaxy_id)?;
        let progenitors: BTreeSet<_> = self
            .galaxies
            .iter()
            .filter_map(|galaxy| {
                (self.surviving_galaxy_id(galaxy.id()) == Some(remnant)).then_some(galaxy.id())
            })
            .collect();
        Some(crate::galaxy::estimate_central_black_hole(
            &self.galaxies,
            &progenitors,
            &self.galaxy_encounters,
            self.clock.elapsed(),
        ))
    }

    pub fn active_galaxy_count(&self) -> usize {
        self.galaxy_motion
            .values()
            .filter(|motion| motion.merged_into.is_none())
            .count()
    }

    pub fn surviving_galaxy_id(&self, galaxy_id: GalaxyId) -> Option<GalaxyId> {
        let mut current = galaxy_id;
        for _ in 0..self.galaxies.len() {
            match self.galaxy_motion.get(&current)?.merged_into {
                Some(next) => current = next,
                None => return Some(current),
            }
        }
        None
    }

    pub fn primary_galaxy_id(&self) -> GalaxyId {
        self.galaxies[0].id()
    }

    pub fn events(&self) -> &[AstroEvent] {
        &self.events
    }

    pub fn region_objects(&self, coordinate: RegionCoord) -> Option<&ObjectRegistry> {
        self.region_objects_in(self.primary_galaxy_id(), coordinate)
    }

    pub fn region_objects_in(
        &self,
        galaxy_id: GalaxyId,
        coordinate: RegionCoord,
    ) -> Option<&ObjectRegistry> {
        let key = GalaxyRegion::new(galaxy_id, coordinate);
        if galaxy_id == self.primary_galaxy_id() && coordinate == RegionCoord::ORIGIN {
            Some(&self.objects)
        } else {
            self.regions.get(&key).map(MaterializedRegion::objects)
        }
    }

    pub fn known_regions(&self) -> impl Iterator<Item = (RegionCoord, u64, bool)> + '_ {
        let primary = self.primary_galaxy_id();
        self.galaxy_regions()
            .filter(move |(key, _, _)| key.galaxy_id() == primary)
            .map(|(key, tick, materialized)| (key.coordinate(), tick, materialized))
    }

    pub fn galaxy_regions(&self) -> impl Iterator<Item = (GalaxyRegion, u64, bool)> + '_ {
        let primary_origin = GalaxyRegion::new(self.primary_galaxy_id(), RegionCoord::ORIGIN);
        self.known_regions.iter().map(move |(key, tick)| {
            (
                *key,
                *tick,
                *key == primary_origin || self.regions.contains_key(key),
            )
        })
    }

    pub fn materialize_region(&mut self, coordinate: RegionCoord) -> Result<bool, RegionError> {
        self.materialize_region_in(self.primary_galaxy_id(), coordinate)
    }

    /// Ensures every cataloged galaxy has an explicit origin population.
    /// Returns the number of regions newly materialized. Existing regions are
    /// left intact, so this can also upgrade a previously sparse saved state.
    pub fn materialize_all_galaxy_origins(&mut self) -> Result<usize, RegionError> {
        let galaxy_ids: Vec<_> = self.galaxies.iter().map(|galaxy| galaxy.id()).collect();
        let mut created = 0;
        for galaxy_id in galaxy_ids {
            if self.materialize_region_in(galaxy_id, RegionCoord::ORIGIN)? {
                created += 1;
            }
        }
        Ok(created)
    }

    pub fn materialize_region_in(
        &mut self,
        galaxy_id: GalaxyId,
        coordinate: RegionCoord,
    ) -> Result<bool, RegionError> {
        let key = GalaxyRegion::new(galaxy_id, coordinate);
        if let Some(survivor) = self.surviving_galaxy_id(galaxy_id)
            && survivor != galaxy_id
            && !self.known_regions.contains_key(&key)
        {
            return Err(RegionError::MergedGalaxy {
                galaxy_id,
                survivor,
            });
        }
        let galaxy = self
            .galaxies
            .iter()
            .find(|galaxy| galaxy.id() == galaxy_id)
            .copied()
            .ok_or(RegionError::UnknownGalaxy(galaxy_id))?;
        if (galaxy_id == self.primary_galaxy_id() && coordinate == RegionCoord::ORIGIN)
            || self.regions.contains_key(&key)
        {
            return Ok(false);
        }
        let region_seed = scoped_region_seed(self.seed(), galaxy_id, coordinate);
        let mut objects =
            ObjectRegistry::generate_in_galaxy(region_seed, self.config.generation, galaxy)
                .map_err(RegionError::Generation)?;
        let mut transitions = objects
            .form_stars_between_with_encounters(
                region_seed,
                self.config.generation,
                SimDuration::ZERO,
                self.clock.elapsed(),
                self.clock.tick_duration().duration(),
                &self.galaxy_encounters,
            )
            .map_err(RegionError::Generation)?;
        transitions.extend(
            objects
                .advance_stars_to(self.clock.elapsed(), self.clock.tick_duration().duration())
                .map_err(RegionError::Time)?,
        );
        objects
            .advance_planetary_formation_to(
                region_seed,
                self.config.generation,
                self.clock.elapsed(),
            )
            .map_err(RegionError::Generation)?;
        objects
            .advance_orbits_to(self.clock.elapsed())
            .map_err(RegionError::Generation)?;
        self.regions.insert(
            key,
            MaterializedRegion::new(coordinate, self.clock.tick(), objects),
        );
        self.known_regions.insert(key, self.clock.tick());
        self.record_transitions(key, transitions);
        Ok(true)
    }

    pub fn unload_region(&mut self, coordinate: RegionCoord) -> Result<bool, RegionError> {
        self.unload_region_in(self.primary_galaxy_id(), coordinate)
    }

    pub fn unload_region_in(
        &mut self,
        galaxy_id: GalaxyId,
        coordinate: RegionCoord,
    ) -> Result<bool, RegionError> {
        let key = GalaxyRegion::new(galaxy_id, coordinate);
        if galaxy_id == self.primary_galaxy_id() && coordinate == RegionCoord::ORIGIN {
            return Err(RegionError::CannotUnloadOrigin);
        }
        if !self.galaxies.iter().any(|galaxy| galaxy.id() == galaxy_id) {
            return Err(RegionError::UnknownGalaxy(galaxy_id));
        }
        Ok(self.regions.remove(&key).is_some())
    }

    pub fn advance_tick(&mut self) -> Result<(), SimulationAdvanceError> {
        self.advance_ticks(1)
    }

    /// Advances the simulation and returns only the events triggered in this interval.
    pub fn advance_ticks_with_events(
        &mut self,
        ticks: u64,
    ) -> Result<&[AstroEvent], SimulationAdvanceError> {
        self.advance_ticks_with_events_bounded(
            ticks,
            NonZeroUsize::new(MAX_CONCURRENT_GALAXIES).expect("limit is non-zero"),
        )
    }

    pub fn advance_ticks_with_events_bounded(
        &mut self,
        ticks: u64,
        max_concurrent_galaxies: NonZeroUsize,
    ) -> Result<&[AstroEvent], SimulationAdvanceError> {
        let first_new_event = self.events.len();
        self.advance_ticks_bounded(ticks, max_concurrent_galaxies)?;
        Ok(&self.events[first_new_event..])
    }

    pub fn advance_ticks(&mut self, ticks: u64) -> Result<(), SimulationAdvanceError> {
        self.advance_ticks_bounded(
            ticks,
            NonZeroUsize::new(MAX_CONCURRENT_GALAXIES).expect("limit is non-zero"),
        )
    }

    pub fn advance_ticks_bounded(
        &mut self,
        ticks: u64,
        max_concurrent_galaxies: NonZeroUsize,
    ) -> Result<(), SimulationAdvanceError> {
        let mut next_clock = self.clock;
        next_clock
            .advance_by(ticks)
            .map_err(SimulationAdvanceError::Time)?;
        let universe_seed = self.seed();
        let primary_galaxy_id = self.primary_galaxy_id();
        advance_encounters(
            &mut self.galaxy_motion,
            &mut self.galaxy_encounters,
            universe_seed,
            primary_galaxy_id,
            self.clock.elapsed(),
            next_clock.elapsed(),
            next_clock.tick_duration().duration(),
        )
        .map_err(SimulationAdvanceError::Time)?;
        let mut origin_transitions = self
            .objects
            .form_stars_between_with_encounters(
                self.seed(),
                self.config.generation,
                self.clock.elapsed(),
                next_clock.elapsed(),
                next_clock.tick_duration().duration(),
                &self.galaxy_encounters,
            )
            .map_err(SimulationAdvanceError::Objects)?;
        origin_transitions.extend(
            self.objects
                .advance_stars_to(next_clock.elapsed(), next_clock.tick_duration().duration())
                .map_err(SimulationAdvanceError::Time)?,
        );
        self.objects
            .advance_planetary_formation_to(
                self.seed(),
                self.config.generation,
                next_clock.elapsed(),
            )
            .map_err(SimulationAdvanceError::Objects)?;
        self.objects
            .advance_orbits_to(next_clock.elapsed())
            .map_err(SimulationAdvanceError::Objects)?;
        let generation = self.config.generation;
        let universe_seed = self.config.seed;
        let current_elapsed = self.clock.elapsed();
        let galaxy_encounters = &self.galaxy_encounters;
        let update_region = |(key, region): (&GalaxyRegion, &mut MaterializedRegion)| {
            let region_seed = scoped_region_seed(universe_seed, key.galaxy_id(), key.coordinate());
            let mut transitions = region
                .objects_mut()
                .form_stars_between_with_encounters(
                    region_seed,
                    generation,
                    current_elapsed,
                    next_clock.elapsed(),
                    next_clock.tick_duration().duration(),
                    galaxy_encounters,
                )
                .map_err(SimulationAdvanceError::Objects)?;
            transitions.extend(
                region
                    .objects_mut()
                    .advance_stars_to(next_clock.elapsed(), next_clock.tick_duration().duration())
                    .map_err(SimulationAdvanceError::Time)?,
            );
            region
                .objects_mut()
                .advance_planetary_formation_to(region_seed, generation, next_clock.elapsed())
                .map_err(SimulationAdvanceError::Objects)?;
            region
                .objects_mut()
                .advance_orbits_to(next_clock.elapsed())
                .map_err(SimulationAdvanceError::Objects)?;
            region.set_last_updated_tick(next_clock.tick());
            Ok((*key, transitions))
        };
        let active_galaxies: BTreeSet<_> = self
            .regions
            .keys()
            .map(|region| region.galaxy_id())
            .collect();
        let active_galaxies: Vec<_> = active_galaxies.into_iter().collect();
        let galaxy_batch_size = max_concurrent_galaxies.get().min(MAX_CONCURRENT_GALAXIES);
        let mut region_transitions = Vec::with_capacity(self.regions.len());
        for galaxy_batch in active_galaxies.chunks(galaxy_batch_size) {
            let batch: Result<Vec<_>, SimulationAdvanceError> = if self.regions.len() >= 2 {
                self.regions
                    .par_iter_mut()
                    .filter(|(key, _)| galaxy_batch.contains(&key.galaxy_id()))
                    .map(&update_region)
                    .collect()
            } else {
                self.regions
                    .iter_mut()
                    .filter(|(key, _)| galaxy_batch.contains(&key.galaxy_id()))
                    .map(&update_region)
                    .collect()
            };
            region_transitions.extend(batch?);
        }
        self.clock = next_clock;
        let primary_origin = GalaxyRegion::new(self.primary_galaxy_id(), RegionCoord::ORIGIN);
        self.known_regions.insert(primary_origin, self.clock.tick());
        for key in self.regions.keys() {
            self.known_regions.insert(*key, self.clock.tick());
        }
        self.record_transitions(primary_origin, origin_transitions);
        for (key, transitions) in region_transitions {
            self.record_transitions(key, transitions);
        }
        Ok(())
    }

    fn record_transitions(&mut self, key: GalaxyRegion, transitions: Vec<StellarTransition>) {
        let coordinate = key.coordinate();
        let region_seed = scoped_region_seed(self.seed(), key.galaxy_id(), coordinate);
        for transition in transitions {
            let namespace = event_namespace(
                transition.tick,
                transition.to,
                transition.interaction,
                transition.related_star_id,
            );
            let id = derive_object_seed(region_seed, transition.star_id.value(), &namespace);
            if self.events.iter().any(|event| event.id() == id) {
                continue;
            }
            self.events.push(AstroEvent::new(
                id,
                transition.tick,
                AstroEventLocation {
                    galaxy_id: key.galaxy_id(),
                    region: coordinate,
                    star_id: transition.star_id,
                },
                transition.from,
                transition.to,
                AstroEventDetails {
                    impacts: transition.impacts,
                    natal_kick_meters_per_second: transition.natal_kick_meters_per_second,
                    interaction: transition.interaction,
                    related_star_id: transition.related_star_id,
                },
            ));
        }
        self.events.sort_by_key(|event| (event.tick(), event.id()));
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), StateError> {
        atomic_write(path.as_ref(), self.encode().as_bytes())
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self, StateError> {
        let contents = fs::read_to_string(path).map_err(StateError::Io)?;
        Self::decode(&contents)
    }

    fn encode(&self) -> String {
        let generation = self.config.generation;
        let mut output = format!(
            "{STATE_HEADER}\nstate_format_version={STATE_FORMAT_VERSION}\nsimulation_version={SIMULATION_VERSION}\ngeneration_version={GENERATION_VERSION}\nrng_version={RNG_VERSION}\ncrate_version={}\nseed={}\nuniverse_id={}\ntick={}\nelapsed_seconds={}\ntick_duration_seconds={}\n",
            env!("CARGO_PKG_VERSION"),
            self.seed(),
            self.id.value(),
            self.clock.tick(),
            self.clock.elapsed().seconds(),
            self.clock.tick_duration().duration().seconds(),
        );
        writeln!(output, "star_systems={}", generation.star_systems)
            .expect("writing to a String cannot fail");
        writeln!(output, "galaxies={}", generation.galaxies)
            .expect("writing to a String cannot fail");
        writeln!(
            output,
            "planets_per_system_min={}",
            generation.planets_per_system_min
        )
        .expect("writing to a String cannot fail");
        writeln!(
            output,
            "planets_per_system_max={}",
            generation.planets_per_system_max
        )
        .expect("writing to a String cannot fail");
        for (key, value) in [
            ("star_mass_solar_min", generation.star_mass_solar_min),
            ("star_mass_solar_max", generation.star_mass_solar_max),
            ("planet_mass_earth_min", generation.planet_mass_earth_min),
            ("planet_mass_earth_max", generation.planet_mass_earth_max),
            ("orbit_au_min", generation.orbit_au_min),
            ("orbit_au_max", generation.orbit_au_max),
            ("galaxy_radius_parsecs", generation.galaxy_radius_parsecs),
        ] {
            writeln!(output, "{key}={:016x}", value.to_bits())
                .expect("writing to a String cannot fail");
        }
        self.objects.encode(&mut output);
        output.push_str("[regions]\n");
        for (key, last_updated_tick, materialized) in self.galaxy_regions() {
            let coordinate = key.coordinate();
            writeln!(
                output,
                "R|{}|{}|{}|{}|{}|{}",
                key.galaxy_id().value(),
                coordinate.x,
                coordinate.y,
                coordinate.z,
                last_updated_tick,
                u8::from(materialized)
            )
            .expect("writing to a String cannot fail");
        }
        output.push_str("[events]\n");
        for event in &self.events {
            writeln!(
                output,
                "E|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                event.id(),
                event.tick(),
                event.galaxy_id().value(),
                event.region().x,
                event.region().y,
                event.region().z,
                event.star_id().value(),
                event.from().code(),
                event.to().code(),
                encode_impacts(event.impacts()),
                event.natal_kick_meters_per_second(),
                event.interaction().map_or("-", |kind| kind.code()),
                event
                    .related_star_id()
                    .map_or_else(|| "-".to_owned(), |id| id.value().to_string()),
            )
            .expect("writing to a String cannot fail");
        }
        output.push_str("[galaxy_motion]\n");
        for (id, motion) in &self.galaxy_motion {
            let [px, py, pz] = motion.position_at_epoch.meters();
            let [vx, vy, vz] = motion.velocity.meters_per_second();
            writeln!(
                output,
                "D|{}|{:016x}|{:016x}|{:016x}|{:016x}|{:016x}|{:016x}|{}|{}|{:016x}|{}",
                id.value(),
                px.to_bits(),
                py.to_bits(),
                pz.to_bits(),
                vx.to_bits(),
                vy.to_bits(),
                vz.to_bits(),
                motion.epoch.seconds(),
                motion.morphology.code(),
                motion.radius_parsecs.to_bits(),
                motion
                    .merged_into
                    .map_or_else(|| "-".to_owned(), |id| id.value().to_string()),
            )
            .expect("writing to a String cannot fail");
        }
        output.push_str("[galaxy_encounters]\n");
        for encounter in &self.galaxy_encounters {
            writeln!(
                output,
                "C|{}|{}|{}|{}|{}|{}",
                encounter.id(),
                encounter.occurred_at().seconds(),
                encounter.tick(),
                encounter.first().value(),
                encounter.second().value(),
                encounter.kind().code(),
            )
            .expect("writing to a String cannot fail");
        }
        output
    }

    fn decode(contents: &str) -> Result<Self, StateError> {
        let mut lines = contents.lines();
        let header = lines.next();
        if header == Some("COSMOS_STATE_V1") {
            return Self::decode_v1(lines);
        }
        if header != Some(STATE_HEADER) {
            return Err(StateError::Compatibility(
                "unsupported state header; expected COSMOS_STATE".to_owned(),
            ));
        }

        let mut fields = BTreeMap::new();
        let mut found_objects = false;
        for line in lines.by_ref() {
            if line == "[objects]" {
                found_objects = true;
                break;
            }
            let (key, value) = line
                .split_once('=')
                .ok_or_else(|| StateError::Invalid(format!("malformed state line '{line}'")))?;
            if fields.insert(key, value).is_some() {
                return Err(StateError::Invalid(format!(
                    "duplicate state field '{key}'"
                )));
            }
        }
        if !found_objects {
            return Err(StateError::Invalid(
                "state is missing the [objects] section".to_owned(),
            ));
        }

        let state_format_version: u32 = parse_field(&fields, "state_format_version")?;
        if !(2..=STATE_FORMAT_VERSION).contains(&state_format_version) {
            return Err(StateError::Compatibility(format!(
                "unsupported state_format_version {state_format_version}; this build supports 2 through {STATE_FORMAT_VERSION}"
            )));
        }
        require_version(
            &fields,
            "simulation_version",
            if state_format_version >= 19 {
                SIMULATION_VERSION
            } else if state_format_version >= 18 {
                7
            } else if state_format_version >= 17 {
                6
            } else if state_format_version >= 16 {
                5
            } else if state_format_version >= 15 {
                4
            } else if state_format_version >= 14 {
                3
            } else if state_format_version >= 13 {
                2
            } else {
                1
            },
        )?;
        require_version(
            &fields,
            "generation_version",
            if state_format_version >= 19 {
                GENERATION_VERSION
            } else if state_format_version >= 18 {
                14
            } else if state_format_version >= 17 {
                13
            } else if state_format_version >= 15 {
                12
            } else if state_format_version >= 14 {
                11
            } else if state_format_version >= 13 {
                10
            } else if state_format_version >= 12 {
                9
            } else if state_format_version >= 11 {
                8
            } else if state_format_version >= 10 {
                7
            } else if state_format_version >= 9 {
                6
            } else if state_format_version >= 8 {
                5
            } else if state_format_version >= 7 {
                4
            } else if state_format_version >= 6 {
                3
            } else if state_format_version >= 4 {
                2
            } else {
                1
            },
        )?;
        require_version(&fields, "rng_version", RNG_VERSION)?;
        let _crate_version = required(&fields, "crate_version")?;

        let seed = parse_field(&fields, "seed")?;
        let id = UniverseId::new(parse_field(&fields, "universe_id")?);
        let expected_id = UniverseId::new(derive_object_seed(seed, 0, "universe"));
        if id != expected_id {
            return Err(StateError::Invalid(
                "universe ID does not match the stored seed".to_owned(),
            ));
        }

        let tick = parse_field(&fields, "tick")?;
        let elapsed = SimDuration::from_seconds(parse_field(&fields, "elapsed_seconds")?);
        let tick_duration = TickDuration::new(SimDuration::from_seconds(parse_field(
            &fields,
            "tick_duration_seconds",
        )?))
        .map_err(StateError::Time)?;
        let clock =
            SimulationClock::from_parts(tick, elapsed, tick_duration).map_err(StateError::Time)?;
        let generation = GenerationConfig {
            galaxies: if state_format_version >= 9 {
                parse_field(&fields, "galaxies")?
            } else {
                1
            },
            star_systems: parse_field(&fields, "star_systems")?,
            planets_per_system_min: parse_field(&fields, "planets_per_system_min")?,
            planets_per_system_max: parse_field(&fields, "planets_per_system_max")?,
            star_mass_solar_min: parse_float_field(&fields, "star_mass_solar_min")?,
            star_mass_solar_max: parse_float_field(&fields, "star_mass_solar_max")?,
            planet_mass_earth_min: parse_float_field(&fields, "planet_mass_earth_min")?,
            planet_mass_earth_max: parse_float_field(&fields, "planet_mass_earth_max")?,
            orbit_au_min: parse_float_field(&fields, "orbit_au_min")?,
            orbit_au_max: parse_float_field(&fields, "orbit_au_max")?,
            galaxy_radius_parsecs: parse_float_field(&fields, "galaxy_radius_parsecs")?,
        };
        let config = SimulationConfig {
            seed,
            tick_duration,
            generation,
        };
        config.validate().map_err(StateError::Configuration)?;
        let galaxies: Vec<_> = (0..generation.galaxies)
            .map(|index| {
                GalaxyCatalogEntry::generate(seed, index, generation.galaxy_radius_parsecs)
            })
            .collect::<Result<_, _>>()
            .map_err(StateError::Objects)?;
        let primary_galaxy_id = galaxies[0].id();

        let remaining: Vec<_> = lines.collect();
        let (object_lines, region_lines, event_lines, motion_lines, encounter_lines) =
            if state_format_version == 2 {
                (remaining.as_slice(), &[][..], &[][..], &[][..], &[][..])
            } else {
                let regions_index = remaining
                    .iter()
                    .position(|line| *line == "[regions]")
                    .ok_or_else(|| StateError::Invalid("state is missing [regions]".to_owned()))?;
                let events_index = remaining
                    .iter()
                    .position(|line| *line == "[events]")
                    .ok_or_else(|| StateError::Invalid("state is missing [events]".to_owned()))?;
                if events_index <= regions_index {
                    return Err(StateError::Invalid(
                        "state sections are out of order".to_owned(),
                    ));
                }
                if state_format_version >= 17 {
                    let motion_index = remaining
                        .iter()
                        .position(|line| *line == "[galaxy_motion]")
                        .ok_or_else(|| {
                            StateError::Invalid("state is missing [galaxy_motion]".to_owned())
                        })?;
                    let encounter_index = remaining
                        .iter()
                        .position(|line| *line == "[galaxy_encounters]")
                        .ok_or_else(|| {
                            StateError::Invalid("state is missing [galaxy_encounters]".to_owned())
                        })?;
                    if motion_index <= events_index || encounter_index <= motion_index {
                        return Err(StateError::Invalid(
                            "state sections are out of order".to_owned(),
                        ));
                    }
                    (
                        &remaining[..regions_index],
                        &remaining[regions_index + 1..events_index],
                        &remaining[events_index + 1..motion_index],
                        &remaining[motion_index + 1..encounter_index],
                        &remaining[encounter_index + 1..],
                    )
                } else {
                    (
                        &remaining[..regions_index],
                        &remaining[regions_index + 1..events_index],
                        &remaining[events_index + 1..],
                        &[][..],
                        &[][..],
                    )
                }
            };

        let galaxy_motion = if state_format_version >= 17 {
            decode_galaxy_motion(motion_lines, &galaxies, elapsed)?
        } else {
            galaxies
                .iter()
                .map(|galaxy| (galaxy.id(), GalaxyMotion::from_catalog(*galaxy, elapsed)))
                .collect()
        };
        let galaxy_encounters = if state_format_version >= 17 {
            decode_galaxy_encounters(encounter_lines, &galaxies, seed, tick, tick_duration)?
        } else {
            Vec::new()
        };

        let legacy_object_lines: Vec<_> = if state_format_version < 12 {
            object_lines
                .iter()
                .copied()
                .filter(|line| {
                    !line.starts_with("X|") && !line.starts_with("K|") && !line.starts_with("Q|")
                })
                .collect()
        } else {
            object_lines.to_vec()
        };
        let mut objects =
            ObjectRegistry::decode(legacy_object_lines.into_iter()).map_err(StateError::Objects)?;
        if state_format_version < 18 && elapsed != SimDuration::ZERO {
            objects
                .advance_orbits_to(elapsed)
                .map_err(StateError::Objects)?;
        }
        if state_format_version < 19 {
            objects
                .refresh_main_sequence_luminosities()
                .map_err(StateError::Time)?;
            objects
                .refresh_planet_temperatures()
                .map_err(StateError::Time)?;
        }
        if !objects
            .galaxies()
            .all(|galaxy| galaxy.id() == primary_galaxy_id)
        {
            return Err(StateError::Invalid(
                "origin objects do not belong to the primary galaxy".to_owned(),
            ));
        }
        validate_registry(&objects, elapsed, generation)?;

        let mut known_regions = BTreeMap::new();
        let mut materialized_flags = BTreeMap::new();
        if state_format_version == 2 {
            let origin = GalaxyRegion::new(primary_galaxy_id, RegionCoord::ORIGIN);
            known_regions.insert(origin, tick);
            materialized_flags.insert(origin, true);
        } else {
            for line in region_lines {
                let fields: Vec<_> = line.split('|').collect();
                let expected_fields = if state_format_version >= 9 { 7 } else { 6 };
                if fields.len() != expected_fields || fields[0] != "R" {
                    return Err(StateError::Invalid(format!(
                        "malformed region record '{line}'"
                    )));
                }
                let offset = usize::from(state_format_version >= 9);
                let galaxy_id = if state_format_version >= 9 {
                    GalaxyId::new(parse_number(fields[1], "region galaxy ID")?)
                } else {
                    primary_galaxy_id
                };
                if !galaxies.iter().any(|galaxy| galaxy.id() == galaxy_id) {
                    return Err(StateError::Invalid(format!(
                        "region references unknown galaxy {galaxy_id}"
                    )));
                }
                let coordinate = RegionCoord::new(
                    parse_number(fields[1 + offset], "region x")?,
                    parse_number(fields[2 + offset], "region y")?,
                    parse_number(fields[3 + offset], "region z")?,
                );
                let key = GalaxyRegion::new(galaxy_id, coordinate);
                let last_updated_tick = parse_number(fields[4 + offset], "region update tick")?;
                if last_updated_tick > tick {
                    return Err(StateError::Invalid(format!(
                        "region {coordinate} was updated after the universal clock"
                    )));
                }
                let materialized = match fields[5 + offset] {
                    "0" => false,
                    "1" => true,
                    _ => {
                        return Err(StateError::Invalid(format!(
                            "region {coordinate} has an invalid materialized flag"
                        )));
                    }
                };
                if known_regions.insert(key, last_updated_tick).is_some() {
                    return Err(StateError::Invalid(format!(
                        "duplicate region {coordinate}"
                    )));
                }
                materialized_flags.insert(key, materialized);
            }
            let primary_origin = GalaxyRegion::new(primary_galaxy_id, RegionCoord::ORIGIN);
            if known_regions.get(&primary_origin) != Some(&tick)
                || materialized_flags.get(&primary_origin) != Some(&true)
            {
                return Err(StateError::Invalid(
                    "the origin region must be materialized at the universal tick".to_owned(),
                ));
            }
        }

        let mut regions = BTreeMap::new();
        let primary_origin = GalaxyRegion::new(primary_galaxy_id, RegionCoord::ORIGIN);
        for (key, materialized) in &materialized_flags {
            if *key == primary_origin || !materialized {
                continue;
            }
            if known_regions[key] != tick {
                return Err(StateError::Invalid(format!(
                    "materialized region {key} is not current"
                )));
            }
            let galaxy = galaxies
                .iter()
                .find(|galaxy| galaxy.id() == key.galaxy_id())
                .copied()
                .expect("validated region galaxy exists");
            let region_seed = scoped_region_seed(seed, key.galaxy_id(), key.coordinate());
            let mut region_objects =
                ObjectRegistry::generate_in_galaxy(region_seed, generation, galaxy)
                    .map_err(StateError::Objects)?;
            region_objects
                .form_stars_between_with_encounters(
                    region_seed,
                    generation,
                    SimDuration::ZERO,
                    elapsed,
                    tick_duration.duration(),
                    &galaxy_encounters,
                )
                .map_err(StateError::Objects)?;
            region_objects
                .advance_stars_to(elapsed, tick_duration.duration())
                .map_err(StateError::Time)?;
            region_objects
                .advance_planetary_formation_to(region_seed, generation, elapsed)
                .map_err(StateError::Objects)?;
            region_objects
                .advance_orbits_to(elapsed)
                .map_err(StateError::Objects)?;
            validate_registry(&region_objects, elapsed, generation)?;
            regions.insert(
                *key,
                MaterializedRegion::new(key.coordinate(), tick, region_objects),
            );
        }

        let mut events = Vec::new();
        let mut event_ids = BTreeSet::new();
        for line in event_lines {
            let fields: Vec<_> = line.split('|').collect();
            let valid_fields = if state_format_version >= 9 {
                fields.len() == 14
            } else {
                matches!(fields.len(), 9..=11 | 13)
            };
            if !valid_fields || fields[0] != "E" {
                return Err(StateError::Invalid(format!(
                    "malformed event record '{line}'"
                )));
            }
            let event_id = parse_number(fields[1], "event ID")?;
            let event_tick = parse_number(fields[2], "event tick")?;
            let offset = usize::from(state_format_version >= 9);
            let galaxy_id = if state_format_version >= 9 {
                GalaxyId::new(parse_number(fields[3], "event galaxy ID")?)
            } else {
                primary_galaxy_id
            };
            let coordinate = RegionCoord::new(
                parse_number(fields[3 + offset], "event region x")?,
                parse_number(fields[4 + offset], "event region y")?,
                parse_number(fields[5 + offset], "event region z")?,
            );
            let key = GalaxyRegion::new(galaxy_id, coordinate);
            let star_id =
                crate::id::StarId::new(parse_number(fields[6 + offset], "event star ID")?);
            let from = crate::objects::StellarStage::parse_code(fields[7 + offset])
                .map_err(StateError::Objects)?;
            let to = crate::objects::StellarStage::parse_code(fields[8 + offset])
                .map_err(StateError::Objects)?;
            let impacts = if fields.len() >= 10 + offset {
                parse_impacts(fields[9 + offset])?
            } else {
                Vec::new()
            };
            let has_details = if state_format_version >= 9 {
                true
            } else {
                matches!(fields.len(), 11 | 13)
            };
            let natal_kick_meters_per_second = if has_details {
                parse_number(fields[10 + offset], "natal kick")?
            } else {
                0
            };
            let has_interaction = state_format_version >= 9 || fields.len() == 13;
            let interaction = if has_interaction && fields[11 + offset] != "-" {
                Some(
                    crate::objects::StellarInteractionKind::parse_code(fields[11 + offset])
                        .map_err(StateError::Objects)?,
                )
            } else {
                None
            };
            let related_star_id = if has_interaction && fields[12 + offset] != "-" {
                Some(crate::id::StarId::new(parse_number(
                    fields[12 + offset],
                    "related star ID",
                )?))
            } else {
                None
            };
            if matches!(
                interaction,
                Some(
                    crate::objects::StellarInteractionKind::TypeIaSupernova
                        | crate::objects::StellarInteractionKind::Kilonova
                        | crate::objects::StellarInteractionKind::BlackHoleMerger
                        | crate::objects::StellarInteractionKind::CompactObjectMerger
                )
            ) && related_star_id.is_none()
            {
                return Err(StateError::Invalid(
                    "compact-binary event is missing its companion ID".to_owned(),
                ));
            }
            if event_tick > tick || !known_regions.contains_key(&key) {
                return Err(StateError::Invalid(
                    "event lies outside the known simulation history".to_owned(),
                ));
            }
            let region_seed = scoped_region_seed(seed, galaxy_id, coordinate);
            let namespace = event_namespace(event_tick, to, interaction, related_star_id);
            let expected_id = derive_object_seed(region_seed, star_id.value(), &namespace);
            if event_id != expected_id || !event_ids.insert(event_id) {
                return Err(StateError::Invalid(
                    "event ID is invalid or duplicated".to_owned(),
                ));
            }
            events.push(AstroEvent::new(
                event_id,
                event_tick,
                AstroEventLocation {
                    galaxy_id,
                    region: coordinate,
                    star_id,
                },
                from,
                to,
                AstroEventDetails {
                    impacts,
                    natal_kick_meters_per_second,
                    interaction,
                    related_star_id,
                },
            ));
        }
        events.sort_by_key(|event| (event.tick(), event.id()));

        Ok(Self {
            id,
            config,
            clock,
            galaxies,
            galaxy_motion,
            galaxy_encounters,
            objects,
            known_regions,
            regions,
            events,
        })
    }

    fn decode_v1<'a>(lines: impl Iterator<Item = &'a str>) -> Result<Self, StateError> {
        let mut fields = BTreeMap::new();
        for line in lines {
            let (key, value) = line
                .split_once('=')
                .ok_or_else(|| StateError::Invalid(format!("malformed state line '{line}'")))?;
            if fields.insert(key, value).is_some() {
                return Err(StateError::Invalid(format!(
                    "duplicate state field '{key}'"
                )));
            }
        }
        let _version = required(&fields, "version")?;
        let seed = parse_field(&fields, "seed")?;
        let stored_id: u64 = parse_field(&fields, "universe_id")?;
        let expected_id = derive_object_seed(seed, 0, "universe");
        if stored_id != expected_id {
            return Err(StateError::Invalid(
                "universe ID does not match the stored seed".to_owned(),
            ));
        }
        let tick = parse_field(&fields, "tick")?;
        let elapsed = SimDuration::from_seconds(parse_field(&fields, "elapsed_seconds")?);
        let tick_duration = TickDuration::new(SimDuration::from_seconds(parse_field(
            &fields,
            "tick_duration_seconds",
        )?))
        .map_err(StateError::Time)?;
        let clock =
            SimulationClock::from_parts(tick, elapsed, tick_duration).map_err(StateError::Time)?;
        let config = SimulationConfig::new(seed, tick_duration);
        let mut universe = Self::new(config).map_err(|error| match error {
            UniverseCreateError::Configuration(error) => StateError::Configuration(error),
            UniverseCreateError::Objects(error) => StateError::Objects(error),
        })?;
        let transitions = universe
            .objects
            .advance_stars_to(elapsed, tick_duration.duration())
            .map_err(StateError::Time)?;
        universe
            .objects
            .advance_planetary_formation_to(universe.seed(), universe.config().generation, elapsed)
            .map_err(StateError::Objects)?;
        universe
            .objects
            .advance_orbits_to(elapsed)
            .map_err(StateError::Objects)?;
        universe.clock = clock;
        for (galaxy, motion) in &mut universe.galaxy_motion {
            *motion = GalaxyMotion::from_catalog(
                *universe
                    .galaxies
                    .iter()
                    .find(|entry| entry.id() == *galaxy)
                    .expect("legacy galaxy exists"),
                elapsed,
            );
        }
        let origin = GalaxyRegion::new(universe.primary_galaxy_id(), RegionCoord::ORIGIN);
        universe.known_regions.insert(origin, tick);
        universe.record_transitions(origin, transitions);
        Ok(universe)
    }
}

fn parse_float_bits(value: &str, label: &str) -> Result<f64, StateError> {
    let bits = u64::from_str_radix(value, 16)
        .map_err(|_| StateError::Invalid(format!("{label} has invalid float bits")))?;
    let number = f64::from_bits(bits);
    if !number.is_finite() {
        return Err(StateError::Invalid(format!("{label} must be finite")));
    }
    Ok(number)
}

fn decode_galaxy_motion(
    lines: &[&str],
    catalog: &[GalaxyCatalogEntry],
    elapsed: SimDuration,
) -> Result<BTreeMap<GalaxyId, GalaxyMotion>, StateError> {
    let known: BTreeSet<_> = catalog.iter().map(|galaxy| galaxy.id()).collect();
    let mut motions = BTreeMap::new();
    for line in lines {
        let fields: Vec<_> = line.split('|').collect();
        if fields.len() != 12 || fields[0] != "D" {
            return Err(StateError::Invalid(format!(
                "malformed galaxy motion record '{line}'"
            )));
        }
        let id = GalaxyId::new(parse_number(fields[1], "galaxy motion ID")?);
        if !known.contains(&id) || motions.contains_key(&id) {
            return Err(StateError::Invalid(
                "galaxy motion ID is unknown or duplicated".to_owned(),
            ));
        }
        let position_at_epoch = Position::from_meters(
            parse_float_bits(fields[2], "galaxy x")?,
            parse_float_bits(fields[3], "galaxy y")?,
            parse_float_bits(fields[4], "galaxy z")?,
        )
        .map_err(|error| StateError::Invalid(error.to_string()))?;
        let velocity = Velocity::from_meters_per_second(
            parse_float_bits(fields[5], "galaxy vx")?,
            parse_float_bits(fields[6], "galaxy vy")?,
            parse_float_bits(fields[7], "galaxy vz")?,
        )
        .map_err(|error| StateError::Invalid(error.to_string()))?;
        let epoch = SimDuration::from_seconds(parse_number(fields[8], "galaxy epoch")?);
        if epoch > elapsed {
            return Err(StateError::Invalid(
                "galaxy motion epoch exceeds the clock".to_owned(),
            ));
        }
        let morphology =
            crate::objects::GalaxyMorphology::parse_code(fields[9]).map_err(StateError::Objects)?;
        let radius_parsecs = parse_float_bits(fields[10], "galaxy radius")?;
        if radius_parsecs <= 0.0 {
            return Err(StateError::Invalid(
                "galaxy radius must be positive".to_owned(),
            ));
        }
        let merged_into = if fields[11] == "-" {
            None
        } else {
            let survivor = GalaxyId::new(parse_number(fields[11], "merged galaxy ID")?);
            if survivor == id || !known.contains(&survivor) {
                return Err(StateError::Invalid(
                    "galaxy merger target is invalid".to_owned(),
                ));
            }
            Some(survivor)
        };
        let motion = GalaxyMotion {
            position_at_epoch,
            velocity,
            epoch,
            morphology,
            radius_parsecs,
            merged_into,
        };
        motion.position_at(elapsed).map_err(StateError::Time)?;
        motions.insert(id, motion);
    }
    if motions.len() != catalog.len() {
        return Err(StateError::Invalid(
            "galaxy motion catalog is incomplete".to_owned(),
        ));
    }
    for id in &known {
        let mut current = *id;
        for _ in 0..catalog.len() {
            let Some(next) = motions[&current].merged_into else {
                break;
            };
            current = next;
        }
        if motions[&current].merged_into.is_some() {
            return Err(StateError::Invalid(
                "galaxy merger targets form a cycle".to_owned(),
            ));
        }
    }
    Ok(motions)
}

fn decode_galaxy_encounters(
    lines: &[&str],
    catalog: &[GalaxyCatalogEntry],
    seed: u64,
    tick: u64,
    tick_duration: TickDuration,
) -> Result<Vec<GalaxyEncounter>, StateError> {
    let known: BTreeSet<_> = catalog.iter().map(|galaxy| galaxy.id()).collect();
    let mut events = Vec::new();
    let mut ids = BTreeSet::new();
    let mut pairs = BTreeSet::new();
    for line in lines {
        let fields: Vec<_> = line.split('|').collect();
        if fields.len() != 7 || fields[0] != "C" {
            return Err(StateError::Invalid(format!(
                "malformed galaxy encounter record '{line}'"
            )));
        }
        let id: u64 = parse_number(fields[1], "galaxy encounter ID")?;
        let seconds: u128 = parse_number(fields[2], "galaxy encounter time")?;
        let event_tick: u64 = parse_number(fields[3], "galaxy encounter tick")?;
        let first = GalaxyId::new(parse_number(fields[4], "first galaxy ID")?);
        let second = GalaxyId::new(parse_number(fields[5], "second galaxy ID")?);
        let kind = GalaxyEncounterKind::parse_code(fields[6])
            .ok_or_else(|| StateError::Invalid("unknown galaxy encounter kind".to_owned()))?;
        let expected_tick = u64::try_from(
            seconds
                .checked_add(tick_duration.duration().seconds() - 1)
                .ok_or(StateError::Time(TimeError::Overflow))?
                / tick_duration.duration().seconds(),
        )
        .map_err(|_| StateError::Time(TimeError::Overflow))?;
        if first >= second
            || !known.contains(&first)
            || !known.contains(&second)
            || event_tick != expected_tick
            || event_tick > tick
            || id != encounter_id(seed, first, second, kind, seconds)
            || !ids.insert(id)
            || !pairs.insert((first, second, kind))
        {
            return Err(StateError::Invalid(
                "galaxy encounter is invalid or duplicated".to_owned(),
            ));
        }
        events.push(GalaxyEncounter::new(
            id,
            SimDuration::from_seconds(seconds),
            event_tick,
            first,
            second,
            kind,
        ));
    }
    events.sort_by_key(|event| (event.occurred_at(), event.id()));
    Ok(events)
}

fn atomic_write(path: &Path, contents: &[u8]) -> Result<(), StateError> {
    let parent = path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let file_name = path.file_name().ok_or_else(|| {
        StateError::Invalid(format!("'{}' is not a state-file path", path.display()))
    })?;

    let mut temp_path = None;
    let mut temp_file = None;
    for _ in 0..16 {
        let sequence = TEMP_FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let mut temp_name = OsString::from(".");
        temp_name.push(file_name);
        temp_name.push(format!(".{}.{}.tmp", std::process::id(), sequence));
        let candidate = parent.join(temp_name);
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => {
                temp_path = Some(candidate);
                temp_file = Some(file);
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(StateError::Io(error)),
        }
    }

    let temp_path = temp_path.ok_or_else(|| {
        StateError::Invalid("could not reserve a temporary checkpoint file".to_owned())
    })?;
    let mut temp_file = temp_file.expect("temporary path and file are set together");
    let write_result = temp_file
        .write_all(contents)
        .and_then(|()| temp_file.sync_all());
    if let Err(error) = write_result {
        drop(temp_file);
        let _ = fs::remove_file(&temp_path);
        return Err(StateError::Io(error));
    }
    drop(temp_file);

    if let Err(error) = fs::rename(&temp_path, path) {
        let _ = fs::remove_file(&temp_path);
        return Err(StateError::Io(error));
    }
    File::open(parent)
        .and_then(|directory| directory.sync_all())
        .map_err(StateError::Io)
}

fn require_version(
    fields: &BTreeMap<&str, &str>,
    key: &str,
    supported: u32,
) -> Result<(), StateError> {
    let found: u32 = parse_field(fields, key)?;
    if found == supported {
        Ok(())
    } else {
        Err(StateError::Compatibility(format!(
            "unsupported {key} {found}; this build supports {supported}"
        )))
    }
}

fn encode_impacts(impacts: &[PlanetImpact]) -> String {
    impacts
        .iter()
        .map(|impact| format!("{}:{}", impact.planet_id.value(), impact.kind.code()))
        .collect::<Vec<_>>()
        .join(",")
}

fn parse_impacts(value: &str) -> Result<Vec<PlanetImpact>, StateError> {
    if value.is_empty() {
        return Ok(Vec::new());
    }
    let mut ids = BTreeSet::new();
    value
        .split(',')
        .map(|record| {
            let (id, kind) = record.split_once(':').ok_or_else(|| {
                StateError::Invalid(format!("malformed planetary impact '{record}'"))
            })?;
            let planet_id = crate::id::PlanetId::new(parse_number(id, "impact planet ID")?);
            if !ids.insert(planet_id) {
                return Err(StateError::Invalid(format!(
                    "duplicate planetary impact for {planet_id}"
                )));
            }
            Ok(PlanetImpact {
                planet_id,
                kind: PlanetImpactKind::parse_code(kind).map_err(StateError::Objects)?,
            })
        })
        .collect()
}

fn event_namespace(
    tick: u64,
    stage: crate::objects::StellarStage,
    interaction: Option<crate::objects::StellarInteractionKind>,
    related_star_id: Option<crate::id::StarId>,
) -> String {
    if let Some(interaction) = interaction {
        if matches!(
            interaction,
            crate::objects::StellarInteractionKind::TypeIaSupernova
                | crate::objects::StellarInteractionKind::Kilonova
                | crate::objects::StellarInteractionKind::BlackHoleMerger
                | crate::objects::StellarInteractionKind::CompactObjectMerger
        ) {
            return format!(
                "compact-binary-event-v1-{tick}-{}-{}-{}",
                stage.code(),
                interaction.code(),
                related_star_id
                    .expect("compact merger has related star")
                    .value()
            );
        }
        format!(
            "stellar-event-v2-{tick}-{}-{}",
            stage.code(),
            interaction.code()
        )
    } else {
        format!("stellar-event-v1-{tick}-{}", stage.code())
    }
}

fn scoped_region_seed(seed: u64, galaxy_id: GalaxyId, coordinate: RegionCoord) -> u64 {
    let primary_galaxy_id = GalaxyId::new(derive_object_seed(seed, 0, "galaxy-id-v1"));
    if galaxy_id == primary_galaxy_id {
        return derive_region_seed(seed, coordinate.x, coordinate.y, coordinate.z);
    }
    let galaxy_seed = derive_object_seed(seed, galaxy_id.value(), "galaxy-region-seed-v1");
    derive_region_seed(galaxy_seed, coordinate.x, coordinate.y, coordinate.z)
}

fn required<'a>(fields: &'a BTreeMap<&str, &str>, key: &str) -> Result<&'a str, StateError> {
    fields
        .get(key)
        .copied()
        .ok_or_else(|| StateError::Invalid(format!("missing state field '{key}'")))
}

fn parse_field<T>(fields: &BTreeMap<&str, &str>, key: &str) -> Result<T, StateError>
where
    T: std::str::FromStr,
{
    required(fields, key)?.parse().map_err(|_| {
        StateError::Invalid(format!("state field '{key}' has an invalid numeric value"))
    })
}

fn parse_float_field(fields: &BTreeMap<&str, &str>, key: &str) -> Result<f64, StateError> {
    u64::from_str_radix(required(fields, key)?, 16)
        .map(f64::from_bits)
        .map_err(|_| StateError::Invalid(format!("state field '{key}' has invalid float bits")))
}

fn parse_number<T>(value: &str, field: &str) -> Result<T, StateError>
where
    T: std::str::FromStr,
{
    value
        .parse()
        .map_err(|_| StateError::Invalid(format!("invalid {field}")))
}

fn validate_registry(
    objects: &ObjectRegistry,
    elapsed: SimDuration,
    generation: GenerationConfig,
) -> Result<(), StateError> {
    if objects.stars().any(|star| {
        star.epoch_age()
            .checked_add(SimDuration::from_seconds(
                elapsed.seconds().saturating_sub(star.formed_at().seconds()),
            ))
            .is_none_or(|expected| star.age() != expected)
    }) {
        return Err(StateError::Invalid(
            "star ages do not match their epoch age plus universal elapsed time".to_owned(),
        ));
    }
    let counts = objects.counts();
    if counts.galaxies != 1
        || counts.star_systems < generation.star_systems as usize
        || counts.stars < generation.star_systems as usize
        || counts.stars > counts.star_systems * 3
    {
        return Err(StateError::Invalid(
            "object counts do not match the stored generation configuration".to_owned(),
        ));
    }
    if objects.star_systems().any(|system| {
        let primary = objects.star(system.star_id());
        let exempt_from_initial_minimum = primary.is_some_and(|star| {
            star.stage() == crate::objects::StellarStage::Disrupted
                || star.formed_at() > SimDuration::ZERO
        });
        (!exempt_from_initial_minimum
            && system.planet_ids().len() < usize::from(generation.planets_per_system_min))
            || system.planet_ids().len() > usize::from(generation.planets_per_system_max)
    }) {
        return Err(StateError::Invalid(
            "planet counts do not match the stored generation configuration".to_owned(),
        ));
    }
    Ok(())
}

#[derive(Debug)]
pub enum SimulationAdvanceError {
    Time(TimeError),
    Objects(ObjectError),
}

impl fmt::Display for SimulationAdvanceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Time(error) => write!(formatter, "could not advance simulation time: {error}"),
            Self::Objects(error) => {
                write!(formatter, "could not advance simulated objects: {error}")
            }
        }
    }
}

impl Error for SimulationAdvanceError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Time(error) => Some(error),
            Self::Objects(error) => Some(error),
        }
    }
}

#[derive(Debug)]
pub enum RegionError {
    CannotUnloadOrigin,
    UnknownGalaxy(GalaxyId),
    MergedGalaxy {
        galaxy_id: GalaxyId,
        survivor: GalaxyId,
    },
    Generation(ObjectError),
    Time(TimeError),
}

impl fmt::Display for RegionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CannotUnloadOrigin => formatter.write_str("the origin region cannot be unloaded"),
            Self::UnknownGalaxy(id) => write!(formatter, "unknown galaxy {id}"),
            Self::MergedGalaxy {
                galaxy_id,
                survivor,
            } => write!(
                formatter,
                "galaxy {galaxy_id} merged into {survivor}; materialize a region in {survivor}"
            ),
            Self::Generation(error) => write!(formatter, "could not generate region: {error}"),
            Self::Time(error) => write!(formatter, "could not catch region up: {error}"),
        }
    }
}

impl Error for RegionError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Generation(error) => Some(error),
            Self::Time(error) => Some(error),
            Self::CannotUnloadOrigin => None,
            Self::UnknownGalaxy(_) | Self::MergedGalaxy { .. } => None,
        }
    }
}

#[derive(Debug)]
pub enum UniverseCreateError {
    Configuration(ConfigValidationError),
    Objects(ObjectError),
}

impl fmt::Display for UniverseCreateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Configuration(error) => write!(formatter, "invalid configuration: {error}"),
            Self::Objects(error) => write!(formatter, "could not generate objects: {error}"),
        }
    }
}

impl Error for UniverseCreateError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Configuration(error) => Some(error),
            Self::Objects(error) => Some(error),
        }
    }
}

#[derive(Debug)]
pub enum StateError {
    Io(std::io::Error),
    Compatibility(String),
    Configuration(ConfigValidationError),
    Invalid(String),
    Objects(ObjectError),
    Time(TimeError),
}

impl fmt::Display for StateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "could not access simulation state: {error}"),
            Self::Compatibility(message) => write!(formatter, "incompatible state: {message}"),
            Self::Configuration(error) => {
                write!(formatter, "invalid stored configuration: {error}")
            }
            Self::Invalid(message) => write!(formatter, "invalid simulation state: {message}"),
            Self::Objects(error) => write!(formatter, "invalid object registry: {error}"),
            Self::Time(error) => write!(formatter, "invalid simulation clock: {error}"),
        }
    }
}

impl Error for StateError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Configuration(error) => Some(error),
            Self::Objects(error) => Some(error),
            Self::Time(error) => Some(error),
            Self::Compatibility(_) | Self::Invalid(_) => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_config() -> SimulationConfig {
        SimulationConfig::new(918_273_645, "1 kyr".parse().unwrap())
    }

    fn legacy_state_without_galaxy_sections(universe: &Universe) -> String {
        universe
            .encode()
            .split("[galaxy_motion]\n")
            .next()
            .unwrap()
            .to_owned()
    }

    #[test]
    fn version_fifteen_state_loads_under_compact_binary_rules() {
        let universe = Universe::new(test_config()).unwrap();
        let legacy = legacy_state_without_galaxy_sections(&universe)
            .replace("state_format_version=19", "state_format_version=15")
            .replace("simulation_version=8", "simulation_version=4")
            .replace("generation_version=15", "generation_version=12");
        assert_eq!(Universe::decode(&legacy).unwrap(), universe);
    }

    #[test]
    fn version_sixteen_state_starts_galaxy_motion_at_saved_clock() {
        let mut config = test_config();
        config.generation.galaxies = 2;
        let mut universe = Universe::new(config).unwrap();
        universe.advance_ticks(10_000).unwrap();
        let legacy = legacy_state_without_galaxy_sections(&universe)
            .replace("state_format_version=19", "state_format_version=16")
            .replace("simulation_version=8", "simulation_version=5")
            .replace("generation_version=15", "generation_version=12");
        let restored = Universe::decode(&legacy).unwrap();
        let original =
            GalaxyCatalogEntry::generate(config.seed, 1, config.generation.galaxy_radius_parsecs)
                .unwrap();
        assert_eq!(
            restored.galaxies().nth(1).unwrap().position(),
            original.position()
        );
        assert!(restored.galaxy_encounters().is_empty());
    }

    #[test]
    fn version_seventeen_state_gains_cluster_member_motion() {
        let mut universe = Universe::new(test_config()).unwrap();
        universe.advance_ticks(10_000).unwrap();
        let system = universe.objects().star_systems().next().unwrap();
        let cluster = universe
            .objects()
            .star_cluster(system.cluster_id().unwrap())
            .unwrap();
        let expected = system.position();
        let [cx, cy, cz] = cluster.position().meters();
        let system_prefix = format!("Y|{}|", system.id().value());
        let legacy = universe
            .encode()
            .lines()
            .map(|line| {
                if line.starts_with(&system_prefix) {
                    let mut fields: Vec<_> = line.split('|').map(str::to_owned).collect();
                    for (field, coordinate) in [cx, cy, cz].into_iter().enumerate() {
                        fields[field + 3] = format!("{:016x}", coordinate.to_bits());
                    }
                    fields.join("|")
                } else {
                    line.replace("state_format_version=19", "state_format_version=17")
                        .replace("simulation_version=8", "simulation_version=6")
                        .replace("generation_version=15", "generation_version=13")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let restored = Universe::decode(&legacy).unwrap();
        assert_eq!(
            restored
                .objects()
                .star_system(system.id())
                .unwrap()
                .position(),
            expected
        );
        assert!(restored.encode().contains("state_format_version=19"));
    }

    #[test]
    fn version_eighteen_state_refreshes_main_sequence_luminosity() {
        let universe = Universe::new(test_config()).unwrap();
        let star = universe
            .objects()
            .stars()
            .find(|star| star.stage() == crate::objects::StellarStage::MainSequence)
            .unwrap();
        let star_prefix = format!("S|{}|", star.id().value());
        let old = universe
            .encode()
            .lines()
            .map(|line| {
                if line.starts_with(&star_prefix) {
                    let mut fields: Vec<_> = line.split('|').map(str::to_owned).collect();
                    fields[8] = format!("{:016x}", 0.0_f64.to_bits());
                    fields.join("|")
                } else {
                    line.replace("state_format_version=19", "state_format_version=18")
                        .replace("simulation_version=8", "simulation_version=7")
                        .replace("generation_version=15", "generation_version=14")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let restored = Universe::decode(&old).unwrap();
        assert_eq!(
            restored.objects().star(star.id()).unwrap().luminosity(),
            star.luminosity()
        );
        assert!(restored.encode().contains("state_format_version=19"));
    }

    #[test]
    fn compact_merger_events_round_trip_with_distinct_companions() {
        let mut universe = Universe::new(test_config()).unwrap();
        let ids: Vec<_> = universe
            .objects()
            .stars()
            .take(3)
            .map(|star| star.id())
            .collect();
        let location = GalaxyRegion::new(universe.primary_galaxy_id(), RegionCoord::ORIGIN);
        let transitions = ids[1..]
            .iter()
            .map(|companion_id| StellarTransition {
                tick: 0,
                star_id: ids[0],
                from: crate::objects::StellarStage::BlackHole,
                to: crate::objects::StellarStage::BlackHole,
                impacts: Vec::new(),
                natal_kick_meters_per_second: 0,
                interaction: Some(crate::objects::StellarInteractionKind::BlackHoleMerger),
                related_star_id: Some(*companion_id),
            })
            .collect();
        universe.record_transitions(location, transitions);
        assert_eq!(universe.events().len(), 2);
        let restored = Universe::decode(&universe.encode()).unwrap();
        assert_eq!(restored.events(), universe.events());
        assert!(
            restored
                .events()
                .iter()
                .all(|event| event.kind() == crate::AstroEventKind::BlackHoleMerger)
        );
    }

    #[test]
    fn galaxy_encounters_and_mergers_survive_checkpoints() {
        let mut config = SimulationConfig::new(42, "1 year".parse().unwrap());
        config.generation.galaxies = 8;
        config.generation.star_systems = 1;
        config.generation.planets_per_system_min = 0;
        config.generation.planets_per_system_max = 0;
        let mut whole = Universe::new(config).unwrap();
        let mut split = whole.clone();
        whole.advance_ticks(10_000_000_000).unwrap();
        split.advance_ticks(5_000_000_000).unwrap();
        split = Universe::decode(&split.encode()).unwrap();
        split.advance_ticks(5_000_000_000).unwrap();
        assert_eq!(split.galaxy_encounters(), whole.galaxy_encounters());
        assert_eq!(split.galaxy_motion, whole.galaxy_motion);
        assert_eq!(split.active_galaxy_count(), whole.active_galaxy_count());
        assert!(
            whole
                .galaxy_encounters()
                .iter()
                .any(|event| event.kind() == GalaxyEncounterKind::ClosePass)
        );
        assert!(
            whole
                .galaxy_encounters()
                .iter()
                .any(|event| event.kind() == GalaxyEncounterKind::Merger)
        );
        assert!(whole.active_galaxy_count() < 8);
        let absorbed = whole
            .galaxies()
            .find(|galaxy| galaxy.merged_into().is_some())
            .unwrap();
        let survivor = whole.surviving_galaxy_id(absorbed.id()).unwrap();
        assert_ne!(absorbed.id(), survivor);
        assert!(matches!(
            whole.materialize_region_in(absorbed.id(), RegionCoord::new(1, 0, 0)),
            Err(RegionError::MergedGalaxy { .. })
        ));
        assert_eq!(
            Universe::decode(&whole.encode())
                .unwrap()
                .galaxy_encounters(),
            whole.galaxy_encounters()
        );
    }

    #[test]
    fn encountered_progenitor_region_can_be_rematerialized() {
        let mut config = SimulationConfig::new(42, "1 year".parse().unwrap());
        config.generation.galaxies = 8;
        config.generation.star_systems = 1;
        config.generation.planets_per_system_min = 0;
        config.generation.planets_per_system_max = 0;
        let mut forecast = Universe::new(config).unwrap();
        forecast.advance_ticks(10_000_000_000).unwrap();
        let merger = forecast
            .galaxy_encounters()
            .iter()
            .find(|event| event.kind() == GalaxyEncounterKind::Merger)
            .unwrap();
        let progenitor = if merger.first() == forecast.primary_galaxy_id() {
            merger.second()
        } else {
            merger.first()
        };
        let mut universe = Universe::new(config).unwrap();
        universe
            .materialize_region_in(progenitor, RegionCoord::ORIGIN)
            .unwrap();
        universe.advance_ticks(10_000_000_000).unwrap();
        let original = universe
            .region_objects_in(progenitor, RegionCoord::ORIGIN)
            .unwrap()
            .clone();
        let restored = Universe::decode(&universe.encode()).unwrap();
        assert_eq!(
            restored.region_objects_in(progenitor, RegionCoord::ORIGIN),
            Some(&original)
        );
        universe
            .unload_region_in(progenitor, RegionCoord::ORIGIN)
            .unwrap();
        universe
            .materialize_region_in(progenitor, RegionCoord::ORIGIN)
            .unwrap();
        assert_eq!(
            universe.region_objects_in(progenitor, RegionCoord::ORIGIN),
            Some(&original)
        );
    }

    #[test]
    fn identical_configuration_produces_identical_universes() {
        let mut first = Universe::new(test_config()).unwrap();
        let mut second = Universe::new(test_config()).unwrap();
        first.advance_ticks(10_000).unwrap();
        second.advance_ticks(10_000).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn multiple_galaxies_have_independent_lazy_regions() {
        let mut config = test_config();
        config.generation.galaxies = 3;
        config.generation.star_systems = 2;
        config.generation.star_mass_solar_min = 0.1;
        config.generation.star_mass_solar_max = 0.2;
        let mut universe = Universe::new(config).unwrap();
        let galaxies: Vec<_> = universe.galaxies().collect();
        assert_eq!(galaxies.len(), 3);
        assert_eq!(universe.objects().counts().star_systems, 2);
        assert_ne!(galaxies[0].id(), galaxies[1].id());
        assert_ne!(galaxies[0].position(), galaxies[1].position());

        let second = galaxies[1].id();
        assert!(
            universe
                .materialize_region_in(second, RegionCoord::ORIGIN)
                .unwrap()
        );
        let second_objects = universe
            .region_objects_in(second, RegionCoord::ORIGIN)
            .unwrap();
        assert_eq!(second_objects.counts().star_systems, 2);
        assert!(
            second_objects
                .galaxies()
                .all(|galaxy| galaxy.id() == second)
        );
        assert_eq!(universe.galaxy_regions().count(), 2);
        let inventories: Vec<_> = universe.matter_inventories().collect();
        assert_eq!(inventories.len(), 2);
        assert!(
            inventories
                .iter()
                .all(|(_, inventory)| inventory.tracked_solar_masses() > 0.0)
        );
        universe.advance_ticks(10_000).unwrap();
        assert!(
            universe
                .events()
                .iter()
                .any(|event| event.galaxy_id() == second)
        );

        let decoded = Universe::decode(&universe.encode()).unwrap();
        assert_eq!(decoded, universe);
    }

    #[test]
    fn bounded_galaxy_batches_preserve_deterministic_results() {
        let mut config = test_config();
        config.generation.galaxies = 5;
        config.generation.star_systems = 2;
        let mut universe = Universe::new(config).unwrap();
        let secondary_galaxies: Vec<_> = universe
            .galaxies()
            .skip(1)
            .map(|galaxy| galaxy.id())
            .collect();
        for galaxy_id in secondary_galaxies {
            universe
                .materialize_region_in(galaxy_id, RegionCoord::ORIGIN)
                .unwrap();
        }

        let mut serial = universe.clone();
        let mut four_at_once = universe;
        serial
            .advance_ticks_bounded(10_000, NonZeroUsize::new(1).unwrap())
            .unwrap();
        four_at_once
            .advance_ticks_bounded(10_000, NonZeroUsize::new(4).unwrap())
            .unwrap();

        assert_eq!(serial, four_at_once);
    }

    #[test]
    fn all_galaxy_origins_are_explicit_and_survive_advancement_and_save() {
        let mut config = test_config();
        config.generation.galaxies = 3;
        config.generation.star_systems = 3;
        let mut universe = Universe::new(config).unwrap();

        assert_eq!(universe.materialize_all_galaxy_origins().unwrap(), 2);
        assert_eq!(universe.materialize_all_galaxy_origins().unwrap(), 0);
        assert_eq!(universe.materialized_counts().star_systems, 9);
        assert!(universe.materialized_counts().stars >= 9);
        assert!(universe.materialized_counts().planets > 0);
        assert!(universe.galaxy_regions().all(|(_, _, loaded)| loaded));

        universe.advance_ticks(10_000).unwrap();
        assert_eq!(universe.galaxy_regions().count(), 3);
        assert!(
            universe
                .galaxy_regions()
                .all(|(_, tick, loaded)| loaded && tick == 10_000)
        );
        assert_eq!(Universe::decode(&universe.encode()).unwrap(), universe);
    }

    #[test]
    fn aggregate_population_grows_by_materializing_sparse_regions() {
        let mut config = test_config();
        config.generation.galaxies = 2;
        config.generation.star_systems = 3;
        let mut universe = Universe::new(config).unwrap();
        let second = universe.galaxies().nth(1).unwrap().id();
        universe
            .materialize_region_in(second, RegionCoord::new(1, 0, 0))
            .unwrap();
        universe
            .materialize_region_in(second, RegionCoord::new(2, 0, 0))
            .unwrap();
        let total_systems: usize = universe
            .galaxy_regions()
            .filter_map(|(key, _, materialized)| {
                materialized
                    .then(|| {
                        universe
                            .region_objects_in(key.galaxy_id(), key.coordinate())
                            .map(|objects| objects.counts().star_systems)
                    })
                    .flatten()
            })
            .sum();
        assert_eq!(total_systems, 9);
        assert!(total_systems > config.generation.star_systems as usize);
    }

    #[test]
    fn dynamic_star_formation_is_eventful_and_batch_independent() {
        let mut config = test_config();
        config.generation.star_systems = 1;
        config.generation.star_mass_solar_min = 0.1;
        config.generation.star_mass_solar_max = 0.2;
        config.generation.planets_per_system_min = 0;
        config.generation.planets_per_system_max = 0;

        let mut uninterrupted = Universe::new(config).unwrap();
        uninterrupted.advance_ticks(20_000).unwrap();

        let mut resumed = Universe::new(config).unwrap();
        resumed.advance_ticks(10_000).unwrap();
        resumed = Universe::decode(&resumed.encode()).unwrap();
        resumed.advance_ticks(10_000).unwrap();

        assert_eq!(resumed, uninterrupted);
        assert!(uninterrupted.objects().counts().star_systems > 1);
        assert!(
            uninterrupted
                .events()
                .iter()
                .any(|event| event.kind() == crate::AstroEventKind::StarFormation)
        );
        assert!(uninterrupted.objects().stars().any(|star| {
            star.formed_at() > SimDuration::ZERO && star.epoch_age() == SimDuration::ZERO
        }));
    }

    #[test]
    fn interval_advance_returns_only_new_dynamic_events() {
        let mut config = test_config();
        config.generation.star_systems = 1;
        config.generation.planets_per_system_min = 1;
        config.generation.planets_per_system_max = 1;
        config.generation.orbit_au_min = 1.0;
        config.generation.orbit_au_max = 1.0;
        config.generation.star_mass_solar_min = 30.0;
        config.generation.star_mass_solar_max = 30.0;
        let state = Universe::new(config)
            .unwrap()
            .encode()
            .lines()
            .map(|line| {
                if !line.starts_with("S|") {
                    return line.to_owned();
                }
                let mut fields: Vec<_> = line.split('|').map(str::to_owned).collect();
                fields[4] = fields[3].clone();
                fields[5] = "0".to_owned();
                fields[6] = "0".to_owned();
                fields[10] = "main_sequence".to_owned();
                fields.join("|")
            })
            .collect::<Vec<_>>()
            .join("\n");
        let mut universe = Universe::decode(&state).unwrap();

        assert!(
            universe
                .advance_ticks_with_events(1_000)
                .unwrap()
                .is_empty()
        );
        let events = universe.advance_ticks_with_events(4_000).unwrap();
        assert!(!events.is_empty());
        assert!(events.iter().all(|event| event.tick() > 1_000));
        assert!(events.iter().any(|event| !event.impacts().is_empty()));
        assert!(
            universe
                .objects()
                .planets()
                .any(|planet| planet.status() == crate::objects::PlanetStatus::Destroyed)
        );
        assert_eq!(Universe::decode(&universe.encode()).unwrap(), universe);
    }

    #[test]
    fn legacy_version_twelve_reconstructs_birth_composition_without_replaying_feedback() {
        let mut universe = Universe::new(test_config()).unwrap();
        universe.advance_ticks(1_000).unwrap();
        let legacy = legacy_state_without_galaxy_sections(&universe)
            .lines()
            .filter(|line| !line.starts_with("H|") && !line.starts_with("L|"))
            .map(|line| {
                if line.starts_with("S|") {
                    line.split('|').take(19).collect::<Vec<_>>().join("|")
                } else if line.starts_with("Y|") {
                    line.split('|').take(20).collect::<Vec<_>>().join("|")
                } else {
                    line.replace("state_format_version=19", "state_format_version=12")
                        .replace("simulation_version=8", "simulation_version=1")
                        .replace("generation_version=15", "generation_version=9")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let mut restored = Universe::decode(&legacy).unwrap();
        assert_eq!(
            restored.objects().stars().count(),
            universe.objects().stars().count()
        );
        assert_eq!(
            restored.objects().gas_reservoir_solar_masses(),
            universe.objects().gas_reservoir_solar_masses()
        );
        restored.objects().validate().unwrap();
        restored.advance_ticks(0).unwrap();
        assert_eq!(
            restored.objects().gas_composition(),
            universe.objects().gas_composition()
        );
        assert_eq!(
            restored.objects().gas_reservoir_solar_masses(),
            universe.objects().gas_reservoir_solar_masses()
        );
    }

    #[test]
    fn version_thirteen_state_gains_cloud_and_cluster_ancestry() {
        let universe = Universe::new(test_config()).unwrap();
        let legacy = legacy_state_without_galaxy_sections(&universe)
            .lines()
            .filter(|line| !line.starts_with("H|") && !line.starts_with("L|"))
            .map(|line| {
                if line.starts_with("Y|") {
                    line.split('|').take(20).collect::<Vec<_>>().join("|")
                } else {
                    line.replace("state_format_version=19", "state_format_version=13")
                        .replace("simulation_version=8", "simulation_version=2")
                        .replace("generation_version=15", "generation_version=10")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let restored = Universe::decode(&legacy).unwrap();
        restored.objects().validate().unwrap();
        assert_eq!(
            restored.objects().gas_reservoir_solar_masses(),
            universe.objects().gas_reservoir_solar_masses()
        );
        assert_eq!(
            restored.objects().star_clusters().count(),
            restored.objects().star_systems().count()
        );
        assert_eq!(Universe::decode(&restored.encode()).unwrap(), restored);
    }

    #[test]
    fn version_fourteen_cloud_state_loads_without_shock_fields() {
        let universe = Universe::new(test_config()).unwrap();
        let legacy = legacy_state_without_galaxy_sections(&universe)
            .lines()
            .map(|line| {
                if line.starts_with("H|") {
                    line.split('|').take(14).collect::<Vec<_>>().join("|")
                } else if line.starts_with("L|") {
                    line.split('|').take(8).collect::<Vec<_>>().join("|")
                } else {
                    line.replace("state_format_version=19", "state_format_version=14")
                        .replace("simulation_version=8", "simulation_version=3")
                        .replace("generation_version=15", "generation_version=11")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let restored = Universe::decode(&legacy).unwrap();
        assert_eq!(restored, universe);
    }

    #[test]
    fn enriched_galaxies_match_across_checkpoints_and_lazy_materialization() {
        let mut config = test_config();
        config.generation.galaxies = 2;
        config.generation.star_systems = 16;
        config.generation.star_mass_solar_min = 8.0;
        config.generation.star_mass_solar_max = 30.0;
        config.generation.planets_per_system_min = 0;
        config.generation.planets_per_system_max = 0;
        let mut whole = Universe::new(config).unwrap();
        let galaxy_id = whole.galaxies().nth(1).unwrap().id();
        let coordinate = RegionCoord::new(2, 3, 4);
        whole.materialize_region_in(galaxy_id, coordinate).unwrap();
        let mut split = whole.clone();
        whole.advance_ticks(100_000).unwrap();
        for ticks in [3_333, 6_667, 8_123, 21_877, 60_000] {
            split.advance_ticks(ticks).unwrap();
            split = Universe::decode(&split.encode()).unwrap();
        }
        // Compare feedback and birth properties; orbital integration has its own fidelity rules.
        let feedback = |objects: &ObjectRegistry| {
            (
                objects.gas_reservoir_solar_masses(),
                objects.gas_composition(),
                objects
                    .stars()
                    .map(|star| {
                        (
                            star.id(),
                            star.current_mass(),
                            star.birth_composition(),
                            star.formed_at(),
                            star.stage(),
                        )
                    })
                    .collect::<Vec<_>>(),
            )
        };
        assert_eq!(feedback(split.objects()), feedback(whole.objects()));
        assert_eq!(
            feedback(split.region_objects_in(galaxy_id, coordinate).unwrap()),
            feedback(whole.region_objects_in(galaxy_id, coordinate).unwrap())
        );
        let expected = feedback(whole.region_objects_in(galaxy_id, coordinate).unwrap());
        whole.unload_region_in(galaxy_id, coordinate).unwrap();
        whole.materialize_region_in(galaxy_id, coordinate).unwrap();
        assert_eq!(
            feedback(whole.region_objects_in(galaxy_id, coordinate).unwrap()),
            expected
        );
        assert!(whole.objects().gas_composition().metallicity() > 0.014);
    }

    #[test]
    fn state_round_trips_without_changing_history() {
        let mut universe = Universe::new(test_config()).unwrap();
        universe.advance_ticks(123).unwrap();
        assert_eq!(Universe::decode(&universe.encode()).unwrap(), universe);
    }

    #[test]
    fn save_load_continue_matches_uninterrupted_simulation() {
        let mut uninterrupted = Universe::new(test_config()).unwrap();
        uninterrupted.advance_ticks(1_000).unwrap();

        let mut resumed = Universe::new(test_config()).unwrap();
        resumed.advance_ticks(400).unwrap();
        let mut resumed = Universe::decode(&resumed.encode()).unwrap();
        resumed.advance_ticks(600).unwrap();
        assert_eq!(resumed, uninterrupted);
    }

    #[test]
    fn state_rejects_an_inconsistent_clock() {
        let state = Universe::new(test_config())
            .unwrap()
            .encode()
            .replace("elapsed_seconds=0", "elapsed_seconds=1");
        assert!(matches!(
            Universe::decode(&state),
            Err(StateError::Time(TimeError::InconsistentClock))
        ));
    }

    #[test]
    fn state_rejects_unsupported_versions() {
        let state = Universe::new(test_config())
            .unwrap()
            .encode()
            .replace("rng_version=1", "rng_version=99");
        assert!(matches!(
            Universe::decode(&state),
            Err(StateError::Compatibility(_))
        ));
    }

    #[test]
    fn state_rejects_broken_relationships() {
        let universe = Universe::new(test_config()).unwrap();
        let first_planet = universe.objects().planets().next().unwrap().id().value();
        let state = universe
            .encode()
            .lines()
            .filter(|line| !line.starts_with(&format!("P|{first_planet}|")))
            .collect::<Vec<_>>()
            .join("\n");
        assert!(matches!(
            Universe::decode(&state),
            Err(StateError::Objects(_))
        ));
    }

    #[test]
    fn massive_stars_produce_ordered_lifecycle_events() {
        let mut config = test_config();
        config.generation.star_systems = 1;
        config.generation.planets_per_system_min = 0;
        config.generation.planets_per_system_max = 0;
        config.generation.star_mass_solar_min = 30.0;
        config.generation.star_mass_solar_max = 30.0;
        let state = Universe::new(config)
            .unwrap()
            .encode()
            .lines()
            .map(|line| {
                if !line.starts_with("S|") {
                    return line.to_owned();
                }
                let mut fields: Vec<_> = line.split('|').map(str::to_owned).collect();
                fields[4] = fields[3].clone();
                fields[5] = "0".to_owned();
                fields[6] = "0".to_owned();
                fields[10] = "main_sequence".to_owned();
                fields.join("|")
            })
            .collect::<Vec<_>>()
            .join("\n");
        let mut universe = Universe::decode(&state).unwrap();
        let primary_id = universe.objects().star_systems().next().unwrap().star_id();
        universe.advance_ticks(5_000).unwrap();

        let primary_events: Vec<_> = universe
            .events()
            .iter()
            .filter(|event| event.star_id() == primary_id)
            .collect();
        assert!(primary_events.len() >= 3);
        assert!(
            primary_events
                .windows(2)
                .all(|events| events[0].tick() <= events[1].tick())
        );
        assert!(
            primary_events
                .iter()
                .any(|event| event.kind() == crate::AstroEventKind::CoreCollapseSupernova)
        );
        assert_eq!(
            primary_events.last().unwrap().kind(),
            crate::AstroEventKind::BlackHoleFormation
        );
        assert_eq!(
            universe.objects().stars().next().unwrap().stage(),
            crate::objects::StellarStage::BlackHole
        );
    }

    #[test]
    fn unloaded_region_rematerializes_to_the_same_current_state() {
        let coordinate = RegionCoord::new(18, -4, 91);
        let mut universe = Universe::new(test_config()).unwrap();
        universe.advance_ticks(20_000).unwrap();
        assert!(universe.materialize_region(coordinate).unwrap());
        let expected = universe.region_objects(coordinate).unwrap().clone();
        assert!(universe.unload_region(coordinate).unwrap());
        assert!(universe.region_objects(coordinate).is_none());

        let mut universe = Universe::decode(&universe.encode()).unwrap();
        assert!(universe.materialize_region(coordinate).unwrap());
        assert_eq!(universe.region_objects(coordinate).unwrap(), &expected);
    }

    #[test]
    fn version_two_state_is_migrated_into_the_origin_region() {
        let encoded = Universe::new(test_config()).unwrap().encode();
        let companion_ids: BTreeSet<_> = encoded
            .lines()
            .filter(|line| line.starts_with("Y|"))
            .flat_map(|line| {
                line.split('|')
                    .nth(12)
                    .unwrap_or_default()
                    .split(',')
                    .filter(|value| !value.is_empty())
            })
            .collect();
        let mut migrated = Vec::new();
        for line in encoded.lines() {
            if line == "[regions]" {
                break;
            }
            if line.starts_with("N|")
                || line.starts_with("A|")
                || line.starts_with("C|")
                || line.starts_with("D|")
                || line.starts_with("H|")
                || line.starts_with("L|")
            {
                continue;
            }
            if line.starts_with("Y|") {
                migrated.push(line.split('|').take(11).collect::<Vec<_>>().join("|"));
            } else if line.starts_with("S|") {
                let fields: Vec<_> = line.split('|').collect();
                if companion_ids.contains(fields[1]) {
                    continue;
                }
                migrated.push(
                    [
                        fields[0], fields[1], fields[2], fields[3], "0", fields[7], fields[8],
                        fields[9], fields[10],
                    ]
                    .join("|"),
                );
            } else if line.starts_with("P|") {
                migrated.push(line.split('|').take(12).collect::<Vec<_>>().join("|"));
            } else {
                migrated.push(
                    line.replace("state_format_version=19", "state_format_version=2")
                        .replace("generation_version=15", "generation_version=1")
                        .replace("simulation_version=8", "simulation_version=1"),
                );
            }
        }
        let decoded = Universe::decode(&migrated.join("\n")).unwrap();
        assert_eq!(decoded.known_regions().count(), 1);
        assert!(decoded.region_objects(RegionCoord::ORIGIN).is_some());
    }

    #[test]
    fn version_six_state_is_accepted_and_upgraded() {
        let encoded = legacy_state_without_galaxy_sections(&Universe::new(test_config()).unwrap());
        let encoded = encoded
            .lines()
            .filter(|line| {
                !line.starts_with("N|")
                    && !line.starts_with("A|")
                    && !line.starts_with("C|")
                    && !line.starts_with("D|")
                    && !line.starts_with("H|")
                    && !line.starts_with("L|")
            })
            .filter(|line| !line.starts_with("galaxies="))
            .map(|line| {
                if line.starts_with("M|") {
                    line.split('|').take(3).collect::<Vec<_>>().join("|")
                } else if line.starts_with("Y|") {
                    line.split('|').take(13).collect::<Vec<_>>().join("|")
                } else if line.starts_with("R|") || line.starts_with("E|") {
                    let mut fields: Vec<_> = line.split('|').collect();
                    fields.remove(if line.starts_with("R|") { 1 } else { 3 });
                    fields.join("|")
                } else {
                    line.replace("state_format_version=19", "state_format_version=6")
                        .replace("generation_version=15", "generation_version=3")
                        .replace("simulation_version=8", "simulation_version=1")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let decoded = Universe::decode(&encoded).unwrap();
        let upgraded = decoded.encode();
        assert!(upgraded.contains("state_format_version=19"));
        assert!(upgraded.contains("generation_version=15"));
        assert_eq!(decoded.objects().gas_composition(), Default::default());
    }

    #[test]
    fn clock_only_version_one_state_is_migrated() {
        let config = test_config();
        let universe = Universe::new(config).unwrap();
        let elapsed = config.tick_duration.duration().checked_mul(10_000).unwrap();
        let state = format!(
            "COSMOS_STATE_V1\nversion=0.1.0\nseed={}\nuniverse_id={}\ntick=10000\nelapsed_seconds={}\ntick_duration_seconds={}\n",
            universe.seed(),
            universe.id().value(),
            elapsed.seconds(),
            config.tick_duration.duration().seconds(),
        );
        let decoded = Universe::decode(&state).unwrap();
        assert_eq!(decoded.clock().tick(), 10_000);
        assert_eq!(decoded.objects().counts().star_systems, 4);
        assert!(
            decoded
                .objects()
                .stars()
                .all(|star| { star.epoch_age().checked_add(elapsed) == Some(star.age()) })
        );
    }
}
