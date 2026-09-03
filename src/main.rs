use std::error::Error;
use std::fmt;
use std::io::{self, Write as _};
use std::num::{NonZeroU64, NonZeroUsize};
use std::path::PathBuf;
use std::process::ExitCode;
use std::thread;
use std::time::Duration;

use clap::{ArgAction, Args, Parser, Subcommand};
use universe_sim::{
    GalaxyId, GasComposition, MAX_CONCURRENT_GALAXIES, ObjectId, ObjectRef, RegionCoord,
    SimDuration, SimulationConfig, TickDuration, Universe,
};

#[derive(Debug, Parser)]
#[command(
    name = "universe-sim",
    version,
    about = "Deterministic astrophysical simulation foundation",
    after_help = "Durations use integer values such as 1s, '1 hour', 1000y, '10 kyr', or '1 Myr'."
)]
struct Cli {
    /// Emit structured diagnostic records to standard error.
    #[arg(short, long, action = ArgAction::Count, global = true)]
    verbose: u8,

    /// Worker threads to use; defaults to the CPUs available to this process.
    #[arg(long, value_name = "COUNT", global = true)]
    threads: Option<NonZeroUsize>,

    #[command(subcommand)]
    command: CliCommand,
}

#[derive(Debug, Subcommand)]
enum CliCommand {
    /// Create and save a new seeded universe.
    New(NewArgs),
    /// Advance a saved universe by a number of universal ticks.
    Simulate(SimulateArgs),
    /// Print the metadata in a saved universe.
    Inspect(InspectArgs),
    /// Generate a procedural region and keep it active in the saved state.
    Materialize(MaterializeArgs),
    /// Unload a non-origin region while retaining its sparse history.
    Unload(UnloadArgs),
    /// List known regions and their materialization state.
    Regions(RegionsArgs),
    /// List the deterministic galaxy catalog.
    Galaxies(GalaxiesArgs),
    /// List recorded stellar lifecycle events.
    Events(EventsArgs),
}

#[derive(Debug, Args)]
struct NewArgs {
    /// Seed used to generate the universe.
    #[arg(long, required_unless_present = "config", conflicts_with = "config")]
    seed: Option<u64>,

    /// Duration represented by one universal tick; defaults to one year.
    #[arg(
        long,
        value_name = "DURATION",
        value_parser = parse_tick_duration,
        conflicts_with = "config"
    )]
    tick_length: Option<TickDuration>,

    /// Number of star systems to generate initially; stars and planets vary procedurally.
    #[arg(long, value_name = "COUNT")]
    star_systems: Option<u32>,

    /// Number of deterministic galaxies in the sparse universe catalog.
    #[arg(long, value_name = "COUNT")]
    galaxies: Option<u16>,

    /// Read the seed and tick duration from a configuration file.
    #[arg(long, value_name = "PATH", conflicts_with_all = ["seed", "tick_length"])]
    config: Option<PathBuf>,

    /// Destination for the initial universe state.
    #[arg(long, value_name = "PATH", default_value = "universe.state")]
    output: PathBuf,
}

#[derive(Debug, Args)]
struct SimulateArgs {
    /// Existing universe state to advance.
    #[arg(long, value_name = "PATH")]
    state: PathBuf,

    /// Number of universal ticks to advance.
    #[arg(
        long,
        required_unless_present = "continuous",
        conflicts_with = "continuous"
    )]
    ticks: Option<u64>,

    /// Continue advancing and checkpointing until interrupted.
    #[arg(long, conflicts_with = "ticks")]
    continuous: bool,

    /// Save after each batch of this many ticks.
    #[arg(long)]
    checkpoint_every: Option<NonZeroU64>,

    /// Pause after each continuous checkpoint; defaults to 1000 milliseconds.
    #[arg(long, requires = "continuous")]
    pace_ms: Option<u64>,

    /// Destination state; defaults to updating --state in place.
    #[arg(long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Maximum galaxies processed concurrently; all materialized galaxies still advance.
    #[arg(
        long,
        value_name = "COUNT",
        default_value = "4",
        value_parser = parse_galaxies_at_once
    )]
    galaxies_at_once: NonZeroUsize,
}

#[derive(Debug, Args)]
struct InspectArgs {
    /// Universe state to inspect.
    #[arg(long, value_name = "PATH")]
    state: PathBuf,

    /// Galaxy containing the region; defaults to the primary galaxy.
    #[arg(long, value_name = "G-ID")]
    galaxy: Option<GalaxyId>,

    /// Region containing the object.
    #[arg(long, default_value = "0,0,0", allow_hyphen_values = true)]
    region: RegionCoord,

    /// Optional object to inspect (G-, SYS-, S-, or P-).
    object: Option<ObjectId>,
}

#[derive(Debug, Args)]
struct MaterializeArgs {
    #[arg(long, value_name = "PATH")]
    state: PathBuf,

    /// Galaxy to materialize within; defaults to the primary galaxy.
    #[arg(long, value_name = "G-ID")]
    galaxy: Option<GalaxyId>,

    /// Procedural region coordinate in x,y,z form.
    #[arg(long, allow_hyphen_values = true)]
    region: RegionCoord,

    #[arg(long, value_name = "PATH")]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct UnloadArgs {
    #[arg(long, value_name = "PATH")]
    state: PathBuf,

    /// Galaxy containing the region; defaults to the primary galaxy.
    #[arg(long, value_name = "G-ID")]
    galaxy: Option<GalaxyId>,

    /// Region coordinate in x,y,z form.
    #[arg(long, allow_hyphen_values = true)]
    region: RegionCoord,

    #[arg(long, value_name = "PATH")]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct RegionsArgs {
    #[arg(long, value_name = "PATH")]
    state: PathBuf,

    /// Show only regions in this galaxy.
    #[arg(long, value_name = "G-ID")]
    galaxy: Option<GalaxyId>,
}

#[derive(Debug, Args)]
struct GalaxiesArgs {
    #[arg(long, value_name = "PATH")]
    state: PathBuf,
}

#[derive(Debug, Args)]
struct EventsArgs {
    #[arg(long, value_name = "PATH")]
    state: PathBuf,

    /// Maximum number of most-recent events to display.
    #[arg(long, default_value_t = 50)]
    limit: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    if let Some(threads) = cli.threads
        && let Err(error) = rayon::ThreadPoolBuilder::new()
            .num_threads(threads.get())
            .thread_name(|index| format!("universe-worker-{index}"))
            .build_global()
    {
        eprintln!("error: could not configure worker threads: {error}");
        return ExitCode::FAILURE;
    }
    let diagnostics = Diagnostics::new(cli.verbose > 0);
    match run(cli.command, &diagnostics) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(command: CliCommand, diagnostics: &Diagnostics) -> Result<(), CliError> {
    match command {
        CliCommand::New(args) => new_universe(args, diagnostics),
        CliCommand::Simulate(args) => simulate(args, diagnostics),
        CliCommand::Inspect(args) => inspect(args, diagnostics),
        CliCommand::Materialize(args) => materialize(args, diagnostics),
        CliCommand::Unload(args) => unload(args, diagnostics),
        CliCommand::Regions(args) => list_regions(args, diagnostics),
        CliCommand::Galaxies(args) => list_galaxies(args, diagnostics),
        CliCommand::Events(args) => list_events(args, diagnostics),
    }
}

fn new_universe(args: NewArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let mut config = if let Some(path) = args.config {
        SimulationConfig::from_file(path)?
    } else {
        SimulationConfig::new(
            args.seed
                .expect("clap requires --seed when --config is absent"),
            args.tick_length.unwrap_or(default_tick_duration()),
        )
    };
    if let Some(star_systems) = args.star_systems {
        config.generation.star_systems = star_systems;
    }
    if let Some(galaxies) = args.galaxies {
        config.generation.galaxies = galaxies;
    }

    let universe = Universe::new(config)?;
    diagnostics.emit("universe_generated", &universe, None);
    universe.save(&args.output)?;
    diagnostics.emit("state_saved", &universe, Some(&args.output));
    let counts = universe.objects().counts();
    println!("Universe created.");
    println!("  ID: {}", universe.id());
    println!("  Seed: {}", universe.seed());
    println!("  Tick length: {}", universe.clock().tick_duration());
    println!("  Galaxies: {}", universe.galaxies().count());
    println!("  Star systems: {}", counts.star_systems);
    println!("  Stars: {}", counts.stars);
    println!("  Planets: {}", counts.planets);
    println!(
        "  Star-forming gas: {:.3} solar masses",
        universe.objects().gas_reservoir_solar_masses()
    );
    print_gas_composition(universe.objects().gas_composition());
    println!("  Total bodies: {}", counts.stars + counts.planets);
    println!(
        "  Total objects: {}",
        universe.galaxies().count()
            + counts.star_systems
            + counts.stars
            + counts.planets
            + counts.gas_reservoirs
    );
    println!("  State: {}", args.output.display());
    Ok(())
}

fn simulate(args: SimulateArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let output = args.output.unwrap_or_else(|| args.state.clone());
    let mut universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    if args.continuous {
        let batch_size = args.checkpoint_every.map(NonZeroU64::get).unwrap_or(1_000);
        let pace = Duration::from_millis(args.pace_ms.unwrap_or(1_000));
        if diagnostics.enabled {
            println!("Continuous simulation started.");
            println!("  Checkpoint every: {batch_size} ticks");
            println!("  Pace: {} ms per checkpoint", pace.as_millis());
            println!("  Stop with Ctrl+C; the last completed checkpoint remains valid.");
            io::stdout().flush().map_err(CliError::from)?;
        }
        loop {
            let new_events = universe
                .advance_ticks_with_events_bounded(batch_size, args.galaxies_at_once)?
                .to_vec();
            diagnostics.emit("simulation_advanced", &universe, None);
            universe.save(&output)?;
            diagnostics.emit("checkpoint_saved", &universe, Some(&output));
            if diagnostics.enabled {
                println!(
                    "Checkpoint: tick {}, elapsed {}, new events {}",
                    universe.clock().tick(),
                    universe.clock().elapsed(),
                    new_events.len()
                );
            }
            for event in &new_events {
                print_event(event, "  ", diagnostics.enabled);
            }
            io::stdout().flush().map_err(CliError::from)?;
            if !pace.is_zero() {
                thread::sleep(pace);
            }
        }
    }

    let ticks = args
        .ticks
        .expect("clap requires --ticks when --continuous is absent");
    let starting_tick = universe.clock().tick();
    let batch_size = args
        .checkpoint_every
        .map(NonZeroU64::get)
        .unwrap_or(ticks.max(1));
    let mut remaining = ticks;
    while remaining > 0 {
        let batch = remaining.min(batch_size);
        universe.advance_ticks_bounded(batch, args.galaxies_at_once)?;
        remaining -= batch;
        diagnostics.emit("simulation_advanced", &universe, None);
        universe.save(&output)?;
        diagnostics.emit("checkpoint_saved", &universe, Some(&output));
    }
    if ticks == 0 {
        universe.save(&output)?;
        diagnostics.emit("checkpoint_saved", &universe, Some(&output));
    }

    let new_events = events_since(&universe, starting_tick);
    if diagnostics.enabled {
        println!("Simulation complete.");
        println!("  Universal tick: {}", universe.clock().tick());
        println!("  Elapsed time: {}", universe.clock().elapsed());
        println!("  New events: {}", new_events.len());
    }
    for event in &new_events {
        print_event(
            event,
            if diagnostics.enabled { "    " } else { "" },
            diagnostics.enabled,
        );
    }
    if diagnostics.enabled {
        println!("  State: {}", output.display());
    }
    Ok(())
}

fn events_since(universe: &Universe, starting_tick: u64) -> Vec<&universe_sim::AstroEvent> {
    universe
        .events()
        .iter()
        .filter(|event| event.tick() > starting_tick)
        .collect()
}

fn print_event(event: &universe_sim::AstroEvent, indentation: &str, detailed: bool) {
    if event.from() == event.to() {
        println!(
            "{indentation}Tick {}: {} ({}) occurred in {} region {}",
            event.tick(),
            event.star_id(),
            event.kind(),
            event.galaxy_id(),
            event.region()
        );
    } else {
        println!(
            "{indentation}Tick {}: {} ({}) changed from {} to {} in {} region {}",
            event.tick(),
            event.star_id(),
            event.kind(),
            event.from(),
            event.to(),
            event.galaxy_id(),
            event.region()
        );
    }
    if detailed {
        if event.natal_kick_meters_per_second() > 0 {
            println!(
                "{indentation}  Natal kick: {:.1} km/s",
                f64::from(event.natal_kick_meters_per_second()) / 1_000.0
            );
        }
        for impact in event.impacts() {
            println!(
                "{indentation}  Impact: {} was {}",
                impact.planet_id, impact.kind
            );
        }
        if let Some(related_star_id) = event.related_star_id() {
            println!("{indentation}  Related star: {related_star_id}");
        }
    } else if !event.impacts().is_empty() {
        println!(
            "{indentation}  Planetary impacts: {} (use --verbose for details)",
            event.impacts().len()
        );
    }
}

fn inspect(args: InspectArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    let galaxy_id = args.galaxy.unwrap_or_else(|| universe.primary_galaxy_id());
    if let Some(id) = args.object {
        let objects = universe
            .region_objects_in(galaxy_id, args.region)
            .ok_or_else(|| {
                CliError::message(format!(
                    "region {} in {} is not materialized; run the materialize command first",
                    args.region, galaxy_id
                ))
            })?;
        let dynamics = match id {
            ObjectId::StarSystem(system_id) => objects.dynamics_mode(system_id),
            _ => None,
        };
        let object = objects
            .find(id)
            .ok_or_else(|| CliError::message(format!("object {id} was not found")))?;
        print_object(object);
        if let Some(dynamics) = dynamics {
            println!("  Dynamics: {dynamics}");
        }
        return Ok(());
    }

    let counts = universe.objects().counts();
    println!("UNIVERSE {}", universe.id());
    println!("  Seed: {}", universe.seed());
    println!("  Universal tick: {}", universe.clock().tick());
    println!("  Tick length: {}", universe.clock().tick_duration());
    println!("  Elapsed time: {}", universe.clock().elapsed());
    println!("  Galaxies: {}", universe.galaxies().count());
    println!("  Star systems: {}", counts.star_systems);
    println!("  Stars: {}", counts.stars);
    println!("  Planets: {}", counts.planets);
    println!(
        "  Star-forming gas: {:.3} solar masses",
        universe.objects().gas_reservoir_solar_masses()
    );
    print_gas_composition(universe.objects().gas_composition());
    println!("  Recorded events: {}", universe.events().len());
    println!("  Known regions: {}", universe.galaxy_regions().count());
    Ok(())
}

fn materialize(args: MaterializeArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let output = args.output.unwrap_or_else(|| args.state.clone());
    let mut universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    let galaxy_id = args.galaxy.unwrap_or_else(|| universe.primary_galaxy_id());
    let created = universe.materialize_region_in(galaxy_id, args.region)?;
    universe.save(&output)?;
    diagnostics.emit("state_saved", &universe, Some(&output));
    let objects = universe
        .region_objects_in(galaxy_id, args.region)
        .expect("materialized region is available");
    let counts = objects.counts();
    println!(
        "Region {} {} in {}.",
        args.region,
        if created {
            "materialized"
        } else {
            "was already materialized"
        },
        galaxy_id
    );
    println!("  Star systems: {}", counts.star_systems);
    println!("  Stars: {}", counts.stars);
    println!("  Planets: {}", counts.planets);
    println!(
        "  Star-forming gas: {:.3} solar masses",
        objects.gas_reservoir_solar_masses()
    );
    print_gas_composition(objects.gas_composition());
    if let Some(star) = objects.stars().next() {
        println!("  Example star: {}", star.id());
    }
    println!("  State: {}", output.display());
    Ok(())
}

fn unload(args: UnloadArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let output = args.output.unwrap_or_else(|| args.state.clone());
    let mut universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    let galaxy_id = args.galaxy.unwrap_or_else(|| universe.primary_galaxy_id());
    let removed = universe.unload_region_in(galaxy_id, args.region)?;
    universe.save(&output)?;
    diagnostics.emit("state_saved", &universe, Some(&output));
    println!(
        "Region {} {} in {}.",
        args.region,
        if removed {
            "unloaded"
        } else {
            "was already unloaded"
        },
        galaxy_id
    );
    println!("  Sparse history retained: yes");
    println!("  State: {}", output.display());
    Ok(())
}

fn list_regions(args: RegionsArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    println!("KNOWN REGIONS");
    for (key, last_updated_tick, materialized) in universe.galaxy_regions() {
        if args.galaxy.is_some_and(|galaxy| galaxy != key.galaxy_id()) {
            continue;
        }
        println!(
            "  {}: {}, last updated at tick {}",
            key,
            if materialized {
                "materialized"
            } else {
                "unloaded"
            },
            last_updated_tick
        );
    }
    Ok(())
}

fn list_galaxies(args: GalaxiesArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    println!("GALAXIES");
    for galaxy in universe.galaxies() {
        let [x, y, z] = galaxy.position().meters();
        println!(
            "  {}: {}, radius {:.0} pc, position [{:.3e}, {:.3e}, {:.3e}] m",
            galaxy.id(),
            galaxy.morphology(),
            galaxy.radius_parsecs(),
            x,
            y,
            z
        );
    }
    Ok(())
}

fn list_events(args: EventsArgs, diagnostics: &Diagnostics) -> Result<(), CliError> {
    let universe = Universe::load(&args.state)?;
    diagnostics.emit("state_loaded", &universe, Some(&args.state));
    println!("STELLAR EVENTS ({})", universe.events().len());
    let skip = universe.events().len().saturating_sub(args.limit);
    for event in universe.events().iter().skip(skip) {
        print_event(event, "  ", true);
    }
    Ok(())
}

fn print_object(object: ObjectRef<'_>) {
    match object {
        ObjectRef::Galaxy(galaxy) => {
            let [x, y, z] = galaxy.position().meters();
            println!("GALAXY {}", galaxy.id());
            println!("  Position: [{x:.3e}, {y:.3e}, {z:.3e}] m");
            println!("  Star systems: {}", galaxy.star_system_ids().len());
        }
        ObjectRef::StarSystem(system) => {
            let [x, y, z] = system.position().meters();
            let [vx, vy, vz] = system.velocity().meters_per_second();
            println!("STAR SYSTEM {}", system.id());
            println!("  Galaxy: {}", system.galaxy_id());
            println!("  Galactic component: {}", system.component());
            println!("  Star: {}", system.star_id());
            if !system.companion_star_ids().is_empty() {
                println!(
                    "  Companion stars: {}",
                    system
                        .companion_star_ids()
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            }
            println!("  Planets: {}", system.planet_ids().len());
            println!("  Position: [{x:.3e}, {y:.3e}, {z:.3e}] m");
            println!("  Velocity: [{vx:.3e}, {vy:.3e}, {vz:.3e}] m/s");
        }
        ObjectRef::Star(star) => {
            let [x, y, z] = star.position().meters();
            let [vx, vy, vz] = star.velocity().meters_per_second();
            println!("STAR {}", star.id());
            println!("  System: {}", star.system_id());
            println!(
                "  Initial mass: {:.6} solar masses",
                star.initial_mass().solar_masses()
            );
            println!(
                "  Current mass: {:.6} solar masses",
                star.current_mass().solar_masses()
            );
            println!("  Age at epoch: {}", star.epoch_age());
            println!("  Age: {}", star.age());
            println!("  Stage: {}", star.stage());
            println!("  Metallicity: {:.6}", star.metallicity());
            println!(
                "  Luminosity: {:.6} solar luminosities",
                star.luminosity().solar_luminosities()
            );
            println!("  Temperature: {:.2} K", star.temperature().kelvin());
            println!("  Position: [{x:.3e}, {y:.3e}, {z:.3e}] m");
            println!("  Velocity: [{vx:.3e}, {vy:.3e}, {vz:.3e}] m/s");
            if let Some(orbit) = star.orbital_elements() {
                println!(
                    "  Companion orbit: {:.6} AU, e={:.4}",
                    orbit.semi_major_axis().astronomical_units(),
                    orbit.eccentricity()
                );
            }
        }
        ObjectRef::Planet(planet) => {
            let [x, y, z] = planet.position().meters();
            let [vx, vy, vz] = planet.velocity().meters_per_second();
            println!("PLANET {}", planet.id());
            println!("  System: {}", planet.system_id());
            println!("  Star: {}", planet.star_id());
            println!("  Mass: {:.6} Earth masses", planet.mass().earth_masses());
            println!("  Status: {}", planet.status());
            println!(
                "  Orbital radius: {:.6} AU",
                planet.orbital_radius().astronomical_units()
            );
            println!(
                "  Eccentricity: {:.6}",
                planet.orbital_elements().eccentricity()
            );
            println!(
                "  Inclination: {:.3} degrees",
                planet
                    .orbital_elements()
                    .inclination()
                    .radians()
                    .to_degrees()
            );
            println!("  Position: [{x:.3e}, {y:.3e}, {z:.3e}] m");
            println!("  Velocity: [{vx:.3e}, {vy:.3e}, {vz:.3e}] m/s");
        }
    }
}

fn parse_tick_duration(value: &str) -> Result<TickDuration, String> {
    value
        .parse::<TickDuration>()
        .map_err(|error| error.to_string())
}

fn print_gas_composition(composition: GasComposition) {
    println!(
        "  Gas composition: H {:.3}%, He {:.3}%, C {:.4}%, O {:.4}%, Fe {:.4}%, other metals {:.4}%",
        composition.hydrogen_fraction() * 100.0,
        composition.helium_fraction() * 100.0,
        composition.carbon_fraction() * 100.0,
        composition.oxygen_fraction() * 100.0,
        composition.iron_fraction() * 100.0,
        composition.other_metals_fraction() * 100.0,
    );
}

fn parse_galaxies_at_once(value: &str) -> Result<NonZeroUsize, String> {
    let value = value
        .parse::<usize>()
        .map_err(|_| "galaxies-at-once must be an integer from 1 through 4".to_owned())?;
    if !(1..=MAX_CONCURRENT_GALAXIES).contains(&value) {
        return Err(format!(
            "galaxies-at-once must be between 1 and {MAX_CONCURRENT_GALAXIES}"
        ));
    }
    Ok(NonZeroUsize::new(value).expect("validated as non-zero"))
}

fn default_tick_duration() -> TickDuration {
    TickDuration::new(SimDuration::from_years(1).expect("one year is representable"))
        .expect("one year is nonzero")
}

#[derive(Debug)]
struct CliError {
    message: String,
    source: Option<Box<dyn Error>>,
}

impl CliError {
    fn message(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            source: None,
        }
    }
}

impl fmt::Display for CliError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl Error for CliError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.source.as_deref()
    }
}

macro_rules! cli_error_from {
    ($error:ty) => {
        impl From<$error> for CliError {
            fn from(error: $error) -> Self {
                Self {
                    message: error.to_string(),
                    source: Some(Box::new(error)),
                }
            }
        }
    };
}

cli_error_from!(universe_sim::config::ConfigLoadError);
cli_error_from!(universe_sim::simulation::StateError);
cli_error_from!(universe_sim::simulation::SimulationAdvanceError);
cli_error_from!(universe_sim::time::TimeError);
cli_error_from!(universe_sim::simulation::UniverseCreateError);
cli_error_from!(universe_sim::simulation::RegionError);
cli_error_from!(std::io::Error);

#[derive(Debug, Clone, Copy)]
struct Diagnostics {
    enabled: bool,
}

impl Diagnostics {
    const fn new(enabled: bool) -> Self {
        Self { enabled }
    }

    fn emit(&self, event: &str, universe: &Universe, path: Option<&std::path::Path>) {
        if !self.enabled {
            return;
        }
        eprint!(
            "level=info event={event} universe_id={} tick={} threads={}",
            universe.id(),
            universe.clock().tick(),
            rayon::current_num_threads(),
        );
        if let Some(path) = path {
            eprint!(" path={}", path.display());
        }
        eprintln!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_seeded_new_command() {
        let cli = Cli::try_parse_from([
            "universe-sim",
            "new",
            "--seed",
            "42",
            "--tick-length",
            "10 kyr",
        ])
        .unwrap();

        let CliCommand::New(args) = cli.command else {
            panic!("expected new command");
        };
        assert_eq!(args.seed, Some(42));
        assert_eq!(args.tick_length.unwrap().to_string(), "10 kyr");
        assert_eq!(args.star_systems, None);
        assert_eq!(args.galaxies, None);
        assert_eq!(args.output, PathBuf::from("universe.state"));
    }

    #[test]
    fn parses_a_global_worker_count() {
        let cli =
            Cli::try_parse_from(["universe-sim", "--threads", "4", "new", "--seed", "42"]).unwrap();
        assert_eq!(cli.threads.map(NonZeroUsize::get), Some(4));
    }

    #[test]
    fn accepts_configuration_as_the_generation_source() {
        let cli = Cli::try_parse_from([
            "universe-sim",
            "new",
            "--config",
            "cosmos.toml",
            "--output",
            "test.state",
        ])
        .unwrap();

        let CliCommand::New(args) = cli.command else {
            panic!("expected new command");
        };
        assert_eq!(args.config, Some(PathBuf::from("cosmos.toml")));
        assert_eq!(args.output, PathBuf::from("test.state"));
    }

    #[test]
    fn parses_an_initial_star_system_count() {
        let cli = Cli::try_parse_from([
            "universe-sim",
            "new",
            "--seed",
            "42",
            "--star-systems",
            "250",
        ])
        .unwrap();

        let CliCommand::New(args) = cli.command else {
            panic!("expected new command");
        };
        assert_eq!(args.star_systems, Some(250));
    }

    #[test]
    fn parses_a_sparse_galaxy_count() {
        let cli = Cli::try_parse_from(["universe-sim", "new", "--seed", "42", "--galaxies", "12"])
            .unwrap();
        let CliCommand::New(args) = cli.command else {
            panic!("expected new command");
        };
        assert_eq!(args.galaxies, Some(12));
    }

    #[test]
    fn rejects_conflicting_generation_sources() {
        let result = Cli::try_parse_from([
            "universe-sim",
            "new",
            "--seed",
            "42",
            "--config",
            "cosmos.toml",
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn requires_a_generation_source() {
        assert!(Cli::try_parse_from(["universe-sim", "new"]).is_err());
    }

    #[test]
    fn parses_negative_region_coordinates() {
        let cli = Cli::try_parse_from([
            "universe-sim",
            "materialize",
            "--state",
            "universe.state",
            "--region",
            "-18,4,-91",
        ])
        .unwrap();
        let CliCommand::Materialize(args) = cli.command else {
            panic!("expected materialize command");
        };
        assert_eq!(args.region, RegionCoord::new(-18, 4, -91));
    }

    #[test]
    fn parses_continuous_simulation_mode() {
        let cli = Cli::try_parse_from([
            "universe-sim",
            "simulate",
            "--state",
            "universe.state",
            "--continuous",
            "--checkpoint-every",
            "500",
            "--pace-ms",
            "10",
        ])
        .unwrap();
        let CliCommand::Simulate(args) = cli.command else {
            panic!("expected simulate command");
        };
        assert!(args.continuous);
        assert_eq!(args.ticks, None);
        assert_eq!(args.checkpoint_every.unwrap().get(), 500);
        assert_eq!(args.pace_ms, Some(10));
        assert_eq!(args.galaxies_at_once.get(), 4);
    }

    #[test]
    fn parses_bounded_galaxy_parallelism() {
        let cli = Cli::try_parse_from([
            "universe-sim",
            "simulate",
            "--state",
            "universe.state",
            "--ticks",
            "10",
            "--galaxies-at-once",
            "2",
        ])
        .unwrap();
        let CliCommand::Simulate(args) = cli.command else {
            panic!("expected simulate command");
        };
        assert_eq!(args.galaxies_at_once.get(), 2);
        assert!(
            Cli::try_parse_from([
                "universe-sim",
                "simulate",
                "--state",
                "universe.state",
                "--ticks",
                "10",
                "--galaxies-at-once",
                "5",
            ])
            .is_err()
        );
    }

    #[test]
    fn rejects_ticks_with_continuous_mode() {
        assert!(
            Cli::try_parse_from([
                "universe-sim",
                "simulate",
                "--state",
                "universe.state",
                "--ticks",
                "10",
                "--continuous",
            ])
            .is_err()
        );
    }
}
