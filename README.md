# Universe Sim

Universe Sim is an early, deterministic astrophysical simulation engine written
in Rust. It creates a sparse catalog of galaxies containing deterministic
regions, star systems, stars, and planets from a seed, advances one canonical
simulation clock, and atomically saves or restores the resulting state.

The simulation core uses Rayon for deterministic data parallelism. The
command-line executable uses Clap for typed argument parsing, validation, and
generated help.

## Try it

Create a universe from the tracked configuration:

```console
cargo run -- new --config cosmos.toml --output universe.state
```

Override the initial number of star systems from either the defaults or the
configuration file:

```console
cargo run -- new --seed 42 --star-systems 1000 --output universe.state
```

Create a larger sparse universe without materializing every galaxy immediately:

```console
cargo run -- new --seed 42 --galaxies 64 --star-systems 1000 --output universe.state
cargo run -- galaxies --state universe.state
```

`--star-systems` controls the population batch in each materialized region, not
a universe-wide ceiling. Only the primary galaxy's origin is generated during
`new`; other galaxy regions are generated when requested. The per-region safety
limit remains 100,000 systems, while aggregate potential population grows with
the number of materialized galaxy regions.

The generated star, planet, and moon totals are reported separately because
binary or triple companions, planet occurrence, and moon occurrence are
procedural. Asteroid belts, comet clouds, and young-system debris disks are
stored as aggregate populations with estimated body counts and total mass, so
they do not require materializing billions of individual small bodies. The CLI
also reports the resulting total number of physical bodies and all stored
objects.
Initial stellar masses may extend to 300 solar masses, the simulator's enforced
upper physical limit. The mass function strongly favors smaller stars, so very
massive stars remain rare unless a large population is generated.

Advance it by 10,000 universal ticks:

```console
cargo run -- simulate --state universe.state --ticks 10000 --checkpoint-every 1000
```

The simulator reports stellar stage changes that occur during the interval.
Inspect the retained history with:

```console
cargo run -- events --state universe.state
```

Run continuously in visible checkpoint batches until Ctrl+C:

```console
cargo run -- simulate \
  --state universe.state \
  --continuous \
  --checkpoint-every 1000 \
  --pace-ms 1000
```

Continuous mode defaults to 1,000 ticks and a one-second pause per checkpoint.
Ctrl+C does not perform an extra shutdown save, but every completed checkpoint is
written atomically and remains valid. Events are evaluated for each simulated
interval and printed immediately after the checkpoint in which they occur.
Normal simulation output contains only events triggered during that run. Tick,
elapsed-time, checkpoint, completion, and state-path messages are shown only
with `--verbose`, which also expands planetary impact details. The explicit
`events` command displays saved history.

Large finite jumps such as `--ticks 9999999999` intentionally use analytical
catch-up rather than executing billions of empty loop iterations. Runtime scales
with materialized objects and events, not merely with the numeric tick count.

Inspect the saved state:

```console
cargo run -- inspect --state universe.state
```

Inspect a particular generated object using an ID shown in the state or another
inspection result:

```console
cargo run -- inspect --state universe.state S-123456789
```

Add `--verbose` to any command to emit structured diagnostic records to standard
error.

CPU-heavy generation and advancement use Rayon when the workload is large
enough. The default worker count follows the CPUs available to the process and
can be overridden without changing deterministic state:

```console
cargo run -- --threads 8 new --seed 42 --star-systems 10000
cargo run -- --threads 8 simulate --state universe.state --ticks 1000000
```

Small workloads remain serial to avoid thread-pool overhead. Use `--threads 1`
as a deterministic performance baseline.

Simulation advances every materialized galaxy in deterministic batches of at
most four galaxies. Select a smaller batch when individual galaxies are large:

```console
cargo run -- --threads 8 simulate --state universe.state --ticks 1000000 --galaxies-at-once 2
```

`--galaxies-at-once` accepts 1 through 4 and defaults to 4. It limits how many
galaxies are processed concurrently; it does not skip or freeze the remaining
materialized galaxies.

## Procedural regions

The coordinate space is unbounded. Regions are generated deterministically only
when requested, so the saved state remains finite:

```console
cargo run -- materialize --state universe.state --region 18,-4,91
cargo run -- regions --state universe.state
cargo run -- inspect --state universe.state --region 18,-4,91 S-123456789
cargo run -- unload --state universe.state --region 18,-4,91
```

Regions are scoped by galaxy, so the same coordinate may be independently
materialized in multiple galaxies:

```console
cargo run -- galaxies --state universe.state
cargo run -- materialize --state universe.state --galaxy G-123 --region 18,-4,91
cargo run -- inspect --state universe.state --galaxy G-123 --region 18,-4,91
cargo run -- regions --state universe.state --galaxy G-123
```

Unloading removes the detailed object graph while retaining the region's sparse
history. Materializing it again regenerates the same baseline and catches its
stars and planetary orbits up to the universal clock.

Run the validation suite:

```console
cargo test
```

## Stellar feedback and chemical evolution

Stellar mass lost during lifecycle transitions returns immediately to the local
region's gas reservoir, including winds, envelope shedding, and explosive
disruption. For core-collapse supernovae, the ejected mass is deposited at the
transition to the neutron-star or black-hole remnant, where the stellar mass
prescription applies the loss. Binary mergers and non-retained mass transfer
also return their lost mass. Gas mixes by mass across hydrogen, helium, carbon,
oxygen, iron, and other metals. Stellar stage changes cannot recreate mass
previously lost or transferred.

New stars inherit all six gas fractions at their formation time, with metallicity
derived from those fractions. Their birth composition is saved and shown by
`inspect --state universe.state S-…`; later planets use that stored composition.
Star formation retains unused gas and subtracts the mass placed in the star and
its disk. Feedback is processed chronologically before each formation interval,
so a large time jump can recycle ejecta into subsequent generations. Empty
formation intervals are skipped until stellar evolution can replenish the gas.

Enrichment uses an intentionally simple, uncalibrated population prescription:
giant/envelope ejecta convert 10% of their initial hydrogen to helium, then 5%
of the resulting helium to metals; massive-star winds use 20% and 8%; explosive
remnant/disruption ejecta use 35% and 20%. Added metals follow a fixed C/O/Fe/other
split of 18/44/10/28%. Other losses retain the birth composition. These conversions
preserve the modeled ejecta mass; nuclear binding-energy losses are omitted.
This is instantaneous local mixing, not an expanding supernova-remnant or
hydrodynamics model. Whole-universe matter accounting for the existing planetary
formation, disk dissipation, and destruction prescriptions remains future work.

State format 13 stores stellar birth compositions and versions the new evolution
rules. Older states reconstruct birth compositions from stored metallicity and
continue under the new rules; historical mass loss is not deposited again.

## Current boundaries

- The configuration reader supports the settings demonstrated in
  `cosmos.toml`. It is deliberately not a complete TOML implementation.
- State files carry explicit format, simulation, generation, and RNG versions.
  The earlier clock-only version 1 and object versions 2 through 12 are migrated;
  unknown versions are rejected rather than silently misinterpreted.
- Galactic populations use disk, bulge, and halo spatial/velocity profiles,
  correlated age and metallicity, a broken-power-law stellar mass function,
  mass-dependent binary and occasional triple occurrence, metallicity-weighted planet occurrence,
  and geometrically stable planetary spacing.
- Galaxy catalogs are deterministic and sparse. Spiral, elliptical, irregular,
  and dwarf entries have independent IDs, positions, sizes, population mixes,
  gas budgets, and region seeds. Morphology biases stellar ages, metallicities,
  and ongoing star formation without materializing the full galaxy.
- Planets store full orbital elements and move through analytical Kepler
  propagation, making even enormous time jumps proportional to object count
  rather than tick count. Companion stars have persisted hierarchical orbital
  elements. Systems whose companion orbits intrude on planetary stability zones
  use the symplectic leapfrog N-body path for short resolvable windows; stable
  systems and very long catch-up intervals remain analytical.
- Stars now follow a complete single-star lifecycle from protostar and
  pre-main-sequence contraction through hydrogen burning and mass-dependent
  post-main-sequence phases. Lower-mass stars pass through subgiant, giant,
  core-helium-burning, asymptotic-giant-branch, planetary-nebula, and white-dwarf
  stages. Massive stars become supergiants and, where winds strip the envelope,
  Wolf-Rayet stars before core-collapse supernovae, neutron stars, black holes,
  pulsational pair instability, complete pair-instability disruption, or direct
  collapse. Lifetimes, winds, temperatures, natal kicks, and remnant masses are
  deterministic population prescriptions rather than a stellar-interior or
  hydrodynamics solver.
- Stellar expansion, explosive collapse, and rapid mass loss dynamically
  irradiate, destroy, or unbind planets. Those outcomes are attached to the
  triggering event, persisted as planet status, and stop destroyed or unbound
  planets from being propagated as stable Kepler orbits.
- Every generated region contains a finite persisted star-forming gas reservoir.
  Its normalized hydrogen, helium, carbon, oxygen, iron, and other-metal mass
  fractions vary with galaxy morphology and survive save/load cycles.
  At deterministic collapse intervals, a fixed formation efficiency limits the
  available stellar mass; only material placed in the new star and disk is
  removed. New stars begin as age-zero protostars and inherit the local gas
  composition. Stellar ejecta replenish and enrich the reservoir over time.
- Evolved stars in binaries are checked against their periastron Roche lobes.
  Stable overflow transfers mass with partial retention; contact systems and
  sufficiently unstable mass ratios merge, remove the consumed companion, and
  rejuvenate the merger product.
- Physical radii are stage- and mass-dependent. A planet whose periapsis enters
  its host's radius is destroyed as a collision, while binary components whose
  radii overlap at periastron merge. These are population prescriptions rather
  than hydrodynamic collision calculations.
- Planets can have deterministic individual moons with persisted parent
  relationships and planet-centered Keplerian orbits. Systems may also contain
  aggregate asteroid belts, comet clouds, and debris disks with radial ranges,
  estimated populations, and total masses. Small-body aggregates are currently
  descriptive populations rather than individually integrated bodies.
- Planets include deterministic radius, broad atmosphere classification,
  equilibrium temperature, and a deliberately conservative potential-
  habitability heuristic. Planet temperatures are refreshed when stellar
  luminosity changes; this is not yet a climate or atmospheric-escape solver.
