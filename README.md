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

Stellar mass lost during lifecycle transitions returns immediately to its source
cloud, including winds, envelope shedding, and explosive
disruption. For core-collapse supernovae, the ejected mass is deposited at the
transition to the neutron-star or black-hole remnant, where the stellar mass
prescription applies the loss. Binary mergers and non-retained mass transfer
also return their lost mass. Gas mixes by mass within each cloud across hydrogen,
helium, carbon, oxygen, iron, and other metals. Stellar stage changes cannot recreate mass
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

State format 13 stores stellar birth compositions. Older states reconstruct
birth compositions from stored metallicity; historical mass loss is not
deposited again.

## Compact-binary outcomes

White-dwarf, neutron-star, and black-hole companions now receive a deterministic
gravitational-wave coalescence time once both stars are compact. The model uses
the quadrupole circular-orbit inspiral time with an approximate eccentricity
factor. It schedules the merger chronologically alongside stellar transitions,
including during large time jumps, and records the companion ID in the event.
`simulate` prints new events and `events` displays them from saved history.

- White-dwarf pairs above 1.44 solar masses produce a Type Ia supernova;
  lighter pairs leave a white dwarf. A white dwarf also detonates if binary
  mass transfer takes it above that threshold.
- Neutron-star pairs and neutron-star/black-hole pairs produce a kilonova and
  leave a neutron star or black hole according to the retained mass.
- Black-hole pairs leave a black hole and produce a black-hole merger event.
  Other mixed compact pairs produce a compact-object merger event.

These are population prescriptions. Type Ia ejecta return to the source cloud
with a fixed iron-rich mixture and can trigger a cloud shock; kilonova ejecta
use a fixed heavy-element mixture. Black-hole pairs lose 5% of their mass as
gravitational-wave energy. The model schedules the merger from stored orbital
elements but does not continuously shrink the displayed orbit or resolve the
transient light curve. Generated binaries are usually wider than the compact
orbits needed for a merger within a typical run.

State format 16 versions the compact-binary rules. Version 15 states load and
continue with those outcomes enabled.

## Galaxy encounters

Galaxy catalog entries now move on deterministic straight-line trajectories.
When two active galaxies enter a close-passage distance of three times their
combined radii, the simulator records an encounter. A later approach within
half their combined radii merges them. The primary galaxy survives a merger
involving it; otherwise the lower-ID galaxy survives. The remnant gains a
volume-combined radius, a weighted velocity, and an elliptical morphology.

Close passages temporarily raise the share of available cloud gas eligible for
star formation by 50%; mergers raise it by 100%. The boost lasts 200 million
years, works across large time jumps, and is replayed when a previously known
region is materialized. It draws from stored gas rather than creating mass.

Each galaxy also has a deterministic central black hole whose seed mass depends
on galaxy size and morphology. Black holes accrete steadily at a low rate and
enter brighter, fading AGN episodes after close passages and mergers. The
`galaxies` command and galaxy inspection report black-hole mass and spin,
accretion rate, Eddington ratio, luminosity, and AGN phase. Active AGN reduce
cloud star-formation efficiency. Accretion and feedback use population-level
prescriptions; accreted mass is not drawn from a spatial gas reservoir.

```console
cargo run -- galaxies --state universe.state
cargo run -- encounters --state universe.state
```

`simulate` also prints new encounters. Absorbed galaxies remain in the catalog
as historical progenitors; their existing regions and objects retain their
original IDs and remain inspectable. Previously known regions can be unloaded
and rematerialized; new regions must use the surviving galaxy's ID. Galaxy
trajectories describe catalog-scale motion; the
stored coordinates of individual systems remain in their original local
galaxy frames. The straight-line encounter thresholds and starburst factors are
population prescriptions, not a gravitational N-body model.

State format 18 stores galaxy motion, merger status, and encounter history.
It also introduces cluster member motion; older saved positions are updated to
the new prescription when loaded. Older states start galaxy motion at their
saved clock time.

## Current boundaries

- The configuration reader supports the settings demonstrated in
  `cosmos.toml`. It is deliberately not a complete TOML implementation.
- State files carry explicit format, simulation, generation, and RNG versions.
  The earlier clock-only version 1 and object versions 2 through 16 are migrated;
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
- Every generated region contains persistent diffuse and molecular gas clouds
  with independent mass and six-component composition. The initial mixture
  varies with galaxy morphology. At deterministic collapse intervals, a fixed
  efficiency limits how much gas from one molecular cloud becomes stars and
  disks. Up to four coeval systems form in a compact, persistent cluster;
  their age-zero protostars inherit that cloud's composition. Stellar ejecta
  replenish and enrich the source cloud. Existing generated systems receive
  individual historical clusters so their ancestry remains inspectable.
  Follow a planet's system to its cluster and source cloud with `inspect`.
  `inspect` also reports a cluster's modeled phase, characteristic radius,
  and estimated number of bound systems at the current simulation time. The
  association stays embedded for its first 3 Myr, then expands and gradually
  loses bound members on a 50 Myr evaporation timescale. Member systems move
  outward from the birth center and drift faster after escaping. Their planets
  and companions follow the updated system center during orbital propagation;
  birth-cluster ancestry remains inspectable. This is a deterministic
  population prescription rather than an N-body cluster simulation.
  Core-collapse and pair-instability supernova ejecta can compress the source
  cloud and nearby clouds. A pending shock increases the next collapse's gas
  budget and sibling-system limit; the resulting cluster records its triggering
  star. Shock state and its causal link survive save/load cycles.
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

## Validation and performance baselines

`ObjectRegistry::matter_inventory()` reports gas, stars, active protoplanetary
disks, planets, moons, individual asteroids and comets, and aggregate small-body
mass for a materialized region. `Universe::matter_inventories()` covers every
materialized region, including the primary origin. Aggregate populations are
listed separately because their masses can overlap the sampled individual
bodies. Gas metal mass is exact within the stored six-component reservoir;
stellar and disk metal values are birth-composition proxies. The inventory is
an observable budget, not a claim of whole-universe conservation: planetary
formation, disk dissipation, binary-accreted chemistry, and the composition of
stellar interiors still need fully closed material ledgers.

The validation suite checks that stellar feedback conserves stellar-plus-gas
mass over 12 billion years, that a new star and its disk are funded by gas,
and that the modeled metal proxy grows when ejecta are enriched. It also
checks orbital energy and angular momentum over a billion-year analytical
Kepler orbit and a thousand-year leapfrog integration. A fixed-seed sampling
test checks the two *conditional* stellar mass-function slopes against the
approximate low-mass and high-mass indices summarized in
[Kroupa, *The Local Stellar Initial Mass Function*](https://arxiv.org/abs/astro-ph/0011328).
It does not establish normalization across the mass-function break or validate
supernova rates, remnant ratios, and chemical distributions against observations.

Run the reproducible workload suite with:

```console
cargo run --release --locked --example benchmark
```

It emits CSV with generation and advancement times and a full-state fingerprint.
Each fixed workload runs with one and four Rayon threads; a fingerprint mismatch
fails the run. A shorter `--smoke` case is available for a quick local check.
The reference below was measured on macOS 27.0.1, arm64, Rust 1.98.0 in release mode. Millisecond
values are local observations, not pass/fail thresholds; compare runs on the
same machine and toolchain.

| Case | Seed | Systems | Ticks | Threads | Generation ms | Advancement ms | Fingerprint |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| small | 42 | 64 | 10,000 | 1 | 1.035 | 0.906 | `278bc6bc44c264b6` |
| small | 42 | 64 | 10,000 | 4 | 0.678 | 0.731 | `278bc6bc44c264b6` |
| large | 918273645 | 512 | 1,000,000 | 1 | 5.249 | 36.855 | `8d7598904e05ca3d` |
| large | 918273645 | 512 | 1,000,000 | 4 | 2.714 | 29.564 | `8d7598904e05ca3d` |
