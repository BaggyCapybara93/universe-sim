//! Fixed workloads for reproducible performance comparisons.
//! Run: cargo run --release --example benchmark
use std::time::{Duration, Instant};

use rayon::ThreadPoolBuilder;
use universe_sim::{SimulationConfig, TickDuration, Universe};

#[derive(Clone, Copy)]
struct Case {
    name: &'static str,
    seed: u64,
    systems: u32,
    ticks: u64,
}

fn fingerprint(universe: &Universe) -> u64 {
    // FNV-1a over the complete debug state, calculated outside timed sections.
    // Use only to compare equivalent runs of the same source revision.
    format!("{universe:?}")
        .bytes()
        .fold(0xcbf29ce484222325, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
        })
}

fn run(case: Case, threads: usize) -> (Duration, Duration, u64) {
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    pool.install(|| {
        let mut config = SimulationConfig::new(case.seed, "1 kyr".parse::<TickDuration>().unwrap());
        config.generation.star_systems = case.systems;
        let start = Instant::now();
        let mut universe = Universe::new(config).unwrap();
        let generation = start.elapsed();
        let start = Instant::now();
        universe.advance_ticks(case.ticks).unwrap();
        let advancement = start.elapsed();
        (generation, advancement, fingerprint(&universe))
    })
}

fn main() {
    let smoke = std::env::args().skip(1).any(|arg| arg == "--smoke");
    let cases = if smoke {
        vec![Case {
            name: "smoke",
            seed: 42,
            systems: 8,
            ticks: 10_000,
        }]
    } else {
        vec![
            Case {
                name: "small",
                seed: 42,
                systems: 64,
                ticks: 10_000,
            },
            Case {
                name: "large",
                seed: 918_273_645,
                systems: 512,
                ticks: 1_000_000,
            },
        ]
    };
    println!("case,seed,systems,ticks,threads,generation_ms,advancement_ms,state_fingerprint");
    for case in cases {
        let mut baseline = None;
        for threads in [1, 4] {
            let (generation, advancement, digest) = run(case, threads);
            if let Some(expected) = baseline {
                assert_eq!(
                    digest, expected,
                    "state differs by worker count for {}",
                    case.name
                );
            } else {
                baseline = Some(digest);
            }
            println!(
                "{},{},{},{},{},{:.3},{:.3},{digest:016x}",
                case.name,
                case.seed,
                case.systems,
                case.ticks,
                threads,
                generation.as_secs_f64() * 1_000.0,
                advancement.as_secs_f64() * 1_000.0,
            );
        }
    }
}
