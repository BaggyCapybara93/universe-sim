use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::Duration;

static TEST_SEQUENCE: AtomicU64 = AtomicU64::new(0);

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new() -> Self {
        let sequence = TEST_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "universe-sim-cli-{}-{sequence}",
            std::process::id()
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_universe-sim"))
        .args(arguments)
        .output()
        .unwrap()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8(bytes.to_vec()).unwrap()
}

#[test]
fn creates_advances_and_inspects_generated_objects() {
    let directory = TestDirectory::new();
    let state = directory.path().join("universe.state");
    let state_text = state.to_str().unwrap();

    let created = run(&[
        "--verbose",
        "--threads",
        "2",
        "new",
        "--seed",
        "42",
        "--tick-length",
        "1 kyr",
        "--star-systems",
        "7",
        "--output",
        state_text,
    ]);
    assert!(created.status.success(), "{}", text(&created.stderr));
    assert!(text(&created.stdout).contains("Star systems: 7"));
    assert!(text(&created.stdout).contains("Gas clouds:"));
    assert!(text(&created.stdout).contains("Star clusters: 7"));
    assert!(text(&created.stdout).contains("Gas composition: H"));
    assert!(text(&created.stdout).contains("Total bodies:"));
    assert!(text(&created.stdout).contains("Moons:"));
    assert!(text(&created.stderr).contains("event=state_saved"));
    assert!(text(&created.stderr).contains("threads=2"));

    let advanced = run(&[
        "simulate",
        "--state",
        state_text,
        "--ticks",
        "100",
        "--checkpoint-every",
        "25",
    ]);
    assert!(advanced.status.success(), "{}", text(&advanced.stderr));
    assert!(text(&advanced.stdout).trim().is_empty());

    let universe_inspection = run(&["inspect", "--state", state_text]);
    assert!(universe_inspection.status.success());
    assert!(text(&universe_inspection.stdout).contains("Elapsed time: 100 kyr"));

    let contents = fs::read_to_string(&state).unwrap();
    let star_record = contents
        .lines()
        .find(|line| line.starts_with("S|"))
        .unwrap();
    let star_id = format!("S-{}", star_record.split('|').nth(1).unwrap());
    let inspected = run(&["inspect", "--state", state_text, &star_id]);
    assert!(inspected.status.success(), "{}", text(&inspected.stderr));
    assert!(text(&inspected.stdout).contains(&format!("STAR {star_id}")));
    assert!(text(&inspected.stdout).contains("Age at epoch:"));
    assert!(text(&inspected.stdout).contains("Current mass:"));
    assert!(text(&inspected.stdout).contains("Birth composition: H"));
    assert!(text(&inspected.stdout).contains("Source cluster: CL-"));
    assert!(text(&inspected.stdout).contains("Source cloud: GC-"));
    let cloud_record = contents
        .lines()
        .find(|line| line.starts_with("H|"))
        .unwrap();
    let cloud_id = format!("GC-{}", cloud_record.split('|').nth(1).unwrap());
    let inspected_cloud = run(&["inspect", "--state", state_text, &cloud_id]);
    assert!(
        inspected_cloud.status.success(),
        "{}",
        text(&inspected_cloud.stderr)
    );
    assert!(text(&inspected_cloud.stdout).contains(&format!("GAS CLOUD {cloud_id}")));
    let cluster_record = contents
        .lines()
        .find(|line| line.starts_with("L|"))
        .unwrap();
    let cluster_id = format!("CL-{}", cluster_record.split('|').nth(1).unwrap());
    let inspected_cluster = run(&["inspect", "--state", state_text, &cluster_id]);
    assert!(
        inspected_cluster.status.success(),
        "{}",
        text(&inspected_cluster.stderr)
    );
    assert!(text(&inspected_cluster.stdout).contains(&format!("STAR CLUSTER {cluster_id}")));

    let materialized = run(&["materialize", "--state", state_text, "--region", "18,-4,91"]);
    assert!(
        materialized.status.success(),
        "{}",
        text(&materialized.stderr)
    );
    let materialized_output = text(&materialized.stdout);
    assert!(materialized_output.contains("Region 18,-4,91 materialized"));
    let region_star = materialized_output
        .lines()
        .find_map(|line| line.trim().strip_prefix("Example star: "))
        .unwrap();

    let region_inspection = run(&[
        "inspect",
        "--state",
        state_text,
        "--region",
        "18,-4,91",
        region_star,
    ]);
    assert!(
        region_inspection.status.success(),
        "{}",
        text(&region_inspection.stderr)
    );
    assert!(text(&region_inspection.stdout).contains(region_star));

    let unloaded = run(&["unload", "--state", state_text, "--region", "18,-4,91"]);
    assert!(unloaded.status.success(), "{}", text(&unloaded.stderr));
    let regions = run(&["regions", "--state", state_text]);
    assert!(regions.status.success(), "{}", text(&regions.stderr));
    assert!(text(&regions.stdout).contains("18,-4,91: unloaded"));

    let temporary_files: Vec<_> = fs::read_dir(directory.path())
        .unwrap()
        .filter_map(Result::ok)
        .filter(|entry| entry.file_name().to_string_lossy().ends_with(".tmp"))
        .collect();
    assert!(temporary_files.is_empty());
}

#[test]
fn creates_and_materializes_a_secondary_galaxy() {
    let directory = TestDirectory::new();
    let state = directory.path().join("multi-galaxy.state");
    let state_text = state.to_str().unwrap();
    let created = run(&[
        "new",
        "--seed",
        "42",
        "--galaxies",
        "3",
        "--star-systems",
        "2",
        "--output",
        state_text,
    ]);
    assert!(created.status.success(), "{}", text(&created.stderr));

    let listed = run(&["galaxies", "--state", state_text]);
    assert!(listed.status.success(), "{}", text(&listed.stderr));
    let listed_text = text(&listed.stdout);
    let galaxy_ids: Vec<_> = listed_text
        .lines()
        .filter_map(|line| line.trim().split(':').next())
        .filter(|value| value.starts_with("G-"))
        .collect();
    assert_eq!(galaxy_ids.len(), 3);

    let materialized = run(&[
        "materialize",
        "--state",
        state_text,
        "--galaxy",
        galaxy_ids[1],
        "--region",
        "0,0,0",
    ]);
    assert!(
        materialized.status.success(),
        "{}",
        text(&materialized.stderr)
    );
    assert!(text(&materialized.stdout).contains(galaxy_ids[1]));

    let regions = run(&["regions", "--state", state_text, "--galaxy", galaxy_ids[1]]);
    assert!(regions.status.success(), "{}", text(&regions.stderr));
    assert!(text(&regions.stdout).contains(&format!("{}@0,0,0", galaxy_ids[1])));
}

#[test]
fn lists_persisted_galaxy_encounters_and_mergers() {
    let directory = TestDirectory::new();
    let state = directory.path().join("encounters.state");
    let state_text = state.to_str().unwrap();
    let created = run(&[
        "new",
        "--seed",
        "42",
        "--galaxies",
        "8",
        "--star-systems",
        "1",
        "--tick-length",
        "1 Myr",
        "--output",
        state_text,
    ]);
    assert!(created.status.success(), "{}", text(&created.stderr));

    let advanced = run(&["simulate", "--state", state_text, "--ticks", "10000"]);
    assert!(advanced.status.success(), "{}", text(&advanced.stderr));
    assert!(text(&advanced.stdout).contains("merger"));

    let encounters = run(&["encounters", "--state", state_text]);
    assert!(encounters.status.success(), "{}", text(&encounters.stderr));
    let listed = text(&encounters.stdout);
    assert!(listed.contains("close passage"));
    assert!(listed.contains("merger"));

    let galaxies = run(&["galaxies", "--state", state_text]);
    assert!(galaxies.status.success(), "{}", text(&galaxies.stderr));
    assert!(text(&galaxies.stdout).contains("merged into"));
}

#[test]
fn rejects_corrupted_or_incompatible_state_through_the_cli() {
    let directory = TestDirectory::new();
    let state = directory.path().join("universe.state");
    let state_text = state.to_str().unwrap();
    let created = run(&["new", "--seed", "7", "--output", state_text]);
    assert!(created.status.success(), "{}", text(&created.stderr));

    let contents = fs::read_to_string(&state)
        .unwrap()
        .replace("state_format_version=17", "state_format_version=99");
    fs::write(&state, contents).unwrap();
    let inspected = run(&["inspect", "--state", state_text]);
    assert!(!inspected.status.success());
    assert!(text(&inspected.stderr).contains("incompatible state"));
}

#[test]
fn continuous_mode_leaves_a_valid_completed_checkpoint_when_stopped() {
    let directory = TestDirectory::new();
    let state = directory.path().join("universe.state");
    let state_text = state.to_str().unwrap();
    let created = run(&[
        "new",
        "--seed",
        "99",
        "--tick-length",
        "1 year",
        "--output",
        state_text,
    ]);
    assert!(created.status.success(), "{}", text(&created.stderr));

    let mut child = Command::new(env!("CARGO_BIN_EXE_universe-sim"))
        .args([
            "--verbose",
            "simulate",
            "--state",
            state_text,
            "--continuous",
            "--checkpoint-every",
            "10",
            "--pace-ms",
            "10",
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();

    let mut advanced = false;
    for _ in 0..100 {
        thread::sleep(Duration::from_millis(20));
        let contents = fs::read_to_string(&state).unwrap();
        advanced = contents
            .lines()
            .find_map(|line| line.strip_prefix("tick="))
            .and_then(|tick| tick.parse::<u64>().ok())
            .is_some_and(|tick| tick > 0);
        if advanced {
            break;
        }
    }
    child.kill().unwrap();
    child.wait().unwrap();
    assert!(advanced, "continuous mode did not complete a checkpoint");

    let inspected = run(&["inspect", "--state", state_text]);
    assert!(inspected.status.success(), "{}", text(&inspected.stderr));
}
