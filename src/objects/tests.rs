use super::*;
use crate::config::GenerationConfig;
use crate::id::GalaxyId;
use crate::orbital::OrbitalElements;
use crate::units::{Angle, Length, Mass, Position};

#[test]
fn generation_is_deterministic_and_relationships_are_valid() {
    let first = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
    let second = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
    assert_eq!(first, second);
    first.validate().unwrap();
    assert_eq!(first.counts().galaxies, 1);
    assert_eq!(first.counts().star_systems, 4);
    assert!(first.counts().stars >= 4);
    assert!(first.counts().stars <= 12);
}

#[test]
fn generation_adds_moons_and_aggregate_small_body_populations() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 128;
    let objects = ObjectRegistry::generate(42, config).unwrap();
    let counts = objects.counts();
    assert!(counts.moons > 0);
    assert!(counts.asteroid_belts > 0);
    assert!(counts.comet_clouds > 0);
    assert!(counts.debris_disks > 0);
    assert!(objects.planets().all(|planet| {
        planet.radius().meters() > 0.0 && planet.equilibrium_temperature().kelvin() > 0.0
    }));
    for moon in objects.moons() {
        let planet = objects.planet(moon.planet_id()).unwrap();
        assert_eq!(planet.system_id(), moon.system_id());
        assert!(moon.radius().meters() > 0.0);
        assert!(moon.orbital_radius().meters() > 0.0);
        assert!(
            objects
                .star_system(moon.system_id())
                .unwrap()
                .moon_ids()
                .contains(&moon.id())
        );
    }
    for belt in objects.asteroid_belts() {
        assert!(belt.inner_radius() < belt.outer_radius());
        assert!(belt.estimated_body_count() > 0);
    }
    objects.validate().unwrap();
}

#[test]
fn different_seeds_generate_different_objects() {
    let first = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
    let second = ObjectRegistry::generate(43, GenerationConfig::DEFAULT).unwrap();
    assert_ne!(first, second);
}

#[test]
fn generation_supports_three_hundred_solar_mass_stars() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    config.planets_per_system_min = 0;
    config.planets_per_system_max = 0;
    config.star_mass_solar_min = 300.0;
    config.star_mass_solar_max = 300.0;
    let objects = ObjectRegistry::generate(42, config).unwrap();
    assert!(
        objects
            .stars()
            .all(|star| (star.initial_mass().solar_masses() - 300.0).abs() < 1e-9)
    );
}

#[test]
fn very_massive_low_metallicity_stars_cover_pair_instability_fates() {
    let mass = |solar| Mass::from_solar_masses(solar).unwrap();
    assert_eq!(
        stellar::stellar_fate(mass(100.0), 0.0001),
        stellar::StellarFate::PulsationalPairInstability
    );
    assert_eq!(
        stellar::stellar_fate(mass(180.0), 0.0001),
        stellar::StellarFate::PairInstability
    );
    assert_eq!(
        stellar::stellar_fate(mass(300.0), 0.0001),
        stellar::StellarFate::DirectCollapse
    );
    assert!(
        stellar::stellar_schedule(mass(100.0), 0.0001)
            .unwrap()
            .iter()
            .any(|(_, stage)| *stage == StellarStage::PulsationalPairInstability)
    );
    assert!(
        stellar::stellar_schedule(mass(180.0), 0.0001)
            .unwrap()
            .iter()
            .any(|(_, stage)| *stage == StellarStage::Disrupted)
    );
}

#[test]
fn solar_mass_star_has_a_complete_low_mass_lifecycle() {
    let initial_mass = Mass::from_solar_masses(1.0).unwrap();
    let schedule = stellar::stellar_schedule(initial_mass, 0.014).unwrap();
    let stages: Vec<_> = schedule.iter().map(|(_, stage)| *stage).collect();
    assert_eq!(
        stages,
        vec![
            StellarStage::PreMainSequence,
            StellarStage::MainSequence,
            StellarStage::Subgiant,
            StellarStage::Giant,
            StellarStage::CoreHeliumBurning,
            StellarStage::AsymptoticGiantBranch,
            StellarStage::PlanetaryNebula,
            StellarStage::WhiteDwarf,
        ]
    );
    assert!(
        schedule
            .windows(2)
            .all(|entries| entries[0].0 < entries[1].0)
    );
    assert_eq!(
        stellar::stellar_stage_at(initial_mass, 0.014, SimDuration::ZERO).unwrap(),
        StellarStage::Protostar
    );
    assert_eq!(
        stellar::stellar_stage_at(initial_mass, 0.014, schedule.last().unwrap().0).unwrap(),
        StellarStage::WhiteDwarf
    );
}

#[test]
fn massive_solar_metallicity_star_explodes_after_a_wolf_rayet_phase() {
    let schedule =
        stellar::stellar_schedule(Mass::from_solar_masses(30.0).unwrap(), 0.014).unwrap();
    let stages: Vec<_> = schedule.iter().map(|(_, stage)| *stage).collect();
    assert!(stages.contains(&StellarStage::Supergiant));
    assert!(stages.contains(&StellarStage::WolfRayet));
    assert!(stages.contains(&StellarStage::CoreCollapseSupernova));
    assert_eq!(stages.last(), Some(&StellarStage::BlackHole));
    assert!(
        schedule
            .windows(2)
            .all(|entries| entries[0].0 < entries[1].0)
    );
}

#[test]
fn gas_reservoir_forms_a_persisted_protostar_system() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    config.star_mass_solar_min = 0.1;
    config.star_mass_solar_max = 0.2;
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    let initial_systems = objects.counts().star_systems;
    let initial_gas = objects.gas_reservoir_solar_masses();
    let initial_composition = objects.gas_composition();
    let ten_myr = SimDuration::from_years(10_000_000).unwrap();
    let events = objects
        .form_stars_between(
            42,
            config,
            SimDuration::ZERO,
            ten_myr,
            SimDuration::from_years(1_000).unwrap(),
        )
        .unwrap();
    assert_eq!(events.len(), 1);
    assert_eq!(
        events[0].interaction,
        Some(StellarInteractionKind::StarFormation)
    );
    assert_eq!(objects.counts().star_systems, initial_systems + 1);
    assert!(objects.gas_reservoir_solar_masses() < initial_gas);
    assert_eq!(objects.gas_composition(), initial_composition);
    assert_eq!(
        objects.star(events[0].star_id).unwrap().formed_at(),
        ten_myr
    );
    let mut encoded = String::new();
    objects.encode(&mut encoded);
    assert_eq!(
        ObjectRegistry::decode(encoded.lines().skip(1)).unwrap(),
        objects
    );
}

#[test]
fn galaxy_morphology_sets_the_initial_gas_metallicity() {
    let expected = [
        (GalaxyMorphology::Spiral, 0.014),
        (GalaxyMorphology::Elliptical, 0.020),
        (GalaxyMorphology::Irregular, 0.006),
        (GalaxyMorphology::Dwarf, 0.002),
    ];
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    for (index, (morphology, metallicity)) in expected.into_iter().enumerate() {
        let galaxy = GalaxyCatalogEntry {
            id: GalaxyId::new(index as u64 + 1),
            position: Position::default(),
            morphology,
            radius_parsecs: config.galaxy_radius_parsecs,
        };
        let objects = ObjectRegistry::generate_in_galaxy(index as u64, config, galaxy).unwrap();
        assert!((objects.gas_composition().metallicity() - metallicity).abs() < 1.0e-12);
    }
}

#[test]
fn roche_overflow_transfers_mass_and_contact_merges() {
    let mut selected = (0..100)
        .find_map(|seed| {
            let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
            objects
                .star_systems()
                .find(|system| !system.companion_star_ids().is_empty())
                .map(|system| (objects.clone(), system.id()))
        })
        .expect("a deterministic binary system");
    let system = selected.0.star_system(selected.1).unwrap().clone();
    let primary_id = system.star_id;
    let companion_id = system.companion_star_ids[0];
    for id in [primary_id, companion_id] {
        let star = selected.0.stars.get_mut(&id).unwrap();
        star.current_mass = Mass::from_solar_masses(1.0).unwrap();
        star.initial_mass = star.current_mass;
        star.stage = StellarStage::MainSequence;
    }
    selected.0.stars.get_mut(&primary_id).unwrap().stage = StellarStage::Giant;
    selected
        .0
        .stars
        .get_mut(&companion_id)
        .unwrap()
        .orbital_elements = Some(test_orbit(1.0));
    let before = selected.0.star(primary_id).unwrap().current_mass();
    let outcomes = selected
        .0
        .resolve_binary_interactions(&BTreeMap::from([(primary_id, 1)]), SimDuration::ZERO)
        .unwrap();
    assert_eq!(
        outcomes[0].interaction,
        Some(StellarInteractionKind::MassTransfer)
    );
    assert!(selected.0.star(primary_id).unwrap().current_mass() < before);

    selected
        .0
        .stars
        .get_mut(&companion_id)
        .unwrap()
        .orbital_elements = Some(test_orbit(0.1));
    let outcomes = selected
        .0
        .resolve_binary_interactions(&BTreeMap::from([(primary_id, 2)]), SimDuration::ZERO)
        .unwrap();
    assert_eq!(
        outcomes[0].interaction,
        Some(StellarInteractionKind::StellarMerger)
    );
    assert!(selected.0.star(companion_id).is_none());
    selected.0.validate().unwrap();
}

#[test]
fn expanded_star_physically_consumes_intersecting_planet() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    config.planets_per_system_min = 1;
    config.planets_per_system_max = 1;
    config.orbit_au_min = 0.1;
    config.orbit_au_max = 0.1;
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    let star_id = objects.star_systems().next().unwrap().star_id();
    {
        let star = objects.stars.get_mut(&star_id).unwrap();
        star.current_mass = Mass::from_solar_masses(1.0).unwrap();
        star.stage = StellarStage::Giant;
    }
    let impacts = objects.apply_direct_planet_collisions(star_id, StellarStage::Giant);
    assert_eq!(impacts.len(), 1);
    assert_eq!(impacts[0].kind, PlanetImpactKind::Destroyed);
    assert_eq!(
        objects.planets().next().unwrap().status(),
        PlanetStatus::Destroyed
    );
}

fn test_orbit(semi_major_axis_au: f64) -> OrbitalElements {
    OrbitalElements::new(
        Length::from_astronomical_units(semi_major_axis_au).unwrap(),
        0.0,
        Angle::from_radians(0.0).unwrap(),
        Angle::from_radians(0.0).unwrap(),
        Angle::from_radians(0.0).unwrap(),
        Angle::from_radians(0.0).unwrap(),
        SimDuration::ZERO,
    )
    .unwrap()
}

#[test]
fn generated_companions_have_hierarchical_orbits() {
    let objects = (0..100)
        .find_map(|seed| {
            let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
            let has_companion = objects
                .star_systems()
                .any(|system| !system.companion_star_ids().is_empty());
            has_companion.then_some(objects)
        })
        .expect("at least one deterministic seed generates a companion");
    for system in objects
        .star_systems()
        .filter(|system| !system.companion_star_ids().is_empty())
    {
        let mut axes: Vec<_> = system
            .companion_star_ids()
            .iter()
            .map(|id| {
                objects
                    .star(*id)
                    .unwrap()
                    .orbital_elements()
                    .unwrap()
                    .semi_major_axis()
                    .astronomical_units()
            })
            .collect();
        axes.sort_by(f64::total_cmp);
        if axes.len() == 2 {
            assert!(axes[1] / axes[0] >= 10.0);
        }
    }
}

#[test]
fn interacting_systems_use_the_selective_nbody_path_for_short_windows() {
    let elapsed = SimDuration::from_years(1).unwrap();
    let mut selected = None;
    for seed in 0..200 {
        let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
        let interacting_system = objects.star_systems().find_map(|system| {
            (objects.dynamics_mode(system.id()) == Some(DynamicsMode::Interacting))
                .then_some(system.id())
        });
        if let Some(system_id) = interacting_system {
            selected = Some((objects, system_id));
            break;
        }
    }
    let (mut objects, system_id) = selected.expect("an interacting generated system");
    assert!(
        objects
            .integrate_system_from_epoch(system_id, elapsed)
            .unwrap()
    );
}

#[test]
fn object_records_round_trip_exactly() {
    let objects = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
    let mut encoded = String::new();
    objects.encode(&mut encoded);
    let decoded = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
    assert_eq!(decoded, objects);
}

#[test]
fn orbital_catch_up_changes_positions_without_changing_elements() {
    let mut objects = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
    let before: Vec<_> = objects
        .planets()
        .map(|planet| (planet.id(), planet.position(), planet.orbital_elements()))
        .collect();
    objects
        .advance_orbits_to(SimDuration::from_years(1_000_000).unwrap())
        .unwrap();
    assert!(
        before
            .iter()
            .any(|(id, position, _)| { objects.planet(*id).unwrap().position() != *position })
    );
    assert!(
        before.iter().all(|(id, _, elements)| {
            objects.planet(*id).unwrap().orbital_elements() == *elements
        })
    );
}

#[test]
fn parallel_generation_and_advancement_match_single_threaded_results() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 512;
    let run = || {
        let mut objects = ObjectRegistry::generate(42, config).unwrap();
        let elapsed = SimDuration::from_years(1_000_000_000).unwrap();
        let tick = SimDuration::from_years(1).unwrap();
        let transitions = objects.advance_stars_to(elapsed, tick).unwrap();
        objects.advance_orbits_to(elapsed).unwrap();
        (objects, transitions)
    };
    let single_thread = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(run);
    let four_threads = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap()
        .install(run);
    assert_eq!(single_thread, four_threads);
}
