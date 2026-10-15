use super::*;
use crate::config::GenerationConfig;
use crate::id::{GalaxyId, StarId};
use crate::orbital::OrbitalElements;
use crate::units::{Angle, Length, Mass, Position};

fn set_test_gas(objects: &mut ObjectRegistry, mass: f64, composition: crate::GasComposition) {
    let cloud = objects
        .gas_clouds
        .values_mut()
        .next()
        .expect("test cloud exists");
    cloud.mass_solar = mass;
    cloud.composition = composition;
    objects.gas_reservoir_solar_masses = mass;
    objects.gas_composition = composition;
}

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
fn clusters_expand_and_disperse_without_losing_birth_ancestry() {
    let objects = ObjectRegistry::generate(42, GenerationConfig::DEFAULT).unwrap();
    let cluster = objects.star_clusters().next().unwrap();
    let system_ids = cluster.system_ids().to_vec();
    let newborn = cluster.evolution_at(SimDuration::ZERO);
    let mature = cluster.evolution_at(SimDuration::from_years(20_000_000).unwrap());
    let old = cluster.evolution_at(SimDuration::from_years(500_000_000).unwrap());
    let mut association = cluster.clone();
    association.system_ids = objects.star_systems().map(StarSystem::id).collect();

    assert_eq!(newborn.phase(), ClusterPhase::Embedded);
    assert_eq!(newborn.bound_system_count(), system_ids.len());
    assert!(newborn.characteristic_radius_parsecs() < mature.characteristic_radius_parsecs());
    assert!(mature.characteristic_radius_parsecs() < old.characteristic_radius_parsecs());
    assert_eq!(old.phase(), ClusterPhase::Dispersed);
    assert_eq!(old.bound_system_count(), 0);
    assert_eq!(
        association
            .evolution_at(SimDuration::from_years(3_000_000).unwrap())
            .phase(),
        ClusterPhase::Bound
    );
    assert_eq!(
        association
            .evolution_at(SimDuration::from_years(20_000_000).unwrap())
            .phase(),
        ClusterPhase::Dispersing
    );
    assert_eq!(cluster.system_ids(), system_ids);
    assert_eq!(
        objects.star_system(system_ids[0]).unwrap().cluster_id(),
        Some(cluster.id())
    );

    let mut state = String::new();
    objects.encode(&mut state);
    let restored = ObjectRegistry::decode(state.lines().skip(1)).unwrap();
    assert_eq!(
        restored
            .star_cluster(cluster.id())
            .unwrap()
            .evolution_at(SimDuration::from_years(500_000_000).unwrap()),
        old
    );
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
    assert!(counts.asteroids > 0);
    assert!(counts.comets > 0);
    assert!(counts.protoplanetary_disks > 0);
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
    assert!(objects.asteroids().all(|asteroid| {
        objects
            .asteroid_belt(asteroid.belt_id())
            .is_some_and(|belt| belt.system_id() == asteroid.system_id())
    }));
    assert!(objects.comets().all(|comet| {
        objects
            .comet_cloud(comet.cloud_id())
            .is_some_and(|cloud| cloud.system_id() == comet.system_id())
    }));
    assert!(objects.planets().all(|planet| {
        let composition = planet.bulk_composition();
        (composition.rocky_fraction()
            + composition.icy_fraction()
            + composition.iron_fraction()
            + composition.hydrogen_helium_fraction()
            - 1.0)
            .abs()
            < 1.0e-12
    }));
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
    assert!(events.len() > 1);
    assert_eq!(
        events[0].interaction,
        Some(StellarInteractionKind::StarFormation)
    );
    assert_eq!(
        objects.counts().star_systems,
        initial_systems + events.len()
    );
    let cluster_id = objects
        .star_system(objects.star(events[0].star_id).unwrap().system_id())
        .unwrap()
        .cluster_id()
        .unwrap();
    let cluster = objects.star_cluster(cluster_id).unwrap();
    assert_eq!(cluster.system_ids().len(), events.len());
    assert!(events.iter().all(|event| {
        objects
            .star_system(objects.star(event.star_id).unwrap().system_id())
            .unwrap()
            .cluster_id()
            == Some(cluster_id)
    }));
    assert_eq!(
        objects.gas_cloud(cluster.cloud_id()).unwrap().kind(),
        GasCloudKind::Molecular
    );
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
fn cluster_members_spread_after_birth_and_match_checkpointed_positions() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    config.star_mass_solar_min = 0.1;
    config.star_mass_solar_max = 0.2;
    let mut birth = ObjectRegistry::generate(42, config).unwrap();
    let formed = SimDuration::from_years(10_000_000).unwrap();
    let events = birth
        .form_stars_between(
            42,
            config,
            SimDuration::ZERO,
            formed,
            SimDuration::from_years(1_000).unwrap(),
        )
        .unwrap();
    let system_id = birth.star(events[0].star_id).unwrap().system_id();
    let cluster_id = birth.star_system(system_id).unwrap().cluster_id().unwrap();
    let cluster = birth.star_cluster(cluster_id).unwrap().clone();
    assert!(cluster.system_ids().len() > 1);
    let center = cluster.position().meters();
    let birth_position = birth.star_system(system_id).unwrap().position().meters();

    let later = SimDuration::from_years(30_000_000).unwrap();
    let mut whole = birth.clone();
    whole.advance_orbits_to(later).unwrap();
    let evolved = whole.star_system(system_id).unwrap().position().meters();
    let distance = |point: [f64; 3]| {
        point
            .into_iter()
            .zip(center)
            .map(|(coordinate, origin)| (coordinate - origin).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    assert!(distance(evolved) > distance(birth_position));
    assert_eq!(
        whole.star(events[0].star_id).unwrap().position().meters(),
        evolved
    );
    assert!(
        whole
            .star_cluster(cluster_id)
            .unwrap()
            .evolution_at(later)
            .bound_system_count()
            < cluster.system_ids().len()
    );

    let mut resumed = birth;
    resumed
        .advance_orbits_to(SimDuration::from_years(20_000_000).unwrap())
        .unwrap();
    let mut state = String::new();
    resumed.encode(&mut state);
    let mut resumed = ObjectRegistry::decode(state.lines().skip(1)).unwrap();
    resumed.advance_orbits_to(later).unwrap();
    for system_id in cluster.system_ids() {
        assert_eq!(
            resumed.star_system(*system_id).unwrap().position(),
            whole.star_system(*system_id).unwrap().position()
        );
    }
}

#[test]
fn clouds_enrich_independently_and_planets_retain_cluster_ancestry() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 32;
    config.star_mass_solar_min = 0.1;
    config.star_mass_solar_max = 0.2;
    config.planets_per_system_min = 1;
    config.planets_per_system_max = 1;
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    assert_eq!(objects.counts().gas_clouds, 4);
    assert!(
        objects
            .gas_clouds()
            .any(|cloud| cloud.kind() == GasCloudKind::Diffuse)
    );
    let cloud_id = objects
        .gas_clouds()
        .find(|cloud| cloud.kind() == GasCloudKind::Molecular)
        .unwrap()
        .id();
    let original_other_clouds: Vec<_> = objects
        .gas_clouds()
        .filter(|cloud| cloud.id() != cloud_id)
        .cloned()
        .collect();
    let source_star_id = objects
        .star_systems()
        .find(|system| {
            objects
                .star_cluster(system.cluster_id().unwrap())
                .unwrap()
                .cloud_id()
                == cloud_id
        })
        .unwrap()
        .star_id();
    let ejecta = crate::GasComposition::from_metallicity(0.05).unwrap();
    objects.return_gas(source_star_id, 1.0, ejecta).unwrap();
    assert!(
        objects
            .gas_cloud(cloud_id)
            .unwrap()
            .composition()
            .metallicity()
            > 0.014
    );
    for original in original_other_clouds {
        assert_eq!(objects.gas_cloud(original.id()).unwrap(), &original);
    }
    let birth_composition = objects.gas_cloud(cloud_id).unwrap().composition();
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
    assert!(events.len() > 1);
    let cluster_id = objects
        .star_system(objects.star(events[0].star_id).unwrap().system_id())
        .unwrap()
        .cluster_id()
        .unwrap();
    let cluster = objects.star_cluster(cluster_id).unwrap();
    assert_eq!(cluster.cloud_id(), cloud_id);
    assert_eq!(cluster.system_ids().len(), events.len());
    for event in &events {
        let star = objects.star(event.star_id).unwrap();
        assert_eq!(star.birth_composition(), birth_composition);
        assert_eq!(
            objects.star_system(star.system_id()).unwrap().cluster_id(),
            Some(cluster_id)
        );
    }
    objects
        .advance_planetary_formation_to(42, config, SimDuration::from_years(15_000_000).unwrap())
        .unwrap();
    for system_id in objects.star_cluster(cluster_id).unwrap().system_ids() {
        let system = objects.star_system(*system_id).unwrap();
        assert_eq!(system.planet_ids().len(), 1);
        let planet = objects.planet(system.planet_ids()[0]).unwrap();
        assert_eq!(planet.system_id(), *system_id);
        assert_eq!(planet.star_id(), system.star_id());
    }
    objects.validate().unwrap();
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
            velocity: crate::units::Velocity::default(),
            morphology,
            radius_parsecs: config.galaxy_radius_parsecs,
            merged_into: None,
        };
        let objects = ObjectRegistry::generate_in_galaxy(index as u64, config, galaxy).unwrap();
        assert!((objects.gas_composition().metallicity() - metallicity).abs() < 1.0e-12);
    }
}

#[test]
fn galaxy_close_passage_boosts_star_formation_from_existing_gas() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    config.star_mass_solar_min = 1.0;
    config.star_mass_solar_max = 10.0;
    let mut baseline = ObjectRegistry::generate(42, config).unwrap();
    set_test_gas(&mut baseline, 10.0, crate::GasComposition::default());
    let mut tidal = baseline.clone();
    let galaxy_id = tidal.galaxies().next().unwrap().id();
    let encounter = crate::GalaxyEncounter::new(
        1,
        SimDuration::ZERO,
        0,
        galaxy_id,
        GalaxyId::new(999),
        crate::GalaxyEncounterKind::ClosePass,
    );
    let end = SimDuration::from_years(10_000_000).unwrap();
    let tick = SimDuration::from_years(1).unwrap();
    baseline
        .form_stars_between(42, config, SimDuration::ZERO, end, tick)
        .unwrap();
    tidal
        .form_stars_between_with_encounters(42, config, SimDuration::ZERO, end, tick, &[encounter])
        .unwrap();
    assert!(tidal.gas_reservoir_solar_masses() < baseline.gas_reservoir_solar_masses());
    tidal.validate().unwrap();
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
    let total_before = stellar_and_gas_mass(&selected.0);
    let outcomes = selected
        .0
        .resolve_binary_interactions(&BTreeMap::from([(primary_id, 1)]), SimDuration::ZERO)
        .unwrap();
    assert_eq!(
        outcomes[0].interaction,
        Some(StellarInteractionKind::MassTransfer)
    );
    assert!(selected.0.star(primary_id).unwrap().current_mass() < before);
    assert!((stellar_and_gas_mass(&selected.0) - total_before).abs() < 1e-12);

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
    assert!((stellar_and_gas_mass(&selected.0) - total_before).abs() < 1e-12);
    selected.0.validate().unwrap();
}

fn compact_binary_fixture(
    primary_stage: StellarStage,
    primary_mass: f64,
    companion_stage: StellarStage,
    companion_mass: f64,
) -> (ObjectRegistry, StarId, StarId) {
    let mut objects = (0..100)
        .find_map(|seed| {
            let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
            let has_companion = objects
                .star_systems()
                .any(|system| !system.companion_star_ids().is_empty());
            has_companion.then_some(objects)
        })
        .unwrap();
    let system = objects
        .star_systems()
        .find(|system| !system.companion_star_ids().is_empty())
        .unwrap()
        .clone();
    let companion_id = system.companion_star_ids[0];
    for (id, stage, mass) in [
        (system.star_id, primary_stage, primary_mass),
        (companion_id, companion_stage, companion_mass),
    ] {
        let star = objects.stars.get_mut(&id).unwrap();
        star.initial_mass = Mass::from_solar_masses(mass).unwrap();
        star.current_mass = star.initial_mass;
        star.stage = stage;
        star.formed_at = SimDuration::ZERO;
        star.epoch_age = stellar::stellar_schedule(star.initial_mass, star.metallicity)
            .unwrap()
            .into_iter()
            .find(|(_, next)| {
                matches!(
                    next,
                    StellarStage::WhiteDwarf | StellarStage::NeutronStar | StellarStage::BlackHole
                )
            })
            .unwrap()
            .0;
        star.age = star.epoch_age;
    }
    objects
        .stars
        .get_mut(&companion_id)
        .unwrap()
        .orbital_elements = Some(test_orbit(0.001));
    (objects, system.star_id, companion_id)
}

#[test]
fn compact_binary_inspiral_produces_distinct_outcomes() {
    for (primary_stage, primary_mass, companion_stage, companion_mass, outcome, remnant) in [
        (
            StellarStage::WhiteDwarf,
            0.8,
            StellarStage::WhiteDwarf,
            0.8,
            StellarInteractionKind::TypeIaSupernova,
            StellarStage::Disrupted,
        ),
        (
            StellarStage::NeutronStar,
            1.4,
            StellarStage::NeutronStar,
            1.4,
            StellarInteractionKind::Kilonova,
            StellarStage::BlackHole,
        ),
        (
            StellarStage::NeutronStar,
            1.2,
            StellarStage::NeutronStar,
            1.2,
            StellarInteractionKind::Kilonova,
            StellarStage::NeutronStar,
        ),
        (
            StellarStage::NeutronStar,
            1.4,
            StellarStage::BlackHole,
            5.0,
            StellarInteractionKind::Kilonova,
            StellarStage::BlackHole,
        ),
        (
            StellarStage::BlackHole,
            5.0,
            StellarStage::BlackHole,
            5.0,
            StellarInteractionKind::BlackHoleMerger,
            StellarStage::BlackHole,
        ),
        (
            StellarStage::WhiteDwarf,
            0.6,
            StellarStage::WhiteDwarf,
            0.6,
            StellarInteractionKind::CompactObjectMerger,
            StellarStage::WhiteDwarf,
        ),
        (
            StellarStage::WhiteDwarf,
            0.6,
            StellarStage::BlackHole,
            5.0,
            StellarInteractionKind::CompactObjectMerger,
            StellarStage::BlackHole,
        ),
    ] {
        let (mut objects, primary_id, companion_id) =
            compact_binary_fixture(primary_stage, primary_mass, companion_stage, companion_mass);
        let merger = objects
            .compact_merger_at(primary_id, companion_id)
            .unwrap()
            .unwrap();
        assert!(merger > 0 && merger < SimDuration::from_years(10_000_000).unwrap().seconds());
        let before = objects
            .advance_stars_to(
                SimDuration::from_seconds(merger - 1),
                SimDuration::from_years(1).unwrap(),
            )
            .unwrap();
        assert!(
            before
                .iter()
                .all(|event| event.related_star_id != Some(companion_id))
        );
        let events = objects
            .advance_stars_to(
                SimDuration::from_seconds(merger),
                SimDuration::from_years(1).unwrap(),
            )
            .unwrap();
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].interaction, Some(outcome));
        assert_eq!(events[0].related_star_id, Some(companion_id));
        assert_eq!(objects.star(primary_id).unwrap().stage(), remnant);
        assert!(objects.star(companion_id).is_none());
        if outcome == StellarInteractionKind::TypeIaSupernova {
            objects
                .advance_orbits_to(SimDuration::from_seconds(merger))
                .unwrap();
        }
        objects.validate().unwrap();
        assert!(
            objects
                .advance_stars_to(
                    SimDuration::from_years(10_000_000).unwrap(),
                    SimDuration::from_years(1).unwrap()
                )
                .unwrap()
                .is_empty()
        );
    }
}

#[test]
fn compact_binary_merger_is_checkpoint_and_save_independent() {
    let (mut whole, primary_id, companion_id) = compact_binary_fixture(
        StellarStage::NeutronStar,
        1.4,
        StellarStage::NeutronStar,
        1.4,
    );
    let mut split = whole.clone();
    let end = SimDuration::from_years(1_000_000).unwrap();
    let tick = SimDuration::from_years(1).unwrap();
    let whole_events = whole.advance_stars_to(end, tick).unwrap();
    split
        .advance_stars_to(SimDuration::from_years(10_000).unwrap(), tick)
        .unwrap();
    let mut encoded = String::new();
    split.encode(&mut encoded);
    split = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
    let split_events = split.advance_stars_to(end, tick).unwrap();
    assert_eq!(split, whole);
    assert_eq!(split_events, whole_events);
    assert_eq!(split_events[0].star_id, primary_id);
    assert_eq!(split_events[0].related_star_id, Some(companion_id));
}

#[test]
fn accreting_white_dwarf_detonates_at_chandrasekhar_mass() {
    let (mut objects, primary_id, companion_id) = compact_binary_fixture(
        StellarStage::WhiteDwarf,
        1.42,
        StellarStage::WhiteDwarf,
        1.0,
    );
    let total_before = stellar_and_gas_mass(&objects);
    let donor = objects.stars.get_mut(&companion_id).unwrap();
    donor.stage = StellarStage::Giant;
    donor.orbital_elements = Some(test_orbit(1.0));
    let events = objects
        .resolve_binary_interactions(&BTreeMap::from([(primary_id, 1)]), SimDuration::ZERO)
        .unwrap();
    assert!(events.iter().any(|event| event.interaction
        == Some(StellarInteractionKind::TypeIaSupernova)
        && event.star_id == primary_id));
    assert_eq!(
        objects.star(primary_id).unwrap().stage(),
        StellarStage::Disrupted
    );
    assert_eq!(objects.star(primary_id).unwrap().mass().solar_masses(), 0.0);
    assert!(objects.gas_composition().iron_fraction() > 0.0);
    assert!((stellar_and_gas_mass(&objects) - total_before).abs() < 1e-12);
    objects.validate().unwrap();
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
        before.iter().any(|(id, _, elements)| {
            objects.planet(*id).unwrap().orbital_elements() != *elements
        })
    );
}

#[test]
fn dynamic_protoplanetary_disk_forms_a_system_over_time() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 1;
    config.planets_per_system_min = 1;
    config.planets_per_system_max = 1;
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    let ten_myr = SimDuration::from_years(10_000_000).unwrap();
    objects
        .form_stars_between(
            42,
            config,
            SimDuration::ZERO,
            ten_myr,
            SimDuration::from_years(1_000).unwrap(),
        )
        .unwrap();
    let formed = SimDuration::from_years(15_000_000).unwrap();
    objects
        .advance_planetary_formation_to(42, config, formed)
        .unwrap();
    let dynamic_star = objects
        .stars()
        .find(|star| star.formed_at() == ten_myr)
        .unwrap();
    let system = objects.star_system(dynamic_star.system_id()).unwrap();
    assert_eq!(system.planet_ids().len(), 1);
    let planet = objects.planet(system.planet_ids()[0]).unwrap();
    assert_eq!(planet.formed_at(), formed);
    assert!(
        objects
            .protoplanetary_disks()
            .any(|disk| disk.system_id() == system.id())
    );
}

#[test]
fn planetary_systems_evolve_migration_tides_and_atmospheric_loss() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 128;
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    let planet_before: Vec<_> = objects
        .planets()
        .map(|planet| {
            (
                planet.id(),
                planet.orbital_radius(),
                planet.atmosphere_loss_fraction(),
            )
        })
        .collect();
    let moon_before: Vec<_> = objects
        .moons()
        .map(|moon| (moon.id(), moon.orbital_radius()))
        .collect();
    let asteroid_before: Vec<_> = objects
        .asteroids()
        .map(|asteroid| (asteroid.id(), asteroid.position()))
        .collect();
    objects
        .advance_orbits_to(SimDuration::from_years(1_000_000_000).unwrap())
        .unwrap();
    assert!(planet_before.iter().any(|(id, radius, loss)| {
        let planet = objects.planet(*id).unwrap();
        planet.orbital_radius() != *radius && planet.atmosphere_loss_fraction() >= *loss
    }));
    assert!(
        moon_before
            .iter()
            .any(|(id, radius)| { objects.moon(*id).unwrap().orbital_radius() != *radius })
    );
    assert!(
        asteroid_before
            .iter()
            .any(|(id, position)| { objects.asteroid(*id).unwrap().position() != *position })
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

// A controlled, isolated star with no pre-existing companions or planets.
fn feedback_fixture(mass_solar: f64, metallicity: f64) -> (ObjectRegistry, GenerationConfig) {
    let config = GenerationConfig {
        star_systems: 1,
        star_mass_solar_min: 0.1,
        star_mass_solar_max: 0.2,
        planets_per_system_min: 0,
        planets_per_system_max: 0,
        ..GenerationConfig::DEFAULT
    };
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    let system = objects.star_systems.values_mut().next().unwrap();
    system.companion_star_ids.clear();
    let star_id = system.star_id;
    let star = generation::generate_star_with_birth(
        star_id,
        system.id,
        SimDuration::ZERO,
        crate::GasComposition::from_metallicity(metallicity).unwrap(),
        mass_solar,
    )
    .unwrap();
    objects.stars.clear();
    objects.stars.insert(star_id, star);
    set_test_gas(
        &mut objects,
        0.0,
        crate::GasComposition::from_metallicity(metallicity).unwrap(),
    );
    (objects, config)
}

fn stellar_and_gas_mass(objects: &ObjectRegistry) -> f64 {
    objects.gas_reservoir_solar_masses()
        + objects
            .stars()
            .map(|star| star.mass().solar_masses())
            .sum::<f64>()
}

#[test]
fn stellar_feedback_conserves_mass_enriches_gas_and_is_not_replayed() {
    for (mass, metallicity) in [
        (1.0, 0.014),
        (20.0, 0.014),
        (60.0, 0.014),
        (180.0, 0.0001),
        (300.0, 0.0001),
    ] {
        let (mut objects, _) = feedback_fixture(mass, metallicity);
        let schedule =
            stellar::stellar_schedule(Mass::from_solar_masses(mass).unwrap(), metallicity).unwrap();
        let mut previous_gas = 0.0;
        for (elapsed, _) in &schedule {
            objects
                .advance_stars_to(*elapsed, SimDuration::from_years(1).unwrap())
                .unwrap();
            assert!((stellar_and_gas_mass(&objects) - mass).abs() < 1e-10);
            assert!(objects.gas_reservoir_solar_masses() >= previous_gas);
            objects.gas_composition().validate().unwrap();
            previous_gas = objects.gas_reservoir_solar_masses();
        }
        assert!(previous_gas > 0.0);
        assert!(objects.gas_composition().metallicity() > metallicity);
        let mut encoded = String::new();
        objects.encode(&mut encoded);
        let mut restored = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
        let events = restored
            .advance_stars_to(
                schedule.last().unwrap().0,
                SimDuration::from_years(1).unwrap(),
            )
            .unwrap();
        assert!(events.is_empty());
        assert_eq!(restored, objects);
    }
}

#[test]
fn supernova_shock_triggers_extra_cluster_formation_across_checkpoints() {
    let (mut shocked, config) = feedback_fixture(20.0, 0.014);
    let tick = SimDuration::from_years(1).unwrap();
    let collapse = stellar::stellar_schedule(Mass::from_solar_masses(20.0).unwrap(), 0.014)
        .unwrap()
        .into_iter()
        .find(|(_, stage)| *stage == StellarStage::NeutronStar)
        .unwrap()
        .0;
    let progenitor_id = shocked.stars().next().unwrap().id();
    shocked.advance_stars_to(collapse, tick).unwrap();
    let cloud_id = shocked.gas_clouds().next().unwrap().id();
    let pending = shocked
        .gas_cloud(cloud_id)
        .unwrap()
        .pending_shock()
        .unwrap();
    assert_eq!(pending.source_star_id(), progenitor_id);
    assert_eq!(pending.arrived_at(), collapse);
    assert!(pending.strength() > 0.0);

    let mut encoded = String::new();
    shocked.encode(&mut encoded);
    assert_eq!(
        ObjectRegistry::decode(encoded.lines().skip(1)).unwrap(),
        shocked
    );

    let interval = SimDuration::from_years(10_000_000).unwrap();
    let next_formation = SimDuration::from_seconds(
        (collapse.seconds() / interval.seconds() + 1) * interval.seconds(),
    );
    let mut control = shocked.clone();
    control.gas_clouds.get_mut(&cloud_id).unwrap().pending_shock = None;
    let baseline_gas = shocked.gas_reservoir_solar_masses();
    control
        .form_stars_between(42, config, collapse, next_formation, tick)
        .unwrap();
    let events = shocked
        .form_stars_between(42, config, collapse, next_formation, tick)
        .unwrap();
    let births: Vec<_> = events
        .iter()
        .filter(|event| event.interaction == Some(StellarInteractionKind::StarFormation))
        .collect();
    assert!(!births.is_empty());
    assert!(shocked.gas_reservoir_solar_masses() < control.gas_reservoir_solar_masses());
    assert!(shocked.gas_reservoir_solar_masses() < baseline_gas);
    let cluster_id = shocked
        .star_system(shocked.star(births[0].star_id).unwrap().system_id())
        .unwrap()
        .cluster_id()
        .unwrap();
    assert_eq!(
        shocked
            .star_cluster(cluster_id)
            .unwrap()
            .triggering_star_id(),
        Some(progenitor_id)
    );
    assert!(
        shocked
            .gas_cloud(cloud_id)
            .unwrap()
            .pending_shock()
            .is_none()
    );

    let mut split = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
    let midpoint = SimDuration::from_seconds((collapse.seconds() + next_formation.seconds()) / 2);
    split
        .form_stars_between(42, config, collapse, midpoint, tick)
        .unwrap();
    split = {
        let mut state = String::new();
        split.encode(&mut state);
        ObjectRegistry::decode(state.lines().skip(1)).unwrap()
    };
    split
        .form_stars_between(42, config, midpoint, next_formation, tick)
        .unwrap();
    assert_eq!(split, shocked);
}

#[test]
fn supernova_shock_compresses_a_nearby_diffuse_cloud() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 32;
    let mut objects = ObjectRegistry::generate(42, config).unwrap();
    let diffuse_id = objects
        .gas_clouds()
        .find(|cloud| cloud.kind() == GasCloudKind::Diffuse)
        .unwrap()
        .id();
    let source_id = objects
        .gas_clouds()
        .find(|cloud| cloud.kind() == GasCloudKind::Molecular)
        .unwrap()
        .id();
    let source_star_id = objects
        .star_systems()
        .find(|system| {
            objects
                .star_cluster(system.cluster_id().unwrap())
                .unwrap()
                .cloud_id()
                == source_id
        })
        .unwrap()
        .star_id();
    let source_position = objects.star(source_star_id).unwrap().position();
    objects.gas_clouds.get_mut(&diffuse_id).unwrap().position = source_position;
    objects
        .apply_supernova_shock(
            source_star_id,
            SimDuration::from_years(1_000_000).unwrap(),
            10.0,
        )
        .unwrap();
    let diffuse = objects.gas_cloud(diffuse_id).unwrap();
    assert_eq!(diffuse.kind(), GasCloudKind::Molecular);
    assert_eq!(
        diffuse.pending_shock().unwrap().source_star_id(),
        source_star_id
    );
    objects.validate().unwrap();
}

#[test]
fn formation_inherits_all_gas_components_and_accounts_for_disk_mass() {
    let (mut objects, config) = feedback_fixture(0.1, 0.014);
    let composition = crate::GasComposition::new(0.68, 0.28, 0.01, 0.015, 0.005, 0.01).unwrap();
    set_test_gas(&mut objects, 10.0, composition);
    let before = stellar_and_gas_mass(&objects)
        + objects
            .protoplanetary_disks()
            .map(|disk| disk.remaining_mass().solar_masses())
            .sum::<f64>();
    let events = objects
        .form_stars_between(
            42,
            config,
            SimDuration::ZERO,
            SimDuration::from_years(10_000_000).unwrap(),
            SimDuration::from_years(1).unwrap(),
        )
        .unwrap();
    let birth = events
        .iter()
        .find(|event| event.interaction == Some(StellarInteractionKind::StarFormation))
        .unwrap();
    let star = objects.star(birth.star_id).unwrap();
    assert_eq!(star.birth_composition(), objects.gas_composition());
    assert_eq!(star.metallicity(), objects.gas_composition().metallicity());
    // The original protostar returns a small amount of unprocessed material.
    assert!(star.birth_composition().carbon_fraction() > 0.009);
    let after = stellar_and_gas_mass(&objects)
        + objects
            .protoplanetary_disks()
            .map(|disk| disk.remaining_mass().solar_masses())
            .sum::<f64>();
    assert!((after - before).abs() < 1e-12);
    let mut encoded = String::new();
    objects.encode(&mut encoded);
    assert_eq!(
        ObjectRegistry::decode(encoded.lines().skip(1)).unwrap(),
        objects
    );
}

#[test]
fn feedback_restarts_formation_in_empty_gas_and_is_checkpoint_independent() {
    let (initial, config) = feedback_fixture(20.0, 0.014);
    let advance = |objects: &mut ObjectRegistry, from_years, to_years| {
        let from = SimDuration::from_years(from_years).unwrap();
        let to = SimDuration::from_years(to_years).unwrap();
        let tick = SimDuration::from_years(1).unwrap();
        let mut events = objects
            .form_stars_between(42, config, from, to, tick)
            .unwrap();
        events.extend(objects.advance_stars_to(to, tick).unwrap());
        events
    };
    let mut whole = initial.clone();
    let mut whole_events = advance(&mut whole, 0, 100_000_000);
    let mut split = initial;
    let mut split_events = Vec::new();
    let mut from = 0;
    for to in [
        333_333,
        7_000_000,
        10_000_000,
        18_123_456,
        33_000_000,
        100_000_000,
    ] {
        split_events.extend(advance(&mut split, from, to));
        let mut encoded = String::new();
        split.encode(&mut encoded);
        split = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
        from = to;
    }
    assert_eq!(split, whole);
    whole_events.sort_by_key(|event| (event.tick, event.star_id, event.to));
    split_events.sort_by_key(|event| (event.tick, event.star_id, event.to));
    assert_eq!(split_events, whole_events);
    assert!(
        whole
            .stars()
            .any(|star| star.formed_at() > SimDuration::ZERO && star.metallicity() > 0.014)
    );
}

#[test]
fn exhausted_reservoir_skips_empty_formation_intervals() {
    let (mut objects, config) = feedback_fixture(0.1, 0.014);
    let elapsed = SimDuration::from_years(10_000_000_000_000).unwrap();
    objects
        .form_stars_between(
            42,
            config,
            SimDuration::ZERO,
            elapsed,
            SimDuration::from_years(1).unwrap(),
        )
        .unwrap();
    objects
        .advance_stars_to(elapsed, SimDuration::from_years(1).unwrap())
        .unwrap();
    assert_eq!(objects.formation_sequence, 1_000_000);
    assert!((stellar_and_gas_mass(&objects) - 0.1).abs() < 1e-12);
}

#[test]
fn merger_feedback_invalidates_old_schedules_and_is_batch_independent() {
    let mut objects = (0..100)
        .find_map(|seed| {
            let objects = ObjectRegistry::generate(seed, GenerationConfig::DEFAULT).unwrap();
            let has_companion = objects
                .star_systems()
                .any(|system| !system.companion_star_ids().is_empty());
            has_companion.then_some(objects)
        })
        .unwrap();
    let system = objects
        .star_systems()
        .find(|system| !system.companion_star_ids().is_empty())
        .unwrap()
        .clone();
    for (id, mass) in [(system.star_id, 1.0), (system.companion_star_ids[0], 0.9)] {
        let mut star = generation::generate_star_with_birth(
            id,
            system.id,
            SimDuration::ZERO,
            crate::GasComposition::default(),
            mass,
        )
        .unwrap();
        if id != system.star_id {
            star.orbital_elements = Some(test_orbit(0.001));
        }
        objects.stars.insert(id, star);
    }
    let total_before = stellar_and_gas_mass(&objects);
    let mut split = objects.clone();
    let tick = SimDuration::from_years(1).unwrap();
    let mut whole_events = objects
        .advance_stars_to(SimDuration::from_years(20_000_000_000).unwrap(), tick)
        .unwrap();
    let mut split_events = Vec::new();
    for years in [100_000, 500_000, 1_000_000, 1_000_000_000, 20_000_000_000] {
        split_events.extend(
            split
                .advance_stars_to(SimDuration::from_years(years).unwrap(), tick)
                .unwrap(),
        );
    }
    assert!(
        whole_events
            .iter()
            .any(|event| event.interaction == Some(StellarInteractionKind::StellarMerger))
    );
    assert!((stellar_and_gas_mass(&objects) - total_before).abs() < 1e-10);
    assert_eq!(split, objects);
    whole_events.sort_by_key(|event| (event.tick, event.star_id, event.to));
    split_events.sort_by_key(|event| (event.tick, event.star_id, event.to));
    assert_eq!(split_events, whole_events);
}

#[test]
fn matter_inventory_accounts_for_star_formation_and_feedback_over_ten_billion_years() {
    let (mut objects, config) = feedback_fixture(1.0, 0.014);
    let before = objects.matter_inventory();
    assert_eq!(before.gas_solar_masses, 0.0);
    assert_eq!(before.stellar_solar_masses, 1.0);
    let tick = SimDuration::from_years(1).unwrap();
    let end = SimDuration::from_years(12_000_000_000).unwrap();
    let mut previous = before;
    for years in [
        10_000_000,
        100_000_000,
        1_000_000_000,
        6_000_000_000,
        12_000_000_000,
    ] {
        let elapsed = SimDuration::from_years(years).unwrap();
        objects.advance_stars_to(elapsed, tick).unwrap();
        let now = objects.matter_inventory();
        assert!((now.tracked_solar_masses() - before.tracked_solar_masses()).abs() < 1e-10);
        assert!(now.gas_solar_masses >= previous.gas_solar_masses);
        assert!(now.gas_metal_solar_masses >= previous.gas_metal_solar_masses);
        previous = now;
    }
    assert!(previous.gas_solar_masses > 0.0);
    assert!(
        previous.tracked_metal_proxy_solar_masses() > before.tracked_metal_proxy_solar_masses()
    );
    let mut encoded = String::new();
    objects.encode(&mut encoded);
    let mut restored = ObjectRegistry::decode(encoded.lines().skip(1)).unwrap();
    restored.advance_stars_to(end, tick).unwrap();
    assert_eq!(restored.matter_inventory(), previous);

    // In a second closed setup, gas and disk metal mass transfer into a new star.
    let (mut forming, _) = feedback_fixture(0.1, 0.014);
    let composition = forming.gas_composition();
    set_test_gas(&mut forming, 10.0, composition);
    let initial = forming.matter_inventory();
    forming
        .form_stars_between(
            42,
            config,
            SimDuration::ZERO,
            SimDuration::from_years(10_000_000).unwrap(),
            tick,
        )
        .unwrap();
    let formed = forming.matter_inventory();
    assert!((formed.tracked_solar_masses() - initial.tracked_solar_masses()).abs() < 1e-10);
    assert!(
        (formed.tracked_metal_proxy_solar_masses() - initial.tracked_metal_proxy_solar_masses())
            .abs()
            < 1e-10
    );
    assert!(formed.protoplanetary_disk_solar_masses >= initial.protoplanetary_disk_solar_masses);
}

#[test]
fn matter_inventory_reports_aggregate_populations_separately() {
    let mut config = GenerationConfig::DEFAULT;
    config.star_systems = 128;
    let objects = ObjectRegistry::generate(42, config).unwrap();
    let budget = objects.matter_inventory();
    assert!(budget.aggregate_solar_masses > 0.0);
    assert!(budget.asteroid_solar_masses > 0.0);
    assert!(budget.comet_solar_masses > 0.0);
    assert!(budget.tracked_solar_masses() > budget.stellar_solar_masses);
    assert!(budget.gas_metal_solar_masses > 0.0);
}

#[test]
fn benchmark_seed_follows_reference_imf_slopes_within_each_branch() {
    // Kroupa, The Local Stellar Initial Mass Function (2000), arXiv:astro-ph/0011328:
    // dN/dM follows approximately M^-1.3 below 0.5 solar masses and M^-2.3 above it.
    // This checks conditional branch slopes, not the normalization across the break.
    let integral = |lo: f64, hi: f64, alpha: f64| {
        (hi.powf(1.0 - alpha) - lo.powf(1.0 - alpha)) / (1.0 - alpha)
    };
    let mut counts = [0_u32; 4];
    let mut rng = StableRng::new(42);
    for _ in 0..100_000 {
        let mass = sample_stellar_mass(&mut rng, 0.1, 16.0);
        let bin = if mass < 0.25 {
            0
        } else if mass < 0.5 {
            1
        } else if mass < 2.0 {
            2
        } else {
            3
        };
        counts[bin] += 1;
    }
    let low_expected = integral(0.1, 0.25, 1.3) / integral(0.25, 0.5, 1.3);
    let high_expected = integral(0.5, 2.0, 2.3) / integral(2.0, 16.0, 2.3);
    let low_observed = f64::from(counts[0]) / f64::from(counts[1]);
    let high_observed = f64::from(counts[2]) / f64::from(counts[3]);
    assert!(
        ((low_observed - low_expected) / low_expected).abs() < 0.05,
        "low-mass ratio {low_observed} differs from {low_expected}"
    );
    assert!(
        ((high_observed - high_expected) / high_expected).abs() < 0.05,
        "high-mass ratio {high_observed} differs from {high_expected}"
    );
}
