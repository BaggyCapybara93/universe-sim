//! Inspectable mass inventory for a materialized region.
//!
//! Aggregate belt/cloud/disk masses describe populations that include sampled
//! individual bodies, so they are reported separately from the tracked total.
//! Stellar and disk metal masses use inherited birth composition as a proxy;
//! stellar interiors, binary-accreted chemistry, and planetary elemental
//! abundances are not resolved by this model.

use crate::objects::ObjectRegistry;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatterInventory {
    pub gas_solar_masses: f64,
    pub stellar_solar_masses: f64,
    pub protoplanetary_disk_solar_masses: f64,
    pub planet_solar_masses: f64,
    pub moon_solar_masses: f64,
    pub asteroid_solar_masses: f64,
    pub comet_solar_masses: f64,
    /// Descriptive aggregate populations; may overlap sampled individual bodies.
    pub aggregate_solar_masses: f64,
    pub gas_metal_solar_masses: f64,
    pub stellar_birth_metal_proxy_solar_masses: f64,
    pub disk_birth_metal_proxy_solar_masses: f64,
}

impl MatterInventory {
    /// All individually represented mass and active disk mass, excluding
    /// descriptive aggregate populations to avoid double counting.
    pub fn tracked_solar_masses(self) -> f64 {
        self.gas_solar_masses
            + self.stellar_solar_masses
            + self.protoplanetary_disk_solar_masses
            + self.planet_solar_masses
            + self.moon_solar_masses
            + self.asteroid_solar_masses
            + self.comet_solar_masses
    }

    /// Metals in gas plus the birth-composition proxy for stars and disks.
    pub fn tracked_metal_proxy_solar_masses(self) -> f64 {
        self.gas_metal_solar_masses
            + self.stellar_birth_metal_proxy_solar_masses
            + self.disk_birth_metal_proxy_solar_masses
    }
}

impl ObjectRegistry {
    pub fn matter_inventory(&self) -> MatterInventory {
        let star_metal = self
            .stars()
            .map(|star| star.mass().solar_masses() * star.birth_composition().metallicity())
            .sum();
        let disk_metal = self
            .protoplanetary_disks()
            .map(|disk| {
                let star = self
                    .star_system(disk.system_id())
                    .and_then(|system| self.star(system.star_id()))
                    .expect("validated disk has a primary star");
                disk.remaining_mass().solar_masses() * star.birth_composition().metallicity()
            })
            .sum();
        MatterInventory {
            gas_solar_masses: self.gas_reservoir_solar_masses(),
            stellar_solar_masses: self.stars().map(|star| star.mass().solar_masses()).sum(),
            protoplanetary_disk_solar_masses: self
                .protoplanetary_disks()
                .map(|disk| disk.remaining_mass().solar_masses())
                .sum(),
            planet_solar_masses: self.planets().map(|body| body.mass().solar_masses()).sum(),
            moon_solar_masses: self.moons().map(|body| body.mass().solar_masses()).sum(),
            asteroid_solar_masses: self
                .asteroids()
                .map(|body| body.mass().solar_masses())
                .sum(),
            comet_solar_masses: self.comets().map(|body| body.mass().solar_masses()).sum(),
            aggregate_solar_masses: self
                .asteroid_belts()
                .map(|body| body.total_mass().solar_masses())
                .sum::<f64>()
                + self
                    .comet_clouds()
                    .map(|body| body.total_mass().solar_masses())
                    .sum::<f64>()
                + self
                    .debris_disks()
                    .map(|body| body.total_mass().solar_masses())
                    .sum::<f64>(),
            gas_metal_solar_masses: self.gas_reservoir_solar_masses()
                * self.gas_composition().metallicity(),
            stellar_birth_metal_proxy_solar_masses: star_metal,
            disk_birth_metal_proxy_solar_masses: disk_metal,
        }
    }
}
