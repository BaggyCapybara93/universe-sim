//! Stable seed derivation and pseudo-random generation.
//!
//! This module deliberately avoids `DefaultHasher`, whose output is not a
//! persistent file-format guarantee. Changing this algorithm is a simulation
//! compatibility change and must use a new generation namespace or version.

const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

pub fn derive_object_seed(universe_seed: u64, object_id: u64, namespace: &str) -> u64 {
    let mut hash = FNV_OFFSET_BASIS;
    hash_bytes(&mut hash, b"cosmos-object-seed-v1");
    hash_bytes(&mut hash, &universe_seed.to_le_bytes());
    hash_bytes(&mut hash, &object_id.to_le_bytes());
    hash_bytes(&mut hash, namespace.as_bytes());
    splitmix64(hash)
}

pub fn derive_region_seed(universe_seed: u64, x: i64, y: i64, z: i64) -> u64 {
    let mut hash = FNV_OFFSET_BASIS;
    hash_bytes(&mut hash, b"cosmos-region-seed-v1");
    hash_bytes(&mut hash, &universe_seed.to_le_bytes());
    hash_bytes(&mut hash, &x.to_le_bytes());
    hash_bytes(&mut hash, &y.to_le_bytes());
    hash_bytes(&mut hash, &z.to_le_bytes());
    splitmix64(hash)
}

fn hash_bytes(hash: &mut u64, bytes: &[u8]) {
    for byte in bytes {
        *hash ^= u64::from(*byte);
        *hash = hash.wrapping_mul(FNV_PRIME);
    }
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StableRng {
    state: u64,
}

impl StableRng {
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        splitmix64(self.state)
    }

    pub fn next_f64(&mut self) -> f64 {
        const DENOMINATOR: f64 = (1_u64 << 53) as f64;
        ((self.next_u64() >> 11) as f64) / DENOMINATOR
    }

    pub fn next_u64_below(&mut self, upper_bound: u64) -> Option<u64> {
        if upper_bound == 0 {
            return None;
        }

        let threshold = upper_bound.wrapping_neg() % upper_bound;
        loop {
            let value = self.next_u64();
            if value >= threshold {
                return Some(value % upper_bound);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn object_seeds_are_reproducible_and_isolated() {
        let first = derive_object_seed(918_273_645, 81, "star");
        assert_eq!(first, derive_object_seed(918_273_645, 81, "star"));
        assert_ne!(first, derive_object_seed(918_273_645, 82, "star"));
        assert_ne!(first, derive_object_seed(918_273_645, 81, "planet"));
    }

    #[test]
    fn object_seed_algorithm_has_a_compatibility_vector() {
        assert_eq!(
            derive_object_seed(918_273_645, 0, "universe"),
            4_808_633_388_830_255_789
        );
    }

    #[test]
    fn random_values_are_in_the_unit_interval() {
        let mut rng = StableRng::new(42);
        for _ in 0..100 {
            assert!((0.0..1.0).contains(&rng.next_f64()));
        }
    }

    #[test]
    fn bounded_values_stay_in_range() {
        let mut rng = StableRng::new(42);
        for _ in 0..100 {
            assert!(rng.next_u64_below(7).unwrap() < 7);
        }
        assert_eq!(rng.next_u64_below(0), None);
    }

    #[test]
    fn region_seeds_are_coordinate_scoped() {
        let origin = derive_region_seed(42, 0, 0, 0);
        assert_eq!(origin, derive_region_seed(42, 0, 0, 0));
        assert_ne!(origin, derive_region_seed(42, 1, 0, 0));
    }
}
