use std::error::Error;
use std::fmt;

/// Normalized bulk composition of a gas reservoir by mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasComposition {
    hydrogen: f64,
    helium: f64,
    carbon: f64,
    oxygen: f64,
    iron: f64,
    other_metals: f64,
}

impl GasComposition {
    pub const COMPONENT_COUNT: usize = 6;

    pub fn new(
        hydrogen: f64,
        helium: f64,
        carbon: f64,
        oxygen: f64,
        iron: f64,
        other_metals: f64,
    ) -> Result<Self, CompositionError> {
        let composition = Self {
            hydrogen,
            helium,
            carbon,
            oxygen,
            iron,
            other_metals,
        };
        composition.validate()?;
        Ok(composition)
    }

    /// Creates a solar-pattern mixture with the requested total metal fraction.
    pub fn from_metallicity(metallicity: f64) -> Result<Self, CompositionError> {
        if !metallicity.is_finite() || !(0.0..=0.1).contains(&metallicity) {
            return Err(CompositionError(
                "gas metallicity must be finite and between 0 and 0.1".to_owned(),
            ));
        }
        let helium = 0.247 + 1.5 * metallicity;
        let hydrogen = 1.0 - helium - metallicity;
        Self::new(
            hydrogen,
            helium,
            metallicity * 0.18,
            metallicity * 0.44,
            metallicity * 0.10,
            metallicity * 0.28,
        )
    }

    pub const fn hydrogen_fraction(self) -> f64 {
        self.hydrogen
    }

    pub const fn helium_fraction(self) -> f64 {
        self.helium
    }

    pub const fn carbon_fraction(self) -> f64 {
        self.carbon
    }

    pub const fn oxygen_fraction(self) -> f64 {
        self.oxygen
    }

    pub const fn iron_fraction(self) -> f64 {
        self.iron
    }

    pub const fn other_metals_fraction(self) -> f64 {
        self.other_metals
    }

    pub const fn metallicity(self) -> f64 {
        self.carbon + self.oxygen + self.iron + self.other_metals
    }

    /// Mix two reservoirs by mass, conserving each component's mass.
    pub fn mixed_with(
        self,
        mass: f64,
        other: Self,
        other_mass: f64,
    ) -> Result<Self, CompositionError> {
        let total = mass + other_mass;
        if !mass.is_finite()
            || !other_mass.is_finite()
            || mass < 0.0
            || other_mass < 0.0
            || !total.is_finite()
            || total <= 0.0
        {
            return Err(CompositionError(
                "mixing requires nonnegative masses and a positive finite total".to_owned(),
            ));
        }
        let weight = other_mass / total;
        let mix = |a: f64, b: f64| a * (mass / total) + b * weight;
        Self::new(
            mix(self.hydrogen, other.hydrogen),
            mix(self.helium, other.helium),
            mix(self.carbon, other.carbon),
            mix(self.oxygen, other.oxygen),
            mix(self.iron, other.iron),
            mix(self.other_metals, other.other_metals),
        )
    }

    /// Toy population yield: convert some H to He, then some He to metals.
    /// Fractions are model parameters, not calibrated nucleosynthesis tables.
    pub(crate) fn processed(self, hydrogen_burned: f64, helium_burned: f64) -> Self {
        let helium_added = self.hydrogen * hydrogen_burned;
        let metals_added = (self.helium + helium_added) * helium_burned;
        Self {
            hydrogen: self.hydrogen - helium_added,
            helium: self.helium + helium_added - metals_added,
            carbon: self.carbon + metals_added * 0.18,
            oxygen: self.oxygen + metals_added * 0.44,
            iron: self.iron + metals_added * 0.10,
            other_metals: self.other_metals + metals_added * 0.28,
        }
    }

    pub(crate) fn validate(self) -> Result<(), CompositionError> {
        let fractions = [
            self.hydrogen,
            self.helium,
            self.carbon,
            self.oxygen,
            self.iron,
            self.other_metals,
        ];
        if fractions
            .iter()
            .any(|fraction| !fraction.is_finite() || !(0.0..=1.0).contains(fraction))
        {
            return Err(CompositionError(
                "gas composition fractions must be finite and between 0 and 1".to_owned(),
            ));
        }
        let total: f64 = fractions.iter().sum();
        if (total - 1.0).abs() > 1.0e-12 {
            return Err(CompositionError(format!(
                "gas composition fractions must sum to 1, found {total:.16}"
            )));
        }
        Ok(())
    }
}

impl Default for GasComposition {
    fn default() -> Self {
        Self::from_metallicity(0.014).expect("solar composition is normalized")
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositionError(String);

impl fmt::Display for CompositionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl Error for CompositionError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_pattern_is_normalized_and_preserves_metallicity() {
        let composition = GasComposition::from_metallicity(0.014).unwrap();
        assert!((composition.metallicity() - 0.014).abs() < 1.0e-12);
        assert!(composition.hydrogen_fraction() > composition.helium_fraction());
        assert!(composition.oxygen_fraction() > composition.iron_fraction());
    }

    #[test]
    fn rejects_unnormalized_composition() {
        assert!(GasComposition::new(0.7, 0.2, 0.0, 0.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn mixing_conserves_each_component_and_handles_empty_reservoirs() {
        let first = GasComposition::default();
        let second = first.processed(0.35, 0.20);
        second.validate().unwrap();
        let mixed = first.mixed_with(3.0, second, 7.0).unwrap();
        for fraction in [
            GasComposition::hydrogen_fraction,
            GasComposition::helium_fraction,
            GasComposition::carbon_fraction,
            GasComposition::oxygen_fraction,
            GasComposition::iron_fraction,
            GasComposition::other_metals_fraction,
        ] {
            assert!(
                (fraction(mixed) * 10.0 - fraction(first) * 3.0 - fraction(second) * 7.0).abs()
                    < 1e-12
            );
        }
        assert_eq!(first.mixed_with(0.0, second, 1.0).unwrap(), second);
        assert!(first.mixed_with(-1.0, second, 2.0).is_err());
        assert!(first.mixed_with(f64::NAN, second, 1.0).is_err());
        assert!(first.mixed_with(0.0, second, 0.0).is_err());
    }
}
