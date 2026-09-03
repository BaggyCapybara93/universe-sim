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
}
