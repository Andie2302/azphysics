use core::fmt;

/// Fehler bei Modellparametern, Zuständen und Integration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PhysicsError {
    /// Ein Zustand benötigt mindestens eine Koordinate.
    InvalidDimension,
    /// Masse muss positiv und endlich sein.
    InvalidMass,
    /// Federkonstante muss nichtnegativ und endlich sein.
    InvalidStiffness,
    /// Schrittweite muss endlich, ungleich null und halbierbar sein.
    InvalidTimeStep,
    /// Mindestens ein Schritt ist erforderlich.
    InvalidStepCount,
    /// Nichtendlicher Zustand, Parameter oder Rechenwert.
    NotFinite,
    /// Eine benötigte positive Modellskala ist im Zahlentyp nicht darstellbar.
    NumericalRange,
    /// Ungültige Bedingung: gleiche oder unbekannte Teilchen, negative oder nichtendliche
    /// Ruhelänge bzw. Compliance.
    InvalidConstraint,
    /// Inverse Masse muss endlich und nichtnegativ sein (0 = fixiertes Teilchen).
    InvalidInverseMass,
    /// Mindestens eine Solver-Iteration ist erforderlich.
    InvalidIterationCount,
    /// Längen von Teilchen-, Geschwindigkeits-, Massen- oder Arbeitsfeldern passen nicht zusammen.
    BufferMismatch,
    /// Fehler des mathematischen Integrators; letzter vollständiger Schritt bleibt erhalten.
    Integration(azmath::ode::OdeError),
}
impl fmt::Display for PhysicsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimension => f.write_str("mindestens eine Koordinate erforderlich"),
            Self::InvalidMass => f.write_str("Masse muss endlich und positiv sein"),
            Self::InvalidStiffness => {
                f.write_str("Federkonstante muss endlich und nichtnegativ sein")
            }
            Self::InvalidTimeStep => {
                f.write_str("Schrittweite muss endlich, ungleich null und halbierbar sein")
            }
            Self::InvalidStepCount => f.write_str("mindestens ein Zeitschritt erforderlich"),
            Self::NotFinite => f.write_str("nichtendlicher Zustand, Parameter oder Rechenwert"),
            Self::NumericalRange => f.write_str("Modellskala im Zahlentyp nicht darstellbar"),
            Self::InvalidConstraint => f.write_str(
                "Bedingung ungültig: gleiche/unbekannte Teilchen oder negative bzw. nichtendliche Ruhelänge/Compliance",
            ),
            Self::InvalidInverseMass => {
                f.write_str("inverse Masse muss endlich und nichtnegativ sein")
            }
            Self::InvalidIterationCount => f.write_str("mindestens eine Iteration erforderlich"),
            Self::BufferMismatch => f.write_str("Feldlängen passen nicht zusammen"),
            Self::Integration(e) => write!(f, "Integration: {e}"),
        }
    }
}
impl core::error::Error for PhysicsError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Integration(e) => Some(e),
            _ => None,
        }
    }
}
impl From<azmath::ode::OdeError> for PhysicsError {
    fn from(e: azmath::ode::OdeError) -> Self {
        Self::Integration(e)
    }
}
