//! XPBD-Distanzbedingungen (AZP-040): Punktmassen mit inversen Massen und nachgiebigen
//! paarweisen Abstandsbedingungen. Kein `alloc`, keine Kontakte, keine räumliche Suche.
//!
//! ## Verfahren (ein Zeitschritt `dt`)
//!
//! 1. **Vorhersage** für bewegliche Teilchen (`w > 0`): `v' = v + dt·g`, `x* = x + dt·v'`.
//!    Fixierte Teilchen (`w = 0`) bleiben an ihrem Platz.
//! 2. **Iterationen** (Gauss–Seidel, Bedingungen in Listenreihenfolge): Für
//!    `C = |xᵢ − xⱼ| − L` mit `α̃ = Compliance / dt²` gilt
//!    `Δλ = (−C − α̃·λ) / (wᵢ + wⱼ + α̃)`, `λ += Δλ`, `xᵢ += wᵢ·Δλ·n`, `xⱼ −= wⱼ·Δλ·n`
//!    mit `n = (xᵢ − xⱼ)/|xᵢ − xⱼ|`.
//! 3. **Geschwindigkeit** aus der Positionsänderung: `v = (x_neu − x_alt)/dt`.
//!
//! ## Verträge
//!
//! - **Einheiten (SI):** Position in m, Zeit in s, Beschleunigung in m/s², inverse Masse in 1/kg,
//!   Compliance in m/N (1/Steifigkeit; `0` = starr), `λ` in kg·m.
//! - **Lebensdauer von `λ`:** pro Zeitschritt; er wird am Schrittanfang auf 0 gesetzt, ein
//!   Warmstart über Zeitschritte findet **nicht** statt. Die Kraft der Bedingung ist
//!   `λ · ∇C / dt²` mit `∇ᵢC = n`.
//! - **Fixierte Teilchen** (`w = 0`) bewegen sich nie; ihre Geschwindigkeit wird 0.
//! - **Koinzidente Punkte** (Abstand exakt 0): die Richtung ist undefiniert, die Bedingung wird in
//!   dieser Iteration übersprungen (keine Division durch null). Ist `wᵢ + wⱼ + α̃ = 0` (zwei
//!   fixierte Teilchen, starr), wird ebenfalls übersprungen.
//! - **Atomarer Schritt:** Bei jedem Fehler bleiben `positions` und `velocities` unverändert.
//! - **Genauigkeit:** Mit nachgiebiger Bedingung ist ein Schritt für lineare Federn identisch mit
//!   dem impliziten Euler-Verfahren (Ordnung 1 in `dt`); starre Bedingungen werden nur bis zum
//!   Residuum nach `iterations` Iterationen erfüllt (siehe [`XpbdStats::max_constraint_error`]).
use crate::{PhysicsError, Real};

/// Abstandsbedingung `C(x) = |xₐ − x_b| − Ruhelänge` zwischen zwei Teilchen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceConstraint<T> {
    a: usize,
    b: usize,
    rest_length: T,
    compliance: T,
}

impl<T: Real> DistanceConstraint<T> {
    /// Geprüfte Bedingung. `a != b`, Ruhelänge und Compliance endlich und `>= 0`
    /// (Compliance 0 = starr, in m/N).
    pub fn new(a: usize, b: usize, rest_length: T, compliance: T) -> Result<Self, PhysicsError> {
        if a == b
            || !rest_length.is_finite()
            || rest_length < T::ZERO
            || !compliance.is_finite()
            || compliance < T::ZERO
        {
            return Err(PhysicsError::InvalidConstraint);
        }
        Ok(Self {
            a,
            b,
            rest_length,
            compliance,
        })
    }
    /// Index des ersten Teilchens.
    pub fn a(&self) -> usize {
        self.a
    }
    /// Index des zweiten Teilchens.
    pub fn b(&self) -> usize {
        self.b
    }
    /// Ruhelänge in m.
    pub fn rest_length(&self) -> T {
        self.rest_length
    }
    /// Compliance in m/N.
    pub fn compliance(&self) -> T {
        self.compliance
    }
}

/// Kennzahlen eines Zeitschritts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XpbdStats<T> {
    /// Größter Betrag `|C|` aller Bedingungen **nach** der letzten Iteration, in m.
    pub max_constraint_error: T,
    /// Anzahl der ausgeführten Iterationen.
    pub iterations: usize,
}

/// Euklidische Länge von `d` und Einheitsrichtung; `None` bei Länge 0 oder nichtendlichem Ergebnis.
/// Überlaufsicher durch Skalierung mit der größten Komponente.
fn length_and_direction<T: Real, const D: usize>(d: &[T; D]) -> Option<(T, [T; D])> {
    let scale = d.iter().fold(T::ZERO, |m, &v| m.max(v.abs()));
    if scale == T::ZERO || !scale.is_finite() {
        return None;
    }
    let sum = d.iter().fold(T::ZERO, |s, &v| {
        let q = v / scale;
        s + q * q
    });
    let len = scale * sum.sqrt();
    if !len.is_finite() || len == T::ZERO {
        return None;
    }
    Some((len, core::array::from_fn(|i| d[i] / len)))
}

fn finite_all<T: Real, const D: usize>(v: &[[T; D]]) -> bool {
    v.iter().all(|p| p.iter().all(|x| x.is_finite()))
}

/// Größter Betrag `|C|` der Bedingungen für die gegebenen Positionen (koinzidente Punkte:
/// `|C| = Ruhelänge`). Prüft Indizes und Endlichkeit.
pub fn max_constraint_error<T: Real, const D: usize>(
    positions: &[[T; D]],
    constraints: &[DistanceConstraint<T>],
) -> Result<T, PhysicsError> {
    if D == 0 {
        return Err(PhysicsError::InvalidDimension);
    }
    if !finite_all(positions) {
        return Err(PhysicsError::NotFinite);
    }
    let mut worst = T::ZERO;
    for c in constraints {
        if c.a >= positions.len() || c.b >= positions.len() {
            return Err(PhysicsError::InvalidConstraint);
        }
        let d: [T; D] = core::array::from_fn(|i| positions[c.a][i] - positions[c.b][i]);
        let len = length_and_direction(&d).map_or(T::ZERO, |(l, _)| l);
        let e = (len - c.rest_length).abs();
        if !e.is_finite() {
            return Err(PhysicsError::NotFinite);
        }
        worst = worst.max(e);
    }
    Ok(worst)
}

/// Ein XPBD-Zeitschritt mit gleichförmiger äußerer Beschleunigung (zum Beispiel Schwerkraft).
///
/// - `positions`, `velocities`: Zustand (`n` Teilchen), wird bei Erfolg fortgeschrieben.
/// - `inverse_masses`: `1/m` je Teilchen, `>= 0`; 0 fixiert das Teilchen.
/// - `acceleration`: äußere Beschleunigung in m/s² (wirkt nur auf bewegliche Teilchen).
/// - `dt`: Zeitschritt in s, endlich und `> 0`; `iterations >= 1`.
/// - `predicted`: Arbeitsfeld mit `n` Einträgen; `lambda`: Arbeitsfeld mit `constraints.len()`
///   Einträgen. Beide werden überschrieben, ihr Inhalt vor dem Aufruf ist egal.
///
/// Fehler ([`PhysicsError`]): `InvalidDimension` (`D = 0`), `BufferMismatch`, `InvalidTimeStep`,
/// `InvalidIterationCount`, `InvalidInverseMass`, `InvalidConstraint` (Index außerhalb),
/// `NotFinite` (Zustand, Beschleunigung oder Rechenwert). Bei jedem Fehler bleibt der Zustand
/// unverändert.
#[allow(clippy::too_many_arguments)]
pub fn xpbd_step<T: Real, const D: usize>(
    positions: &mut [[T; D]],
    velocities: &mut [[T; D]],
    inverse_masses: &[T],
    constraints: &[DistanceConstraint<T>],
    acceleration: [T; D],
    dt: T,
    iterations: usize,
    predicted: &mut [[T; D]],
    lambda: &mut [T],
) -> Result<XpbdStats<T>, PhysicsError> {
    if D == 0 {
        return Err(PhysicsError::InvalidDimension);
    }
    let n = positions.len();
    if velocities.len() != n
        || inverse_masses.len() != n
        || predicted.len() != n
        || lambda.len() != constraints.len()
    {
        return Err(PhysicsError::BufferMismatch);
    }
    if !dt.is_finite() || dt <= T::ZERO {
        return Err(PhysicsError::InvalidTimeStep);
    }
    if iterations == 0 {
        return Err(PhysicsError::InvalidIterationCount);
    }
    if inverse_masses
        .iter()
        .any(|w| !w.is_finite() || *w < T::ZERO)
    {
        return Err(PhysicsError::InvalidInverseMass);
    }
    if constraints
        .iter()
        .any(|c| c.a >= n || c.b >= n || c.a == c.b)
    {
        return Err(PhysicsError::InvalidConstraint);
    }
    if !finite_all(positions)
        || !finite_all(velocities)
        || !acceleration.iter().all(|x| x.is_finite())
    {
        return Err(PhysicsError::NotFinite);
    }

    // 1. Vorhersage
    for i in 0..n {
        for k in 0..D {
            predicted[i][k] = if inverse_masses[i] > T::ZERO {
                positions[i][k] + dt * (velocities[i][k] + dt * acceleration[k])
            } else {
                positions[i][k]
            };
        }
    }
    // 2. Iterationen
    lambda.fill(T::ZERO);
    let dt2 = dt * dt;
    for _ in 0..iterations {
        for (c, lam) in constraints.iter().zip(lambda.iter_mut()) {
            let (wa, wb) = (inverse_masses[c.a], inverse_masses[c.b]);
            let alpha = c.compliance / dt2;
            let denom = wa + wb + alpha;
            if denom == T::ZERO || !denom.is_finite() {
                continue;
            }
            let d: [T; D] = core::array::from_fn(|k| predicted[c.a][k] - predicted[c.b][k]);
            let Some((len, dir)) = length_and_direction(&d) else {
                continue;
            };
            let dl = (-(len - c.rest_length) - alpha * *lam) / denom;
            *lam += dl;
            for k in 0..D {
                predicted[c.a][k] += wa * dl * dir[k];
                predicted[c.b][k] -= wb * dl * dir[k];
            }
        }
    }
    // 3. Prüfung vor dem Schreiben (atomar), dann Geschwindigkeit aus der Positionsänderung
    if !finite_all(predicted) {
        return Err(PhysicsError::NotFinite);
    }
    for i in 0..n {
        for k in 0..D {
            if !((predicted[i][k] - positions[i][k]) / dt).is_finite() {
                return Err(PhysicsError::NotFinite);
            }
        }
    }
    let stats = XpbdStats {
        max_constraint_error: max_constraint_error(predicted, constraints)?,
        iterations,
    };
    for i in 0..n {
        for k in 0..D {
            velocities[i][k] = (predicted[i][k] - positions[i][k]) / dt;
            positions[i][k] = predicted[i][k];
        }
    }
    Ok(stats)
}
