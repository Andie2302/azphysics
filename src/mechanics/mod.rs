//! Punktmassen in kartesischen Koordinaten: SI-Konventionen, keine Laufzeit-Einheitentypen.
//!
//! Position in m, Geschwindigkeit in m/s, Masse in kg, Zeit in s, Energie in J.
//! Die Modelle sind autonom und konservativ. In 3D gilt ein rechtshändiges
//! Koordinatensystem; die Richtung einer Beschleunigung wird explizit angegeben.
//! Alle Zustände und Arbeitsbuffer haben feste Größe, der Kern benötigt kein alloc.
use crate::{PhysicsError, Real};
use azmath::ode::{verlet, verlet_work_len};

/// Position und Geschwindigkeit einer Punktmasse, **kein** Impulszustand.
///
/// `D >= 1`; für räumliche Modelle typischerweise 1, 2 oder 3.
/// Private Felder erhalten die Prüfung auf endliche Anfangsdaten.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State<T, const D: usize> {
    position: [T; D],
    velocity: [T; D],
}
impl<T: Real, const D: usize> State<T, D> {
    /// Geprüfter Anfangszustand in m und m/s.
    pub fn new(position: [T; D], velocity: [T; D]) -> Result<Self, PhysicsError> {
        if D == 0 {
            return Err(PhysicsError::InvalidDimension);
        }
        if !position.iter().chain(&velocity).all(|v| v.is_finite()) {
            return Err(PhysicsError::NotFinite);
        }
        Ok(Self { position, velocity })
    }
    /// Kartesische Position in m.
    pub fn position(&self) -> &[T; D] {
        &self.position
    }
    /// Kartesische Geschwindigkeit in m/s.
    pub fn velocity(&self) -> &[T; D] {
        &self.velocity
    }
    /// Impuls `m*v` in kg*m/s; Masse wird geprüft.
    pub fn momentum(&self, mass: T) -> Result<[T; D], PhysicsError> {
        check_mass(mass)?;
        let p = core::array::from_fn(|i| mass * self.velocity[i]);
        if !p.iter().all(|v: &T| v.is_finite()) {
            return Err(PhysicsError::NotFinite);
        }
        Ok(p)
    }
    /// Kinetische Energie `m*|v|²/2` in J; Masse wird geprüft.
    pub fn kinetic_energy(&self, mass: T) -> Result<T, PhysicsError> {
        check_mass(mass)?;
        let half = T::from_f64(0.5);
        checked(
            self.velocity
                .iter()
                .fold(T::ZERO, |s, &v| s + (mass * (v * half)) * v),
        )
    }
}
fn check_mass<T: Real>(mass: T) -> Result<(), PhysicsError> {
    if !mass.is_finite() || mass <= T::ZERO {
        Err(PhysicsError::InvalidMass)
    } else {
        Ok(())
    }
}
fn checked<T: Real>(value: T) -> Result<T, PhysicsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PhysicsError::NotFinite)
    }
}

// Shared clock/observer adapter. The integration algorithm stays in azmath.
#[allow(clippy::too_many_arguments)]
fn advance<T: Real, const D: usize>(
    accel: impl FnMut(T, &[T], &mut [T]),
    state: &mut State<T, D>,
    t0: T,
    dt: T,
    steps: usize,
    work: &mut [T],
    mut observer: impl FnMut(T, &State<T, D>),
) -> Result<(), PhysicsError> {
    if !t0.is_finite() {
        return Err(PhysicsError::NotFinite);
    }
    if !dt.is_finite() || dt == T::ZERO || dt * T::from_f64(0.5) == T::ZERO {
        return Err(PhysicsError::InvalidTimeStep);
    }
    if steps == 0 {
        return Err(PhysicsError::InvalidStepCount);
    }
    let duration = checked(dt * T::from_usize(steps))?;
    let t1 = checked(t0 + duration)?;
    if (dt > T::ZERO && t1 <= t0) || (dt < T::ZERO && t1 >= t0) {
        return Err(azmath::ode::OdeError::StepSizeTooSmall.into());
    }
    verlet(
        accel,
        t0,
        t1,
        &mut state.position,
        &mut state.velocity,
        steps,
        work,
        |t, q, v| {
            let snapshot = State {
                position: core::array::from_fn(|i| q[i]),
                velocity: core::array::from_fn(|i| v[i]),
            };
            observer(t, &snapshot);
        },
    )
    .map_err(PhysicsError::from)
}

/// Arbeitsspeicher der mechanischen Modelle: der Verlet-Buffer aus azmath (`4*D`).
pub const fn mechanics_work_len(d: usize) -> usize {
    verlet_work_len(d)
}

/// Freie Punktmasse ohne äußere Kräfte.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeParticle<T> {
    mass: T,
}
impl<T: Real> FreeParticle<T> {
    /// Positive endliche Masse in kg.
    pub fn new(mass: T) -> Result<Self, PhysicsError> {
        check_mass(mass)?;
        Ok(Self { mass })
    }
    /// Masse in kg.
    pub fn mass(&self) -> T {
        self.mass
    }
    /// Gesamtenergie ist hier die kinetische Energie.
    pub fn energy<const D: usize>(&self, state: &State<T, D>) -> Result<T, PhysicsError> {
        state.kinetic_energy(self.mass)
    }
    /// Impuls in kg*m/s.
    pub fn momentum<const D: usize>(&self, state: &State<T, D>) -> Result<[T; D], PhysicsError> {
        state.momentum(self.mass)
    }
    /// Integriert `steps` nominelle Schritte `dt` ab `t0` mit azmath-Verlet.
    ///
    /// `dt` darf für Rückwärtsintegration negativ sein. Der Observer sieht den
    /// Anfangszustand und jeden vollständigen Schritt. Zeitraster und Puffer
    /// werden geprüft; bei Integrationsfehlern bleibt der letzte vollständige
    /// Zustand erhalten. Für die Floating-Point-Darstellung des Zeitrasters
    /// gelten die Grenzen von [`azmath::ode::verlet`].
    #[allow(clippy::too_many_arguments)]
    pub fn integrate<const D: usize>(
        &self,
        state: &mut State<T, D>,
        t0: T,
        dt: T,
        steps: usize,
        work: &mut [T],
        observer: impl FnMut(T, &State<T, D>),
    ) -> Result<(), PhysicsError> {
        advance(
            |_, _, a| a.fill(T::ZERO),
            state,
            t0,
            dt,
            steps,
            work,
            observer,
        )
    }
}

/// Punktmasse in einem räumlich und zeitlich konstanten Beschleunigungsfeld.
///
/// Die Beschleunigung in m/s² enthält bereits die Massenteilung. Das Potential
/// wird mit Nullpunkt im Koordinatenursprung als `V=-m*a·q` definiert.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantAcceleration<T, const D: usize> {
    mass: T,
    acceleration: [T; D],
}
impl<T: Real, const D: usize> ConstantAcceleration<T, D> {
    /// Masse in kg und Beschleunigungsvektor in m/s².
    pub fn new(mass: T, acceleration: [T; D]) -> Result<Self, PhysicsError> {
        check_mass(mass)?;
        if D == 0 {
            return Err(PhysicsError::InvalidDimension);
        }
        if !acceleration.iter().all(|v| v.is_finite()) {
            return Err(PhysicsError::NotFinite);
        }
        Ok(Self { mass, acceleration })
    }
    /// Masse in kg.
    pub fn mass(&self) -> T {
        self.mass
    }
    /// Konstante Beschleunigung in m/s².
    pub fn acceleration(&self) -> &[T; D] {
        &self.acceleration
    }
    /// Mechanische Energie `m*|v|²/2 - m*a·q` in J.
    pub fn energy(&self, state: &State<T, D>) -> Result<T, PhysicsError> {
        let potential = (0..D).fold(T::ZERO, |s, i| {
            s - self.mass * self.acceleration[i] * state.position[i]
        });
        checked(state.kinetic_energy(self.mass)? + potential)
    }
    /// Impuls in kg*m/s; unter äußerer Kraft im Allgemeinen nicht konstant.
    pub fn momentum(&self, state: &State<T, D>) -> Result<[T; D], PhysicsError> {
        state.momentum(self.mass)
    }
    /// Wie [`FreeParticle::integrate`], mit dem konstanten Beschleunigungsfeld.
    #[allow(clippy::too_many_arguments)]
    pub fn integrate(
        &self,
        state: &mut State<T, D>,
        t0: T,
        dt: T,
        steps: usize,
        work: &mut [T],
        observer: impl FnMut(T, &State<T, D>),
    ) -> Result<(), PhysicsError> {
        advance(
            |_, _, a| a.copy_from_slice(&self.acceleration),
            state,
            t0,
            dt,
            steps,
            work,
            observer,
        )
    }
}

/// Isotroper harmonischer Oszillator um eine feste Gleichgewichtslage.
///
/// `V=k*|q-q_eq|²/2`, `a=-(k/m)*(q-q_eq)`; jede Koordinate ist entkoppelt.
/// `k=0` erlaubt den freien Grenzfall. Stabiler Verlet-Betrieb benötigt
/// `|dt*omega|<2` und für kleine Fehler deutlich kleinere Schritte. Das Modell
/// ändert die gewählte Schrittweite nicht automatisch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicOscillator<T, const D: usize> {
    mass: T,
    stiffness: T,
    equilibrium: [T; D],
    omega_squared: T,
}
impl<T: Real, const D: usize> HarmonicOscillator<T, D> {
    /// Masse in kg, Federkonstante in N/m und Gleichgewichtslage in m.
    pub fn new(mass: T, stiffness: T, equilibrium: [T; D]) -> Result<Self, PhysicsError> {
        check_mass(mass)?;
        if !stiffness.is_finite() || stiffness < T::ZERO {
            return Err(PhysicsError::InvalidStiffness);
        }
        if D == 0 {
            return Err(PhysicsError::InvalidDimension);
        }
        if !equilibrium.iter().all(|v| v.is_finite()) {
            return Err(PhysicsError::NotFinite);
        }
        let omega_squared = stiffness / mass;
        if !omega_squared.is_finite() || (stiffness > T::ZERO && omega_squared == T::ZERO) {
            return Err(PhysicsError::NumericalRange);
        }
        Ok(Self {
            mass,
            stiffness,
            equilibrium,
            omega_squared,
        })
    }
    /// Masse in kg.
    pub fn mass(&self) -> T {
        self.mass
    }
    /// Federkonstante in N/m.
    pub fn stiffness(&self) -> T {
        self.stiffness
    }
    /// Gleichgewichtslage in m.
    pub fn equilibrium(&self) -> &[T; D] {
        &self.equilibrium
    }
    /// Kreisfrequenz in rad/s, Genauigkeit gemäß dem gewählten Backend.
    pub fn angular_frequency(&self) -> T {
        self.omega_squared.sqrt()
    }
    /// Gesamtenergie in J; Nullpunkt des Potentials liegt im Gleichgewicht.
    pub fn energy(&self, state: &State<T, D>) -> Result<T, PhysicsError> {
        let kinetic = state.kinetic_energy(self.mass)?;
        if self.stiffness == T::ZERO {
            return Ok(kinetic);
        }
        let half = T::from_f64(0.5);
        let potential = (0..D).fold(T::ZERO, |s, i| {
            let q = state.position[i] - self.equilibrium[i];
            s + (self.stiffness * (q * half)) * q
        });
        checked(kinetic + potential)
    }
    /// Impuls in kg*m/s; unter Federkraft im Allgemeinen nicht konstant.
    pub fn momentum(&self, state: &State<T, D>) -> Result<[T; D], PhysicsError> {
        state.momentum(self.mass)
    }
    /// Wie [`FreeParticle::integrate`], mit der Federbeschleunigung.
    #[allow(clippy::too_many_arguments)]
    pub fn integrate(
        &self,
        state: &mut State<T, D>,
        t0: T,
        dt: T,
        steps: usize,
        work: &mut [T],
        observer: impl FnMut(T, &State<T, D>),
    ) -> Result<(), PhysicsError> {
        advance(
            |_, q, a| {
                if self.stiffness == T::ZERO {
                    a.fill(T::ZERO);
                } else {
                    for i in 0..D {
                        a[i] = -self.omega_squared * (q[i] - self.equilibrium[i]);
                    }
                }
            },
            state,
            t0,
            dt,
            steps,
            work,
            observer,
        )
    }
}
