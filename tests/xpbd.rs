//! Tests für XPBD-Distanzbedingungen (AZP-040). Erwartungswerte sind analytisch hergeleitet;
//! die Herleitungen stehen bei den Tests.
#![cfg(any(feature = "std", feature = "libm"))]

use azphysics::{
    PhysicsError,
    xpbd::{DistanceConstraint, XpbdStats, max_constraint_error, xpbd_step},
};

macro_rules! xpbd_tests {
    ($module:ident, $t:ty, $tol:expr) => {
        mod $module {
            use super::*;
            const TOL: $t = $tol;

            fn close(a: $t, b: $t) {
                assert!(
                    (a - b).abs() <= TOL * (1.0 as $t).max(b.abs()),
                    "erwartet {b}, erhalten {a}"
                );
            }

            /// Schritt mit eigenen Arbeitsfeldern.
            fn step<const D: usize>(
                p: &mut [[$t; D]],
                v: &mut [[$t; D]],
                w: &[$t],
                c: &[DistanceConstraint<$t>],
                g: [$t; D],
                dt: $t,
                iterations: usize,
            ) -> Result<XpbdStats<$t>, PhysicsError> {
                let mut predicted = vec![[0.0; D]; p.len()];
                let mut lambda = vec![0.0; c.len()];
                xpbd_step(p, v, w, c, g, dt, iterations, &mut predicted, &mut lambda)
            }

            fn pair(rest: $t, compliance: $t) -> [DistanceConstraint<$t>; 1] {
                [DistanceConstraint::new(0, 1, rest, compliance).unwrap()]
            }

            #[test]
            fn rigid_pair_with_equal_masses_projects_in_one_iteration() {
                // C = 2 − 1 = 1, Δλ = −C/(1+1) = −0.5: beide Teilchen bewegen sich um 0.5 aufeinander zu.
                let mut p = [[0.0, 0.0], [2.0, 0.0]];
                let mut v = [[0.0; 2]; 2];
                let s = step(
                    &mut p,
                    &mut v,
                    &[1.0, 1.0],
                    &pair(1.0, 0.0),
                    [0.0; 2],
                    0.1,
                    1,
                )
                .unwrap();
                close(p[0][0], 0.5);
                close(p[1][0], 1.5);
                close(p[0][1], 0.0);
                close(p[1][1], 0.0);
                close(v[0][0], 5.0);
                close(v[1][0], -5.0);
                close(s.max_constraint_error, 0.0);
                assert_eq!(s.iterations, 1);
            }

            #[test]
            fn displacement_is_inverse_proportional_to_mass() {
                // w = (1, 3): Δλ = −1/4; a: +0.25, b: −0.75, also Abstand 2 − 1 = 1.
                let mut p = [[0.0, 0.0], [2.0, 0.0]];
                let mut v = [[0.0; 2]; 2];
                step(
                    &mut p,
                    &mut v,
                    &[1.0, 3.0],
                    &pair(1.0, 0.0),
                    [0.0; 2],
                    0.1,
                    1,
                )
                .unwrap();
                close(p[0][0], 0.25);
                close(p[1][0], 1.25);
                // Impulserhaltung: Σ v/w = 0
                close(v[0][0] / 1.0 + v[1][0] / 3.0, 0.0);
            }

            #[test]
            fn compliance_leaves_the_expected_residual() {
                // α̃ = c/dt² = 0.02/0.01 = 2, Δλ = −C/(2+α̃) = −0.25: a: 0.25, b: 1.75,
                // Rest |C| = C·α̃/(Σw + α̃) = 0.5.
                let mut p = [[0.0, 0.0], [2.0, 0.0]];
                let mut v = [[0.0; 2]; 2];
                let s = step(
                    &mut p,
                    &mut v,
                    &[1.0, 1.0],
                    &pair(1.0, 0.02),
                    [0.0; 2],
                    0.1,
                    1,
                )
                .unwrap();
                close(p[0][0], 0.25);
                close(p[1][0], 1.75);
                close(s.max_constraint_error, 0.5);

                // Das Residuum wächst streng mit der Compliance; Compliance 0 ist die starre Grenze.
                let mut last = -1.0 as $t;
                for compliance in [0.0, 1e-3, 1e-2, 1e-1] {
                    let mut p = [[0.0, 0.0], [2.0, 0.0]];
                    let mut v = [[0.0; 2]; 2];
                    let s = step(
                        &mut p,
                        &mut v,
                        &[1.0, 1.0],
                        &pair(1.0, compliance),
                        [0.0; 2],
                        0.1,
                        1,
                    )
                    .unwrap();
                    assert!(
                        s.max_constraint_error > last,
                        "{compliance}: {}",
                        s.max_constraint_error
                    );
                    last = s.max_constraint_error;
                }
            }

            #[test]
            fn fixed_particles_never_move_and_have_zero_velocity() {
                for compliance in [0.0, 0.01] {
                    let mut p = [[0.0, 0.0], [2.0, 0.0]];
                    let mut v = [[3.0, 1.0], [-2.0, 4.0]];
                    let s = step(
                        &mut p,
                        &mut v,
                        &[0.0, 0.0],
                        &pair(1.0, compliance),
                        [0.0, -10.0],
                        0.1,
                        5,
                    )
                    .unwrap();
                    assert_eq!(p, [[0.0, 0.0], [2.0, 0.0]]);
                    assert_eq!(v, [[0.0, 0.0], [0.0, 0.0]]);
                    close(s.max_constraint_error, 1.0);
                }
            }

            #[test]
            fn coincident_points_are_skipped_without_division_by_zero() {
                let mut p = [[1.0, 1.0], [1.0, 1.0]];
                let mut v = [[0.0; 2]; 2];
                let s = step(
                    &mut p,
                    &mut v,
                    &[1.0, 1.0],
                    &pair(1.0, 0.0),
                    [0.0; 2],
                    0.1,
                    4,
                )
                .unwrap();
                assert_eq!(p, [[1.0, 1.0], [1.0, 1.0]]);
                assert!(p.iter().flatten().all(|x| x.is_finite()));
                close(s.max_constraint_error, 1.0);
            }

            #[test]
            fn works_in_three_dimensions() {
                let mut p = [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]];
                let mut v = [[0.0; 3]; 2];
                step(
                    &mut p,
                    &mut v,
                    &[1.0, 1.0],
                    &pair(1.0, 0.0),
                    [0.0; 3],
                    0.1,
                    1,
                )
                .unwrap();
                close(p[0][2], 1.0);
                close(p[1][2], 2.0);
                close(max_constraint_error(&p, &pair(1.0, 0.0)).unwrap(), 0.0);
            }

            #[test]
            fn free_fall_matches_the_closed_form_of_the_predictor() {
                // Ohne Bedingungen: y_N = g·dt²·N(N+1)/2, v_N = g·dt·N (exakt, dt = 0.5).
                let mut p = [[0.0, 0.0]];
                let mut v = [[0.0; 2]];
                let none: [DistanceConstraint<$t>; 0] = [];
                for _ in 0..8 {
                    step(&mut p, &mut v, &[1.0], &none, [0.0, -10.0], 0.5, 1).unwrap();
                }
                assert_eq!(p[0], [0.0, -90.0]);
                assert_eq!(v[0], [0.0, -40.0]);
            }

            /// Teilchen (m) an fester Verankerung mit Compliance 1/k: ein XPBD-Schritt entspricht für
            /// diese lineare Feder dem impliziten Euler-Verfahren `u' = (u + h·v)/(1 + h²k/m)`, und
            /// weitere Iterationen ändern das Ergebnis nicht (Δλ = 0 am Fixpunkt).
            fn spring_matches_implicit_euler(iterations: usize) {
                let (k, m, l, h): ($t, $t, $t, $t) = (100.0, 2.0, 1.0, 0.01);
                let c = [DistanceConstraint::new(0, 1, l, 1.0 / k).unwrap()];
                let w = [0.0, 1.0 / m];
                let mut p = [[0.0, 0.0], [l + 0.2, 0.0]];
                let mut v = [[0.0; 2]; 2];
                let (mut u, mut vel) = (0.2 as $t, 0.0 as $t);
                for _ in 0..50 {
                    step(&mut p, &mut v, &w, &c, [0.0, 0.0], h, iterations).unwrap();
                    let u_pred = u + h * vel;
                    let u_new = u_pred / (1.0 + h * h * k / m);
                    vel = (u_new - u) / h;
                    u = u_new;
                    close(p[1][0] - l, u);
                    assert!((v[1][0] - vel).abs() <= TOL * 1000.0 * (1.0 as $t).max(vel.abs()));
                    assert_eq!((p[0], v[0]), ([0.0, 0.0], [0.0, 0.0]));
                }
            }

            #[test]
            fn linear_spring_equals_implicit_euler_for_any_iteration_count() {
                spring_matches_implicit_euler(1);
                spring_matches_implicit_euler(7);
            }

            fn hanging_chain(iterations: usize, dt: $t, steps: usize) -> ([[$t; 2]; 6], $t) {
                let mut p: [[$t; 2]; 6] = core::array::from_fn(|i| [0.0, -(i as $t)]);
                let mut v = [[0.0; 2]; 6];
                let w: [$t; 6] = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
                let c: Vec<_> = (0..5)
                    .map(|i| DistanceConstraint::new(i, i + 1, 1.0, 0.0).unwrap())
                    .collect();
                let mut err = 0.0;
                for _ in 0..steps {
                    err = step(&mut p, &mut v, &w, &c, [0.0, -9.81], dt, iterations)
                        .unwrap()
                        .max_constraint_error;
                    assert!(p.iter().chain(&v).flatten().all(|x| x.is_finite()));
                }
                (p, err)
            }

            #[test]
            fn hanging_chain_stays_finite_symmetric_and_bounded() {
                let (p, err) = hanging_chain(30, 1.0 / 60.0, 300);
                // Reine Vertikalkräfte: die x-Komponenten bleiben exakt 0, das fixierte Ende bleibt.
                assert!(p.iter().all(|q| q[0] == 0.0));
                assert_eq!(p[0], [0.0, 0.0]);
                // Die Kette hängt nach unten; die Dehnung bleibt klein.
                assert!(p[5][1] < -4.9 && p[5][1] > -5.5, "{}", p[5][1]);
                assert!(err < 0.1, "Residuum {err}");
            }

            #[test]
            fn chain_is_reproducible_and_more_iterations_do_not_increase_the_residual() {
                let (p1, e10) = hanging_chain(10, 1.0 / 60.0, 200);
                let (p2, _) = hanging_chain(10, 1.0 / 60.0, 200);
                assert_eq!(p1, p2);
                let (_, e60) = hanging_chain(60, 1.0 / 60.0, 200);
                assert!(
                    e60 <= e10 + TOL,
                    "10 Iterationen: {e10}, 60 Iterationen: {e60}"
                );
            }

            #[test]
            fn invalid_parameters_are_rejected_and_the_state_stays_unchanged() {
                let c = pair(1.0, 0.0);
                let base_p = [[0.0, 0.0], [2.0, 0.0]];
                let base_v = [[0.5, 0.0], [0.0, 0.5]];
                let check =
                    |dt: $t, iterations: usize, w: [$t; 2], g: [$t; 2], expected: PhysicsError| {
                        let (mut p, mut v) = (base_p, base_v);
                        assert_eq!(
                            step(&mut p, &mut v, &w, &c, g, dt, iterations),
                            Err(expected)
                        );
                        assert_eq!((p, v), (base_p, base_v));
                    };
                let ok_w = [1.0, 1.0];
                for dt in [0.0, -0.1, <$t>::NAN, <$t>::INFINITY] {
                    check(dt, 1, ok_w, [0.0; 2], PhysicsError::InvalidTimeStep);
                }
                check(0.1, 0, ok_w, [0.0; 2], PhysicsError::InvalidIterationCount);
                for bad in [-1.0, <$t>::NAN, <$t>::INFINITY] {
                    check(
                        0.1,
                        1,
                        [1.0, bad],
                        [0.0; 2],
                        PhysicsError::InvalidInverseMass,
                    );
                }
                check(0.1, 1, ok_w, [0.0, <$t>::NAN], PhysicsError::NotFinite);
            }

            #[test]
            fn invalid_constraints_and_buffers_are_rejected() {
                assert_eq!(
                    DistanceConstraint::<$t>::new(1, 1, 1.0, 0.0),
                    Err(PhysicsError::InvalidConstraint)
                );
                for (rest, compliance) in [
                    (-1.0, 0.0),
                    (<$t>::NAN, 0.0),
                    (1.0, -1.0),
                    (1.0, <$t>::NAN),
                    (<$t>::INFINITY, 0.0),
                ] {
                    assert_eq!(
                        DistanceConstraint::<$t>::new(0, 1, rest, compliance),
                        Err(PhysicsError::InvalidConstraint)
                    );
                }
                // Index außerhalb des Teilchenfeldes
                let bad = [DistanceConstraint::new(0, 5, 1.0, 0.0).unwrap()];
                let (mut p, mut v) = ([[0.0, 0.0], [1.0, 0.0]], [[0.0; 2]; 2]);
                assert_eq!(
                    step(&mut p, &mut v, &[1.0, 1.0], &bad, [0.0; 2], 0.1, 1),
                    Err(PhysicsError::InvalidConstraint)
                );
                assert_eq!(
                    max_constraint_error(&p, &bad),
                    Err(PhysicsError::InvalidConstraint)
                );
                // Feldlängen
                let c = pair(1.0, 0.0);
                let (mut pred, mut lam) = ([[0.0; 2]; 2], [0.0; 1]);
                assert_eq!(
                    xpbd_step(
                        &mut p,
                        &mut v[..1],
                        &[1.0, 1.0],
                        &c,
                        [0.0; 2],
                        0.1,
                        1,
                        &mut pred,
                        &mut lam
                    ),
                    Err(PhysicsError::BufferMismatch)
                );
                assert_eq!(
                    xpbd_step(
                        &mut p,
                        &mut v,
                        &[1.0],
                        &c,
                        [0.0; 2],
                        0.1,
                        1,
                        &mut pred,
                        &mut lam
                    ),
                    Err(PhysicsError::BufferMismatch)
                );
                assert_eq!(
                    xpbd_step(
                        &mut p,
                        &mut v,
                        &[1.0, 1.0],
                        &c,
                        [0.0; 2],
                        0.1,
                        1,
                        &mut pred[..1],
                        &mut lam
                    ),
                    Err(PhysicsError::BufferMismatch)
                );
                assert_eq!(
                    xpbd_step(
                        &mut p,
                        &mut v,
                        &[1.0, 1.0],
                        &c,
                        [0.0; 2],
                        0.1,
                        1,
                        &mut pred,
                        &mut []
                    ),
                    Err(PhysicsError::BufferMismatch)
                );
                // Dimension 0
                let (mut p0, mut v0, mut pr0): ([[$t; 0]; 1], [[$t; 0]; 1], [[$t; 0]; 1]) =
                    ([[]], [[]], [[]]);
                assert_eq!(
                    xpbd_step(&mut p0, &mut v0, &[1.0], &[], [], 0.1, 1, &mut pr0, &mut []),
                    Err(PhysicsError::InvalidDimension)
                );
            }

            #[test]
            fn non_finite_states_and_overflow_are_errors_without_changing_the_state() {
                let c = pair(1.0, 0.0);
                for (p0, v0) in [
                    ([[0.0, 0.0], [<$t>::NAN, 0.0]], [[0.0; 2]; 2]),
                    ([[0.0, 0.0], [1.0, 0.0]], [[0.0, <$t>::INFINITY], [0.0; 2]]),
                    // Überlauf in der Vorhersage: x + dt·v = x + 2·MAX
                    ([[0.0, 0.0], [1.0, 0.0]], [[0.0; 2], [<$t>::MAX, 0.0]]),
                ] {
                    let (mut p, mut v) = (p0, v0);
                    assert_eq!(
                        step(&mut p, &mut v, &[1.0, 1.0], &c, [0.0; 2], 2.0, 3),
                        Err(PhysicsError::NotFinite)
                    );
                    // NaN != NaN: bitweise vergleichen
                    let same = p
                        .iter()
                        .flatten()
                        .zip(p0.iter().flatten())
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                        && v.iter()
                            .flatten()
                            .zip(v0.iter().flatten())
                            .all(|(a, b)| a.to_bits() == b.to_bits());
                    assert!(same, "Zustand wurde trotz Fehler verändert");
                }
            }
        }
    };
}

xpbd_tests!(single, f32, 1e-5);
xpbd_tests!(double, f64, 1e-12);

/// Implizites Euler ist ein Verfahren erster Ordnung: Halbieren von dt halbiert den Fehler
/// gegenüber der exakten Lösung `u(t) = u₀·cos(ωt)` des Federschwingers.
#[test]
fn spring_error_decreases_with_first_order_in_the_time_step() {
    let (k, m, l, u0, t_end) = (50.0f64, 1.0f64, 1.0f64, 0.1f64, 0.5f64);
    let omega = (k / m).sqrt();
    let error = |steps: usize| {
        let h = t_end / steps as f64;
        let c = [DistanceConstraint::new(0, 1, l, 1.0 / k).unwrap()];
        let (mut p, mut v) = ([[0.0, 0.0], [l + u0, 0.0]], [[0.0; 2]; 2]);
        let mut pred = [[0.0; 2]; 2];
        let mut lam = [0.0; 1];
        for _ in 0..steps {
            xpbd_step(
                &mut p,
                &mut v,
                &[0.0, 1.0 / m],
                &c,
                [0.0; 2],
                h,
                1,
                &mut pred,
                &mut lam,
            )
            .unwrap();
        }
        (p[1][0] - l - u0 * (omega * t_end).cos()).abs()
    };
    let (e1, e2, e3) = (error(200), error(400), error(800));
    assert!(e1 > e2 && e2 > e3, "{e1} {e2} {e3}");
    for ratio in [e1 / e2, e2 / e3] {
        assert!((1.8..2.2).contains(&ratio), "Fehlerverhältnis {ratio}");
    }
}
