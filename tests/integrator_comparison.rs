#![cfg(any(feature = "std", feature = "libm"))]

use azmath::ode::{gauss_legendre4, gauss_legendre4_work_len};
use azmath::optimize::OptimizeOptions;
use azphysics::mechanics::{HarmonicOscillator, State, mechanics_work_len};

macro_rules! comparison {
    ($module:ident, $t:ty, $tol:expr, $roundoff:expr) => {
        mod $module {
            use super::*;
            type T = $t;
            fn gl(steps: usize, end: T) -> ([T; 2], T) {
                let mut y = [1.0, 0.0];
                let mut work = [0.0; gauss_legendre4_work_len(2)];
                let mut energy_error: T = 0.0;
                gauss_legendre4(
                    |_, z, out| {
                        out[0] = z[1];
                        out[1] = -4.0 * z[0];
                    },
                    |_, _, j| j.copy_from_slice(&[0.0, 1.0, -4.0, 0.0]),
                    0.0,
                    end,
                    &mut y,
                    steps,
                    OptimizeOptions {
                        tol_x: 0.0,
                        ..OptimizeOptions::with_tolerance($tol)
                    },
                    &mut work,
                    |_, z| {
                        energy_error =
                            energy_error.max((z[0] * z[0] + z[1] * z[1] / 4.0 - 1.0).abs());
                    },
                )
                .unwrap();
                (y, energy_error)
            }
            fn verlet(steps: usize, end: T) -> ([T; 2], T) {
                let model = HarmonicOscillator::new(2.0 as T, 8.0, [0.0]).unwrap();
                let mut state = State::new([1.0], [0.0]).unwrap();
                let mut work = [0.0; mechanics_work_len(1)];
                let mut energy_error: T = 0.0;
                model
                    .integrate(
                        &mut state,
                        0.0,
                        end / steps as T,
                        steps,
                        &mut work,
                        |_, z| {
                            energy_error =
                                energy_error.max((model.energy(z).unwrap() / 4.0 - 1.0).abs());
                        },
                    )
                    .unwrap();
                ([state.position()[0], state.velocity()[0]], energy_error)
            }
            fn error(y: [T; 2], end: T) -> T {
                (y[0] - (2.0 * end).cos()).hypot((y[1] + 2.0 * (2.0 * end).sin()) / 2.0)
            }
            #[test]
            fn second_vs_fourth_order_against_same_reference() {
                let v_coarse = error(verlet(4, 2.0).0, 2.0);
                let v_fine = error(verlet(8, 2.0).0, 2.0);
                let g_coarse = error(gl(4, 2.0).0, 2.0);
                let g_fine = error(gl(8, 2.0).0, 2.0);
                assert!((3.0..6.0).contains(&(v_coarse / v_fine)));
                assert!(
                    (13.0..18.0).contains(&(g_coarse / g_fine)),
                    "{g_coarse}/{g_fine}"
                );
                assert!(g_fine < v_fine / 10.0);
            }
            #[test]
            fn long_run_energy_at_stable_fixed_step() {
                let (_, gl_energy) = gl(800, 200.0);
                let (_, verlet_energy) = verlet(800, 200.0);
                // omega*dt=0.5: Verlet's relative energy envelope is (omega*dt)^2/4.
                assert!(verlet_energy <= 0.0625 + $roundoff);
                assert!(verlet_energy > 0.05);
                assert!(gl_energy < $roundoff, "GL energy={gl_energy}");
            }
        }
    };
}
comparison!(f64_cases, f64, 2e-14, 2e-10);
comparison!(f32_cases, f32, 6e-7, 2e-4);
