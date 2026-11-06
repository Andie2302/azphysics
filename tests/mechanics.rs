#![cfg(any(feature = "std", feature = "libm"))]
use azmath::ode::OdeError;
use azphysics::{
    PhysicsError,
    mechanics::{
        ConstantAcceleration, FreeParticle, HarmonicOscillator, State, mechanics_work_len,
    },
};

macro_rules! checks {
    ($module:ident,$t:ty) => {
        mod $module {
            use super::*;
            #[test]
            fn free_flight_preserves_velocity_momentum_and_energy() {
                let model = FreeParticle::new(2. as $t).unwrap();
                let mut state = State::new([4., 5., 6.], [1., -2., 3.]).unwrap();
                let p = model.momentum(&state).unwrap();
                let e = model.energy(&state).unwrap();
                assert_eq!(p, [2., -4., 6.]);
                assert_eq!(e, 14.);
                let mut work = [<$t>::NAN; mechanics_work_len(3)];
                model
                    .integrate(&mut state, 0., 0.125, 8, &mut work, |_, _| {})
                    .unwrap();
                assert_eq!(state.position(), &[5., 3., 9.]);
                assert_eq!(state.velocity(), &[1., -2., 3.]);
                assert_eq!(model.momentum(&state).unwrap(), p);
                assert_eq!(model.energy(&state).unwrap(), e);
            }
            #[test]
            fn constant_acceleration_matches_position_velocity_impulse_and_energy() {
                let model = ConstantAcceleration::new(2. as $t, [0., 0., -2.]).unwrap();
                let mut state = State::new([1., 2., 3.], [1., -2., 4.]).unwrap();
                assert_eq!(model.mass(), 2.);
                assert_eq!(model.acceleration(), &[0., 0., -2.]);
                let energy = model.energy(&state).unwrap();
                let momentum = model.momentum(&state).unwrap();
                let mut work = [0.; mechanics_work_len(3)];
                model
                    .integrate(&mut state, 0., 0.125, 16, &mut work, |_, _| {})
                    .unwrap();
                assert_eq!(state.position(), &[3., -2., 7.]);
                assert_eq!(state.velocity(), &[1., -2., 0.]);
                assert_eq!(model.momentum(&state).unwrap()[2] - momentum[2], -8.);
                assert_eq!(model.energy(&state).unwrap(), energy);
            }
            #[test]
            fn oscillator_with_mass_and_shifted_equilibrium_matches_analytic_solution() {
                let model = HarmonicOscillator::new(2. as $t, 8., [0.3, -0.2]).unwrap();
                let mut state = State::new([1., 2.], [0.2, -0.4]).unwrap();
                let initial = state;
                assert_eq!(model.mass(), 2.);
                assert_eq!(model.stiffness(), 8.);
                assert_eq!(model.equilibrium(), &[0.3, -0.2]);
                assert!((model.angular_frequency() - 2.).abs() < 16. * <$t>::EPSILON);
                let mut work = [0.; mechanics_work_len(2)];
                model
                    .integrate(&mut state, 0., 0.0005, 2000, &mut work, |_, _| {})
                    .unwrap();
                let tol = 2048. * <$t>::EPSILON + 2e-6;
                for i in 0..2 {
                    let q = model.equilibrium()[i]
                        + (initial.position()[i] - model.equilibrium()[i]) * (2. as $t).cos()
                        + initial.velocity()[i] / 2. * (2. as $t).sin();
                    let v =
                        -2. * (initial.position()[i] - model.equilibrium()[i]) * (2. as $t).sin()
                            + initial.velocity()[i] * (2. as $t).cos();
                    assert!((state.position()[i] - q).abs() < tol);
                    assert!((state.velocity()[i] - v).abs() < tol * 4.);
                }
            }
            fn oscillator_error(steps: usize) -> $t {
                let model = HarmonicOscillator::new(2. as $t, 8., [0.]).unwrap();
                let mut state = State::new([1.], [0.]).unwrap();
                let mut work = [0.; 4];
                model
                    .integrate(
                        &mut state,
                        0.,
                        1. / steps as $t,
                        steps,
                        &mut work,
                        |_, _| {},
                    )
                    .unwrap();
                let q = state.position()[0] - (2. as $t).cos();
                let v = state.velocity()[0] + 2. * (2. as $t).sin();
                (q * q + v * v).sqrt()
            }
            #[test]
            fn step_halving_has_second_order_convergence() {
                let e1 = oscillator_error(25);
                let e2 = oscillator_error(50);
                let e3 = oscillator_error(100);
                assert!(e1 > 0. && e2 > 0. && e3 > 0.);
                assert!((3.4..4.6).contains(&(e1 / e2)), "ratio {}", e1 / e2);
                assert!((3.4..4.6).contains(&(e2 / e3)), "ratio {}", e2 / e3);
            }
            fn energy_deviation(dt: $t, steps: usize) -> $t {
                let model = HarmonicOscillator::new(2. as $t, 8., [0.]).unwrap();
                let mut state = State::new([1.], [0.]).unwrap();
                let mut work = [0.; 4];
                let e0 = model.energy(&state).unwrap();
                let mut max = 0. as $t;
                model
                    .integrate(&mut state, 0., dt, steps, &mut work, |_, s| {
                        max = max.max((model.energy(s).unwrap() - e0).abs());
                    })
                    .unwrap();
                max
            }
            #[test]
            fn stable_steps_have_bounded_energy_error_and_quadratic_step_dependence() {
                let short = energy_deviation(0.01, 1000);
                let long = energy_deviation(0.01, 10000);
                let coarse = energy_deviation(0.02, 500);
                assert!(short > 0.);
                assert!(
                    long <= 2. * short + 1024. * <$t>::EPSILON,
                    "{short}, {long}"
                );
                assert!(long < 0.001 + 1024. * <$t>::EPSILON, "{long}");
                assert!((3.0..5.0).contains(&(coarse / short)), "{}", coarse / short);
            }
            #[test]
            fn reverse_time_recovers_a_nontrivial_trajectory() {
                let model = HarmonicOscillator::new(2. as $t, 8., [0.]).unwrap();
                let mut state = State::new([1.], [0.3]).unwrap();
                let initial = state;
                let mut work = [0.; 4];
                model
                    .integrate(&mut state, 0., 0.01, 100, &mut work, |_, _| {})
                    .unwrap();
                assert!((state.position()[0] - 1.).abs() > 0.1);
                model
                    .integrate(&mut state, 1., -0.01, 100, &mut work, |_, _| {})
                    .unwrap();
                let tol = 1024. * <$t>::EPSILON;
                assert!((state.position()[0] - initial.position()[0]).abs() < tol);
                assert!((state.velocity()[0] - initial.velocity()[0]).abs() < tol);
            }
            #[test]
            fn observers_see_start_and_completed_steps_and_zero_spring_is_free() {
                let model = HarmonicOscillator::new(2. as $t, 0., [10.]).unwrap();
                let mut state = State::new([1.], [2.]).unwrap();
                let mut work = [0.; 4];
                let mut count = 0;
                let mut last = -1.;
                model
                    .integrate(&mut state, 0., 0.125, 8, &mut work, |t, s| {
                        assert!(t > last);
                        last = t;
                        count += 1;
                        assert_eq!(s.velocity(), &[2.]);
                    })
                    .unwrap();
                assert_eq!(count, 9);
                assert_eq!(last, 1.);
                assert_eq!(state.position(), &[3.]);
                assert_eq!(model.angular_frequency(), 0.);
                assert_eq!(model.energy(&state).unwrap(), 4.);
            }
            #[test]
            fn invalid_parameters_states_and_model_scales_are_rejected() {
                for mass in [0. as $t, -1., <$t>::NAN, <$t>::INFINITY, <$t>::NEG_INFINITY] {
                    assert_eq!(FreeParticle::new(mass), Err(PhysicsError::InvalidMass));
                    assert_eq!(
                        ConstantAcceleration::new(mass, [0.]),
                        Err(PhysicsError::InvalidMass)
                    );
                    assert_eq!(
                        HarmonicOscillator::new(mass, 1., [0.]),
                        Err(PhysicsError::InvalidMass)
                    );
                }
                for k in [-1. as $t, <$t>::NAN, <$t>::INFINITY] {
                    assert_eq!(
                        HarmonicOscillator::new(1., k, [0.]),
                        Err(PhysicsError::InvalidStiffness)
                    );
                }
                assert_eq!(
                    State::<$t, 0>::new([], []),
                    Err(PhysicsError::InvalidDimension)
                );
                assert_eq!(
                    ConstantAcceleration::<$t, 0>::new(1., []),
                    Err(PhysicsError::InvalidDimension)
                );
                assert_eq!(
                    HarmonicOscillator::<$t, 0>::new(1., 1., []),
                    Err(PhysicsError::InvalidDimension)
                );
                assert_eq!(State::new([<$t>::NAN], [0.]), Err(PhysicsError::NotFinite));
                assert_eq!(
                    State::new([0.], [<$t>::INFINITY]),
                    Err(PhysicsError::NotFinite)
                );
                assert_eq!(
                    ConstantAcceleration::new(1., [<$t>::NAN]),
                    Err(PhysicsError::NotFinite)
                );
                assert_eq!(
                    HarmonicOscillator::new(1., 1., [<$t>::INFINITY]),
                    Err(PhysicsError::NotFinite)
                );
                assert_eq!(
                    HarmonicOscillator::new(<$t>::MIN_POSITIVE, <$t>::MAX, [0.]),
                    Err(PhysicsError::NumericalRange)
                );
                assert_eq!(
                    HarmonicOscillator::new(<$t>::MAX, <$t>::MIN_POSITIVE, [0.]),
                    Err(PhysicsError::NumericalRange)
                );
            }
            #[test]
            fn bad_clocks_and_buffers_leave_state_unchanged() {
                let model = FreeParticle::new(1. as $t).unwrap();
                let initial = State::new([1.], [2.]).unwrap();
                let mut work = [0.; 4];
                for dt in [0., <$t>::NAN, <$t>::INFINITY, <$t>::from_bits(1)] {
                    let mut state = initial;
                    assert_eq!(
                        model.integrate(&mut state, 0., dt, 1, &mut work, |_, _| panic!(
                            "invalid clock observed"
                        )),
                        Err(PhysicsError::InvalidTimeStep)
                    );
                    assert_eq!(state, initial);
                }
                let mut state = initial;
                assert_eq!(
                    model.integrate(&mut state, 0., 0.1, 0, &mut work, |_, _| {}),
                    Err(PhysicsError::InvalidStepCount)
                );
                assert_eq!(
                    model.integrate(&mut state, <$t>::INFINITY, 0.1, 1, &mut work, |_, _| {}),
                    Err(PhysicsError::NotFinite)
                );
                assert_eq!(
                    model.integrate(&mut state, 0., <$t>::MAX, 2, &mut work, |_, _| {}),
                    Err(PhysicsError::NotFinite)
                );
                assert_eq!(
                    model.integrate(&mut state, 1e30, 0.001, 1, &mut work, |_, _| {}),
                    Err(PhysicsError::Integration(OdeError::StepSizeTooSmall))
                );
                assert_eq!(
                    model.integrate(&mut state, 0., 0.1, 1, &mut work[..3], |_, _| panic!(
                        "short buffer observed"
                    )),
                    Err(PhysicsError::Integration(OdeError::InvalidInput))
                );
                assert_eq!(state, initial);
            }
            #[test]
            fn overflow_retains_last_complete_state_and_diagnostics_report_errors() {
                let model = HarmonicOscillator::new(1. as $t, 1., [0.]).unwrap();
                let mut state = State::new([1.], [0.]).unwrap();
                let mut last = state;
                let mut count = 0;
                let mut work = [0.; 4];
                let error = model.integrate(&mut state, 0., 3., 1000, &mut work, |_, s| {
                    last = *s;
                    count += 1;
                });
                assert_eq!(error, Err(PhysicsError::Integration(OdeError::NotFinite)));
                assert!(count > 2);
                assert_eq!(last, state);
                assert!(state.position()[0].is_finite() && state.velocity()[0].is_finite());
                let huge = State::new([0. as $t], [<$t>::MAX]).unwrap();
                assert_eq!(huge.kinetic_energy(2.), Err(PhysicsError::NotFinite));
                assert_eq!(huge.momentum(2.), Err(PhysicsError::NotFinite));
                assert_eq!(huge.momentum(0.), Err(PhysicsError::InvalidMass));
            }
        }
    };
}
checks!(single, f32);
checks!(double, f64);
