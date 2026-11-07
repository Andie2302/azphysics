//! Gleiche Schrittzahl bedeutet nicht gleichen Rechenaufwand: GL4 löst implizite Stufen.
#[cfg(any(feature = "std", feature = "libm"))]
fn main() -> Result<(), azphysics::PhysicsError> {
    use azmath::ode::{gauss_legendre4, gauss_legendre4_work_len};
    use azmath::optimize::OptimizeOptions;
    use azphysics::mechanics::{HarmonicOscillator, State, mechanics_work_len};
    let model = HarmonicOscillator::new(2.0_f64, 8.0, [0.0])?;
    let end = 10.0;
    let exact = [(2.0_f64 * end).cos(), -2.0 * (2.0_f64 * end).sin()];
    println!("m=2 kg, k=8 N/m, q(0)=1 m, v(0)=0 m/s, t={end} s");
    println!("Fehler = max(|dq|/1 m, |dv|/(2 m/s)); Energiefehler relativ zu E(0).");
    println!("steps  method       state error     max energy error    RHS / Jacobian calls");
    for steps in [100, 200, 400] {
        let mut state = State::new([1.0], [0.0])?;
        let mut work = [0.0; mechanics_work_len(1)];
        let mut max_energy: f64 = 0.0;
        model.integrate(
            &mut state,
            0.0,
            end / steps as f64,
            steps,
            &mut work,
            |_, z| {
                max_energy = max_energy.max((model.energy(z).unwrap() / 4.0 - 1.0).abs());
            },
        )?;
        let error = (state.position()[0] - exact[0])
            .abs()
            .max((state.velocity()[0] - exact[1]).abs() / 2.0);
        println!(
            "{steps:5}  Verlet    {error:14.6e}  {max_energy:18.6e}    {} / 0",
            steps + 1
        );

        let mut y = [1.0, 0.0];
        let mut gl_work = [0.0; gauss_legendre4_work_len(2)];
        let mut rhs_calls = 0;
        let mut jac_calls = 0;
        max_energy = 0.0;
        gauss_legendre4(
            |_, z, out| {
                rhs_calls += 1;
                out[0] = z[1];
                out[1] = -4.0 * z[0];
            },
            |_, _, j| {
                jac_calls += 1;
                j.copy_from_slice(&[0.0, 1.0, -4.0, 0.0]);
            },
            0.0,
            end,
            &mut y,
            steps,
            OptimizeOptions::with_tolerance(2e-14),
            &mut gl_work,
            |_, z| {
                max_energy = max_energy.max((z[0] * z[0] + z[1] * z[1] / 4.0 - 1.0).abs());
            },
        )?;
        let error = (y[0] - exact[0]).abs().max((y[1] - exact[1]).abs() / 2.0);
        println!(
            "{steps:5}  GL4       {error:14.6e}  {max_energy:18.6e}    {rhs_calls} / {jac_calls}"
        );
    }
    Ok(())
}
#[cfg(not(any(feature = "std", feature = "libm")))]
fn main() {}
