#[cfg(any(feature = "std", feature = "libm"))]
fn main() -> Result<(), azphysics::PhysicsError> {
    use azphysics::mechanics::{HarmonicOscillator, State, mechanics_work_len};
    let model = HarmonicOscillator::new(2.0f64, 8.0, [0.0])?;
    let mut state = State::new([1.0], [0.0])?;
    let mut work = [0.0; mechanics_work_len(1)];
    let initial_energy = model.energy(&state)?;
    model.integrate(&mut state, 0.0, 0.01, 1000, &mut work, |_, _| {})?;
    println!(
        "t=10 s, q={} m, v={} m/s, E={} J (initial {} J)",
        state.position()[0],
        state.velocity()[0],
        model.energy(&state)?,
        initial_energy
    );
    Ok(())
}
#[cfg(not(any(feature = "std", feature = "libm")))]
fn main() {}
