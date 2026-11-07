//! Vorbereitetes, nicht ausgeführtes Beispiel.
#[cfg(any(feature = "std", feature = "libm"))]
fn main() -> Result<(), azphysics::PhysicsError> {
    use azphysics::xpbd::{DistanceConstraint, xpbd_step};
    let mut positions = [[0.0_f64, 0.0], [0.0, -1.0], [0.0, -2.0], [0.0, -3.0]];
    let mut velocities = [[0.0; 2]; 4];
    let inverse_masses = [0.0, 1.0, 1.0, 1.0];
    let constraints = [
        DistanceConstraint::new(0, 1, 1.0, 1e-5)?,
        DistanceConstraint::new(1, 2, 1.0, 1e-5)?,
        DistanceConstraint::new(2, 3, 1.0, 1e-5)?,
    ];
    let mut predicted = positions;
    let mut lambdas = [0.0; 3];
    let mut max_error = 0.0_f64;
    for _ in 0..500 {
        let report = xpbd_step(
            &mut positions, &mut velocities, &inverse_masses, &constraints,
            [0.0, -9.81], 0.01, 32, &mut predicted, &mut lambdas,
        )?;
        max_error = max_error.max(report.max_constraint_error);
    }
    println!("XPBD chain: dt=0.01 s, steps=500, passes=32, compliance=1e-5 s²/kg");
    println!("maximum geometric constraint error: {max_error} m; tip: {:?}", positions[3]);
    println!("Example only; no validation or energy-conservation claim.");
    Ok(())
}
#[cfg(not(any(feature = "std", feature = "libm")))]
fn main() {}
