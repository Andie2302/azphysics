#![cfg(any(feature = "std", feature = "libm", feature = "micromath"))]
use azphysics::Real;

// Exercise the trait exported by azphysics using its selected backend.
fn length<T: Real>(x: T, y: T) -> T {
    (x * x + y * y).sqrt()
}

#[test]
fn selected_backend_supports_generic_math() {
    let result = length(3.0f32, 4.0);
    // Micromath is deliberately approximate: smoke bound at this point, not an accuracy guarantee.
    let tolerance = if cfg!(any(feature = "std", feature = "libm")) {
        1e-6
    } else {
        0.25
    };
    assert!(
        result.is_finite() && (result - 5.0).abs() < tolerance,
        "sqrt backend result: {result}"
    );
    #[cfg(any(feature = "std", feature = "libm"))]
    assert!((length(3.0f64, 4.0) - 5.0).abs() < 1e-14);
}

#[test]
fn dependency_solves_a_linear_system() {
    use azmath::linalg::{Matrix, Vector};
    let a = Matrix::new([[4.0f32, 1.0], [1.0, 3.0]]);
    let b = Vector::new([1.0, 2.0]);
    let x = a.solve(&b).unwrap();
    assert!((4.0 * x[0] + x[1] - 1.0).abs() < 1e-6);
    assert!((x[0] + 3.0 * x[1] - 2.0).abs() < 1e-6);
}

#[cfg(feature = "alloc")]
#[test]
fn alloc_feature_exposes_dynamic_math_types() {
    use azmath::linalg::{DynMatrix, DynVector};
    let a = DynMatrix::new(2, 2, vec![2.0f32, 0.0, 0.0, 4.0]).unwrap();
    let x = a.solve(&DynVector::new(vec![2.0, 8.0])).unwrap();
    assert_eq!(x.as_slice(), &[1.0, 2.0]);
}
