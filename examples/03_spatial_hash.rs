#[cfg(all(
    feature = "alloc",
    any(feature = "std", feature = "libm", feature = "micromath")
))]
fn main() -> Result<(), azphysics::spatial::SpatialError> {
    use azphysics::spatial::{Particle, SpatialHash};
    for n in [64, 256, 1024] {
        let points: Vec<_> = (0..n)
            .map(|id| Particle {
                id,
                position: [(id % 32) as f64 * 0.5, (id / 32) as f64 * 0.5],
            })
            .collect();
        let mut index = SpatialHash::new(0.5)?;
        index.rebuild(&points)?;
        let pairs = index.candidate_pairs(0.51)?;
        let mut exact = 0;
        for (i, a) in points.iter().enumerate() {
            for b in &points[i + 1..] {
                let dx = a.position[0] - b.position[0];
                let dy = a.position[1] - b.position[1];
                if dx * dx + dy * dy <= 0.51 * 0.51 {
                    assert!(pairs.binary_search(&(a.id, b.id)).is_ok());
                    exact += 1;
                }
            }
        }
        println!(
            "n={n}, buckets={}, all_pairs={}, candidates={}, exact_neighbors={exact}",
            index.bucket_count(),
            n * (n - 1) / 2,
            pairs.len()
        );
    }
    println!("Candidate-count comparison; not a runtime benchmark.");
    Ok(())
}
#[cfg(not(all(
    feature = "alloc",
    any(feature = "std", feature = "libm", feature = "micromath")
)))]
fn main() {}
