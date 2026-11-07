#![cfg(all(
    feature = "alloc",
    any(feature = "std", feature = "libm", feature = "micromath")
))]

use azphysics::spatial::{Particle, SpatialError, SpatialHash};

macro_rules! cases {
    ($module:ident, $t:ty, $dim:expr) => {
        mod $module {
            use super::*;
            type T = $t;
            const D: usize = $dim;
            fn point(id: usize, x: T) -> Particle<T, D> {
                Particle {
                    id,
                    position: core::array::from_fn(|axis| if axis == 0 { x } else { 0.0 }),
                }
            }
            fn scene() -> Vec<Particle<T, D>> {
                let mut seed = 0x1234_5678_u64;
                (0..96)
                    .map(|id| Particle {
                        id: 1000 - id * 7,
                        position: core::array::from_fn(|_| {
                            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                            ((seed >> 32) % 2001) as T / 100.0 - 10.0
                        }),
                    })
                    .collect()
            }
            fn inside(a: &[T; D], b: &[T; D], radius: T) -> bool {
                a.iter()
                    .zip(b)
                    .map(|(&x, &y)| (x as f64 - y as f64).powi(2))
                    .sum::<f64>()
                    <= (radius as f64).powi(2)
            }
            #[test]
            fn radius_candidates_cover_brute_force() {
                let points = scene();
                for size in [0.125, 0.7, 3.0] {
                    let mut index = SpatialHash::with_bucket_count(size, 1).unwrap();
                    index.rebuild(&points).unwrap();
                    for center in points.iter().take(32) {
                        for radius in [0.0, 0.15, 1.5, 100.0] {
                            let found = index.query_candidates(center.position, radius).unwrap();
                            assert!(found.windows(2).all(|w| w[0] < w[1]));
                            for expected in &points {
                                if inside(&expected.position, &center.position, radius) {
                                    assert!(
                                        found.binary_search(&expected.id).is_ok(),
                                        "missing {}",
                                        expected.id
                                    );
                                }
                            }
                        }
                    }
                }
            }
            #[test]
            fn pairs_cover_brute_force_and_are_unique() {
                let points = scene();
                let mut index = SpatialHash::with_bucket_count(0.5, 1).unwrap();
                index.rebuild(&points).unwrap();
                for radius in [0.0, 0.5, 2.0, 100.0] {
                    let pairs = index.candidate_pairs(radius).unwrap();
                    assert!(pairs.iter().all(|&(a, b)| a < b));
                    assert!(pairs.windows(2).all(|w| w[0] < w[1]));
                    for (i, a) in points.iter().enumerate() {
                        for b in &points[i + 1..] {
                            if inside(&a.position, &b.position, radius) {
                                let pair = (a.id.min(b.id), a.id.max(b.id));
                                assert!(pairs.binary_search(&pair).is_ok());
                            }
                        }
                    }
                }
            }
            #[test]
            fn deterministic_across_order_buckets_and_updates() {
                let mut points = scene();
                let mut a = SpatialHash::with_bucket_count(0.75, 1).unwrap();
                let mut b = SpatialHash::with_bucket_count(0.75, 256).unwrap();
                a.rebuild(&points).unwrap();
                points.reverse();
                b.rebuild(&points).unwrap();
                assert!(a.bucket_count() > 1);
                assert_eq!(
                    a.candidate_pairs(1.25).unwrap(),
                    b.candidate_pairs(1.25).unwrap()
                );
                let saved = points[0];
                a.remove(saved.id).unwrap();
                a.insert(saved).unwrap();
                assert_eq!(
                    a.candidate_pairs(1.25).unwrap(),
                    b.candidate_pairs(1.25).unwrap()
                );
            }
            #[test]
            fn updates_removal_rebuild_and_clear_have_no_stale_entries() {
                let mut index = SpatialHash::with_bucket_count(1.0, 1).unwrap();
                index.insert(point(5, 0.0)).unwrap();
                index.insert(point(2, 0.1)).unwrap();
                assert_eq!(index.candidate_pairs(0.2).unwrap(), vec![(2, 5)]);
                index.update(2, point(2, 20.0).position).unwrap();
                assert!(index.candidate_pairs(0.2).unwrap().is_empty());
                assert_eq!(
                    index.query_candidates(point(0, 0.0).position, 0.2).unwrap(),
                    vec![5]
                );
                index.update(2, point(2, 20.25).position).unwrap();
                assert_eq!(index.position(2), Some(&point(2, 20.25).position));
                assert_eq!(index.remove(2).unwrap(), point(2, 20.25).position);
                assert_eq!(index.remove(2), Err(SpatialError::MissingId));
                assert_eq!(index.len(), 1);
                index.rebuild(&[point(19, -20.0)]).unwrap();
                assert!(index.position(5).is_none());
                assert_eq!(
                    index.query_candidates(point(0, 0.0).position, 0.2).unwrap(),
                    Vec::<usize>::new()
                );
                index.clear();
                assert!(index.is_empty());
                assert!(index.candidate_pairs(0.0).unwrap().is_empty());
                index.insert(point(19, 0.0)).unwrap();
                assert_eq!(index.len(), 1);
            }
            #[test]
            fn errors_are_atomic() {
                let mut index = SpatialHash::new(1.0).unwrap();
                index.insert(point(7, 0.0)).unwrap();
                assert_eq!(index.insert(point(7, 1.0)), Err(SpatialError::DuplicateId));
                assert_eq!(
                    index.update(8, point(8, 1.0).position),
                    Err(SpatialError::MissingId)
                );
                for bad in [T::NAN, T::INFINITY, T::NEG_INFINITY] {
                    assert_eq!(index.insert(point(8, bad)), Err(SpatialError::NotFinite));
                    assert_eq!(
                        index.update(7, point(7, bad).position),
                        Err(SpatialError::NotFinite)
                    );
                    assert_eq!(
                        index.query_candidates(point(0, bad).position, 0.0),
                        Err(SpatialError::NotFinite)
                    );
                    assert_eq!(
                        index.rebuild(&[point(8, 1.0), point(9, bad)]),
                        Err(SpatialError::NotFinite)
                    );
                }
                assert_eq!(
                    index.rebuild(&[point(8, 1.0), point(8, 2.0)]),
                    Err(SpatialError::DuplicateId)
                );
                for bad in [-1.0, T::NAN, T::INFINITY] {
                    assert_eq!(
                        index.query_candidates(point(0, 0.0).position, bad),
                        Err(SpatialError::InvalidRadius)
                    );
                    assert_eq!(index.candidate_pairs(bad), Err(SpatialError::InvalidRadius));
                }
                assert_eq!(index.len(), 1);
                assert_eq!(index.position(7), Some(&point(7, 0.0).position));
            }
            #[test]
            fn boundaries_negative_cells_coincident_and_signed_zero() {
                let near = T::from_bits((1.0 as T).to_bits() - 1);
                let next = T::from_bits((1.0 as T).to_bits() + 1);
                let mut index = SpatialHash::with_bucket_count(1.0, 1).unwrap();
                let points: Vec<_> = [-next, -1.0, -near, -0.0, 0.0, near, 1.0, next]
                    .into_iter()
                    .enumerate()
                    .map(|(id, x)| point(id, x))
                    .collect();
                index.rebuild(&points).unwrap();
                for a in &points {
                    let found = index.query_candidates(a.position, 0.0).unwrap();
                    assert!(found.contains(&a.id));
                    for b in &points {
                        if inside(&a.position, &b.position, 1.0) {
                            assert!(
                                index
                                    .query_candidates(a.position, 1.0)
                                    .unwrap()
                                    .contains(&b.id)
                            );
                        }
                    }
                }
                assert!(index.candidate_pairs(0.0).unwrap().contains(&(3, 4)));
            }
            #[test]
            fn wide_queries_and_coordinate_limits_are_bounded() {
                let mut index = SpatialHash::new(1.0).unwrap();
                let limit = 1_125_899_906_842_624.0 as T;
                index.insert(point(1, limit)).unwrap();
                index.insert(point(2, -limit)).unwrap();
                assert_eq!(
                    index
                        .query_candidates(point(0, 0.0).position, T::MAX)
                        .unwrap(),
                    vec![1, 2]
                );
                assert_eq!(
                    index.insert(point(3, limit * 2.0)),
                    Err(SpatialError::CoordinateRange)
                );
                assert_eq!(
                    index.update(1, point(1, limit * 2.0).position),
                    Err(SpatialError::CoordinateRange)
                );
                assert!(
                    index
                        .query_candidates(point(0, limit).position, 0.0)
                        .unwrap()
                        .contains(&1)
                );
                let mut small = SpatialHash::new(T::from_bits(1)).unwrap();
                small.insert(point(1, 0.0)).unwrap();
                assert_eq!(
                    small.insert(point(2, 1.0)),
                    Err(SpatialError::CoordinateRange)
                );
            }
            #[test]
            fn hash_collisions_keep_distinct_cells_separate() {
                let mut index = SpatialHash::with_bucket_count(1.0, 1).unwrap();
                index.insert(point(1, 0.0)).unwrap();
                index.insert(point(2, 20.0)).unwrap();
                assert_eq!(index.bucket_count(), 1);
                assert_eq!(index.query_candidates([0.0; D], 0.0).unwrap(), vec![1]);
                assert!(index.candidate_pairs(0.1).unwrap().is_empty());
                index.remove(1).unwrap();
                assert_eq!(
                    index
                        .query_candidates(point(0, 20.0).position, 0.0)
                        .unwrap(),
                    vec![2]
                );
            }
            #[test]
            fn empty_and_singleton() {
                let mut index = SpatialHash::new(0.5).unwrap();
                assert!(index.is_empty());
                assert_eq!(index.cell_size(), 0.5);
                assert!(index.query_candidates([0.0; D], 1.0).unwrap().is_empty());
                assert_eq!(index.remove(0), Err(SpatialError::MissingId));
                index.insert(point(usize::MAX, 0.0)).unwrap();
                assert_eq!(
                    index.query_candidates([0.0; D], 0.0).unwrap(),
                    vec![usize::MAX]
                );
                assert!(index.candidate_pairs(1.0).unwrap().is_empty());
                index.rebuild(&[]).unwrap();
                assert!(index.is_empty());
            }
        }
    };
}
cases!(f32_2d, f32, 2);
cases!(f64_2d, f64, 2);
cases!(f32_3d, f32, 3);
cases!(f64_3d, f64, 3);

#[test]
fn constructor_rejects_invalid_shapes_and_scales() {
    assert_eq!(
        SpatialHash::<f64, 0>::new(1.0).unwrap_err(),
        SpatialError::InvalidDimension
    );
    assert_eq!(
        SpatialHash::<f64, 1>::new(1.0).unwrap_err(),
        SpatialError::InvalidDimension
    );
    assert_eq!(
        SpatialHash::<f64, 4>::new(1.0).unwrap_err(),
        SpatialError::InvalidDimension
    );
    for size in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            SpatialHash::<f64, 2>::new(size).unwrap_err(),
            SpatialError::InvalidCellSize
        );
    }
    for count in [0, 3, 63] {
        assert_eq!(
            SpatialHash::<f64, 2>::with_bucket_count(1.0, count).unwrap_err(),
            SpatialError::InvalidBucketCount
        );
    }
    assert_eq!(
        SpatialHash::<f64, 2>::with_bucket_count(1.0, 1_usize << (usize::BITS - 1)).unwrap_err(),
        SpatialError::CapacityOverflow
    );
}

#[test]
fn sparse_scenes_reduce_candidate_count_at_three_sizes() {
    for n in [64, 256, 1024] {
        let points: Vec<_> = (0..n)
            .map(|id| Particle {
                id,
                position: [(id % 32) as f64 * 0.5, (id / 32) as f64 * 0.5],
            })
            .collect();
        let mut index = SpatialHash::new(0.5).unwrap();
        index.rebuild(&points).unwrap();
        let pairs = index.candidate_pairs(0.51).unwrap();
        let all = n * (n - 1) / 2;
        assert!(pairs.len() < all / 2, "n={n}, candidates={}", pairs.len());
        for (i, a) in points.iter().enumerate() {
            for b in &points[i + 1..] {
                let dx = a.position[0] - b.position[0];
                let dy = a.position[1] - b.position[1];
                if dx * dx + dy * dy <= 0.51 * 0.51 {
                    assert!(pairs.binary_search(&(a.id, b.id)).is_ok());
                }
            }
        }
    }
}
