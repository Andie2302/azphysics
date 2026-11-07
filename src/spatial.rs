//! Dynamischer Partikelindex (Spatial Hash) für 2D/3D, benötigt `alloc`.
//!
//! Der Index besitzt Positionskopien und stabile, vom Aufrufer vergebene IDs.
//! Abfragen liefern deterministisch sortierte **Kandidaten**, keine exakten
//! Kugelnachbarn oder Kontakte. Ein Narrow-Phase-Test muss anschließend folgen.
//! Alle Koordinaten und Abstände verwenden dieselbe Längeneinheit (in azphysics SI).
use alloc::collections::BTreeMap;
use alloc::vec;
use alloc::vec::Vec;
use core::fmt;

use crate::Real;

const CELL_LIMIT: f64 = 1_125_899_906_842_624.0; // 2^50, leaves room for rounding halo.
const MAX_QUERY_CELLS: usize = 4096;

/// Fehler der räumlichen Indizierung. Fehler verändern den Index nicht.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SpatialError {
    /// Unterstützt werden genau zwei oder drei Dimensionen.
    InvalidDimension,
    /// Zellgröße muss endlich und strikt positiv sein.
    InvalidCellSize,
    /// Hash-Bucketzahl muss eine positive Zweierpotenz sein.
    InvalidBucketCount,
    /// Abfrageradius muss endlich und nichtnegativ sein.
    InvalidRadius,
    /// Nichtendliche Position.
    NotFinite,
    /// `|position / cell_size| > 2^50` oder nicht darstellbarer Quotient.
    CoordinateRange,
    /// Eine ID darf nur einmal vorkommen.
    DuplicateId,
    /// Update/Entfernung einer nicht vorhandenen ID.
    MissingId,
    /// Die Bucketzahl kann nicht weiter verdoppelt werden.
    CapacityOverflow,
}

impl fmt::Display for SpatialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidDimension => "Spatial Hash unterstützt nur 2D und 3D",
            Self::InvalidCellSize => "Zellgröße muss endlich und positiv sein",
            Self::InvalidBucketCount => "Bucketzahl muss eine positive Zweierpotenz sein",
            Self::InvalidRadius => "Radius muss endlich und nichtnegativ sein",
            Self::NotFinite => "nichtendliche Position",
            Self::CoordinateRange => "Position/Zellgröße außerhalb des Indexbereichs",
            Self::DuplicateId => "doppelte Partikel-ID",
            Self::MissingId => "Partikel-ID nicht vorhanden",
            Self::CapacityOverflow => "Bucketzahl nicht darstellbar",
        })
    }
}
impl core::error::Error for SpatialError {}

/// Ein Partikel für Einfügen und atomaren Neuaufbau.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Particle<T, const D: usize> {
    /// Stabile externe ID, unabhängig von Einfügefolge/Bucketposition.
    pub id: usize,
    /// Kartesische Position in derselben Einheit wie die Zellgröße.
    pub position: [T; D],
}

#[derive(Clone, Debug)]
struct Record<T, const D: usize> {
    position: [T; D],
    cell: [i64; D],
}

/// Hash-Index eines dünn besetzten gleichmäßigen Gitters.
///
/// Hash-Kollisionen werden durch vollständige Zellkoordinaten aufgelöst.
/// Die Bucketzahl wächst automatisch bei mehr als zwei Partikeln pro Bucket
/// im Mittel. ID-Zugriff über BTreeMap: O(log n). Abfragen enumerieren höchstens
/// 4096 Zellen, sonst scannen sie die gespeicherten Partikel. Dichte Cluster
/// können daher O(n) Arbeit pro Abfrage und O(n²) Kandidatenpaare verursachen.
/// Ergebnislisten/Neuaufbau benötigen zusätzliche Allokationen; gewöhnliches
/// Rust-Allokationsversagen wird nicht in [`SpatialError`] umgewandelt.
///
/// Zellkoordinaten werden unabhängig von Näherungs-Backends in f64 gebildet.
/// `|position / cell_size| <= 2^50` begrenzt die Rundungsunsicherheit. Abfragen
/// erweitern Zellgrenzen um eine Zelle; breite oder nicht darstellbare Grenzen
/// führen zu einem konservativen Scan. Das verhindert verlorene Kandidaten an
/// Float-/Zellgrenzen, kann aber zusätzliche Kandidaten liefern.
///
/// ```
/// # #[cfg(any(feature = "std", feature = "libm", feature = "micromath"))] {
/// use azphysics::spatial::{Particle, SpatialHash};
/// let mut index = SpatialHash::<f32, 2>::new(0.5).unwrap();
/// index.insert(Particle { id: 42, position: [0.0, 0.0] }).unwrap();
/// index.insert(Particle { id: 7, position: [0.2, 0.0] }).unwrap();
/// assert_eq!(index.candidate_pairs(0.3).unwrap(), vec![(7, 42)]);
/// index.update(7, [10.0, 0.0]).unwrap();
/// assert!(index.candidate_pairs(0.3).unwrap().is_empty());
/// # }
/// ```
#[derive(Clone, Debug)]
pub struct SpatialHash<T, const D: usize> {
    cell_size: T,
    records: BTreeMap<usize, Record<T, D>>,
    buckets: Vec<Vec<usize>>,
}

fn cell_coordinate(value: f64, size: f64) -> Result<i64, SpatialError> {
    let q = value / size;
    if !q.is_finite() || q.abs() > CELL_LIMIT {
        return Err(SpatialError::CoordinateRange);
    }
    let truncated = q as i64;
    Ok(if q < truncated as f64 {
        truncated - 1
    } else {
        truncated
    })
}

fn hash<const D: usize>(cell: &[i64; D], count: usize) -> usize {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for coordinate in cell {
        for byte in coordinate.to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h as usize & (count - 1)
}

impl<T: Real, const D: usize> SpatialHash<T, D> {
    /// Leerer Index mit 64 anfänglichen Hash-Buckets.
    pub fn new(cell_size: T) -> Result<Self, SpatialError> {
        Self::with_bucket_count(cell_size, 64)
    }

    /// Wie [`Self::new`], mit positiver Zweierpotenz als anfänglicher Bucketzahl.
    pub fn with_bucket_count(cell_size: T, buckets: usize) -> Result<Self, SpatialError> {
        if D != 2 && D != 3 {
            return Err(SpatialError::InvalidDimension);
        }
        if !cell_size.to_f64().is_finite() || cell_size <= T::ZERO {
            return Err(SpatialError::InvalidCellSize);
        }
        if !buckets.is_power_of_two() {
            return Err(SpatialError::InvalidBucketCount);
        }
        if buckets > isize::MAX as usize / core::mem::size_of::<Vec<usize>>() {
            return Err(SpatialError::CapacityOverflow);
        }
        Ok(Self {
            cell_size,
            records: BTreeMap::new(),
            buckets: vec![Vec::new(); buckets],
        })
    }

    /// Zahl gespeicherter Partikel.
    pub fn len(&self) -> usize {
        self.records.len()
    }
    /// Ob keine Partikel gespeichert sind.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    /// Feste Zellgröße; Änderungen erfordern einen neuen Index.
    pub fn cell_size(&self) -> T {
        self.cell_size
    }
    /// Aktuelle Anzahl Hash-Buckets, keine Kapazitätszusage.
    pub fn bucket_count(&self) -> usize {
        self.buckets.len()
    }
    /// Gespeicherte Positionskopie einer ID.
    pub fn position(&self, id: usize) -> Option<&[T; D]> {
        self.records.get(&id).map(|r| &r.position)
    }

    fn cell(&self, position: &[T; D]) -> Result<[i64; D], SpatialError> {
        let mut cell = [0; D];
        for (axis, coordinate) in position.iter().enumerate() {
            if !coordinate.to_f64().is_finite() {
                return Err(SpatialError::NotFinite);
            }
            cell[axis] = cell_coordinate(coordinate.to_f64(), self.cell_size.to_f64())?;
        }
        Ok(cell)
    }

    fn grow(&mut self) -> Result<(), SpatialError> {
        if self.records.len() / self.buckets.len() < 2 {
            return Ok(());
        }
        let count = self
            .buckets
            .len()
            .checked_mul(2)
            .ok_or(SpatialError::CapacityOverflow)?;
        if count > isize::MAX as usize / core::mem::size_of::<Vec<usize>>() {
            return Err(SpatialError::CapacityOverflow);
        }
        let mut buckets = vec![Vec::new(); count];
        for (&id, record) in &self.records {
            buckets[hash(&record.cell, count)].push(id);
        }
        self.buckets = buckets;
        Ok(())
    }

    /// Fügt eine neue ID ein. Doppelte IDs oder ungültige Positionen ändern nichts.
    pub fn insert(&mut self, particle: Particle<T, D>) -> Result<(), SpatialError> {
        if self.records.contains_key(&particle.id) {
            return Err(SpatialError::DuplicateId);
        }
        let cell = self.cell(&particle.position)?;
        self.grow()?;
        let count = self.buckets.len();
        self.buckets[hash(&cell, count)].push(particle.id);
        self.records.insert(
            particle.id,
            Record {
                position: particle.position,
                cell,
            },
        );
        Ok(())
    }

    /// Versetzt eine vorhandene ID; Prüfung erfolgt vor jeder Zustandsänderung.
    pub fn update(&mut self, id: usize, position: [T; D]) -> Result<(), SpatialError> {
        let old = self.records.get(&id).ok_or(SpatialError::MissingId)?.cell;
        let cell = self.cell(&position)?;
        if old != cell {
            let count = self.buckets.len();
            self.buckets[hash(&old, count)].retain(|&other| other != id);
            self.buckets[hash(&cell, count)].push(id);
        }
        self.records.insert(id, Record { position, cell });
        Ok(())
    }

    /// Entfernt eine ID und gibt ihre bisherige Position zurück.
    pub fn remove(&mut self, id: usize) -> Result<[T; D], SpatialError> {
        let record = self.records.remove(&id).ok_or(SpatialError::MissingId)?;
        let count = self.buckets.len();
        self.buckets[hash(&record.cell, count)].retain(|&other| other != id);
        Ok(record.position)
    }

    /// Leert den Index, behält Zellgröße und Bucketzahl.
    pub fn clear(&mut self) {
        self.records.clear();
        for bucket in &mut self.buckets {
            bucket.clear();
        }
    }

    /// Vollständiger Neuaufbau. Bei Fehlern bleibt der gesamte alte Index erhalten.
    pub fn rebuild(&mut self, particles: &[Particle<T, D>]) -> Result<(), SpatialError> {
        let mut next = Self::with_bucket_count(self.cell_size, self.buckets.len())?;
        for &particle in particles {
            next.insert(particle)?;
        }
        *self = next;
        Ok(())
    }

    fn radius(radius: T) -> Result<f64, SpatialError> {
        let radius = radius.to_f64();
        if !radius.is_finite() || radius < 0.0 {
            return Err(SpatialError::InvalidRadius);
        }
        Ok(radius)
    }

    /// Aufsteigend sortierte, eindeutige IDs konservativer Radius-Kandidaten.
    ///
    /// Enthält alle gespeicherten Punkte im abgeschlossenen Radius, zusätzlich
    /// Punkte in der umgebenden Zellbox und Rundungshalos. Radius null ist erlaubt.
    /// Große Radien führen zu einem begrenzten Scan statt riesiger Zellschleifen.
    /// `center` unterliegt denselben Bereichsgrenzen wie gespeicherte Positionen.
    pub fn query_candidates(&self, center: [T; D], radius: T) -> Result<Vec<usize>, SpatialError> {
        self.cell(&center)?;
        let radius = Self::radius(radius)?;
        let size = self.cell_size.to_f64();
        let (mut low, mut high) = ([0; D], [0; D]);
        for axis in 0..D {
            let lo = cell_coordinate(center[axis].to_f64() - radius, size);
            let hi = cell_coordinate(center[axis].to_f64() + radius, size);
            match (lo, hi) {
                (Ok(lo), Ok(hi)) => {
                    low[axis] = lo - 1;
                    high[axis] = hi + 1;
                }
                _ => return Ok(self.records.keys().copied().collect()),
            }
        }
        let count = (0..D).try_fold(1_usize, |count, axis| {
            usize::try_from(high[axis] - low[axis] + 1)
                .ok()
                .and_then(|width| count.checked_mul(width))
        });
        let mut result = Vec::new();
        if count.is_none_or(|count| count > MAX_QUERY_CELLS) {
            for (&id, record) in &self.records {
                if (0..D)
                    .all(|axis| record.cell[axis] >= low[axis] && record.cell[axis] <= high[axis])
                {
                    result.push(id);
                }
            }
        } else {
            let mut cell = low;
            loop {
                for &id in &self.buckets[hash(&cell, self.buckets.len())] {
                    if self.records[&id].cell == cell {
                        result.push(id);
                    }
                }
                let mut axis = 0;
                while axis < D && cell[axis] == high[axis] {
                    cell[axis] = low[axis];
                    axis += 1;
                }
                if axis == D {
                    break;
                }
                cell[axis] += 1;
            }
            result.sort_unstable();
        }
        Ok(result)
    }

    /// Kandidatenpaare für einen gemeinsamen Abstandsradius.
    ///
    /// Jedes Paar `(a,b)` hat `a < b`; die Liste ist lexikographisch sortiert
    /// und enthält keine Duplikate/Selbstpaare. Es sind Broad-Phase-Kandidaten,
    /// keine exakten Abstandstreffer. Die Anzahl kann quadratisch wachsen.
    pub fn candidate_pairs(&self, radius: T) -> Result<Vec<(usize, usize)>, SpatialError> {
        Self::radius(radius)?;
        let mut pairs = Vec::new();
        for (&id, record) in &self.records {
            for other in self.query_candidates(record.position, radius)? {
                if id < other {
                    pairs.push((id, other));
                }
            }
        }
        Ok(pairs)
    }
}
