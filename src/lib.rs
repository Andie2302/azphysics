//! Physikalische Modelle und Simulationen für das az-Projekt.
//!
//! Punktmassenmodelle für freien Flug, konstante Beschleunigung und harmonische Oszillatoren.
//! Mathematische Werkzeuge werden direkt aus einer festgelegten azmath-Revision
//! bezogen; std/libm/micromath und alloc werden als Features weitergeleitet. Die gemeinsame
//! Planung liegt unter <https://github.com/Andie2302/aztodo>.
#![no_std]
#[cfg(feature = "alloc")]
extern crate alloc;
#[cfg(feature = "std")]
extern crate std;

/// Gemeinsame Zahlengrundlage aus azmath; Backend-Auswahl über die Crate-Features.
pub use azmath::Real;

mod error;
pub mod mechanics;
#[cfg(feature = "alloc")]
pub mod spatial;
pub use error::PhysicsError;
