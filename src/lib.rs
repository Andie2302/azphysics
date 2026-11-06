//! Physikalische Modelle und Simulationen für das az-Projekt.
//!
//! Aktueller Stand: Projektgerüst ohne implementierte Physikverfahren.
//! Mathematische Löser werden künftig aus azmath verwendet. Die gemeinsame
//! Planung liegt unter <https://github.com/Andie2302/aztodo>.
#![no_std]
#[cfg(feature = "alloc")]
extern crate alloc;
#[cfg(feature = "std")]
extern crate std;
