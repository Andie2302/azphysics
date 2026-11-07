# azphysics

Physikalische Modelle und Simulationen auf Grundlage von
[azmath](https://github.com/Andie2302/azmath).

**Stand 0.0.1: experimentelle Punktmassenmodelle.** Freier Flug, konstante
Beschleunigung und harmonischer Oszillator verwenden Verlet aus azmath.
Die mathematische Grundlage wird bereits direkt von GitHub eingebunden;
es werden keine Quellen kopiert. Rust 2024, `no_std`, optionales `alloc` und
Referenztests folgen dem Aufbau von azmath.

## Festgelegter azmath-Stand

```toml
[dependencies]
azmath = { git = "https://github.com/Andie2302/azmath", rev = "75c5dc04c51ddf02568f15f4ab3c46f7ba4c2c7c", default-features = false }
```

Geprüfte Revision: [`75c5dc04`](https://github.com/Andie2302/azmath/commit/75c5dc04c51ddf02568f15f4ab3c46f7ba4c2c7c).
Die zugehörige Gauss–Legendre-PR (#7) wurde nach 18 erfolgreichen CI-Prüfungen gemergt. `Cargo.lock` hält
zusätzlich den aufgelösten Git-Commit und die weiteren Abhängigkeiten fest.
Ein neuer Commit auf `master` ändert diese Abhängigkeit **nicht automatisch**.
Auch `cargo update` zieht bei festem `rev` keinen anderen azmath-Commit heran.
Eine Branch-Umbenennung beeinflusst diese Commit-Bindung ebenfalls nicht.

Ein bewusstes Update besteht aus drei Schritten:

1. Neue geprüfte SHA in `Cargo.toml` unter `rev` eintragen.
2. Mit `cargo update -p azmath` das Lockfile aktualisieren und die CI-Prüfmatrix ausführen.
3. `Cargo.toml` und `Cargo.lock` gemeinsam committen; Revision und Nachweise in
   [aztodo/STATUS.md](https://github.com/Andie2302/aztodo/blob/master/STATUS.md) aktualisieren.

Ohne festes `rev` könnte Cargo beim Aktualisieren einem Branch folgen; auch
dort verhindert ein bestehendes Lockfile normalerweise Änderungen bei jedem
Build. Hier verwenden wir ausdrücklich die feste Revision.

## Erste Modelle

`mechanics::State` speichert Position in m und Geschwindigkeit in m/s.
`FreeParticle`, `ConstantAcceleration` und `HarmonicOscillator` ergänzen Masse,
Beschleunigungsfeld bzw. Federparameter. Getter, Energie in J und Impuls in
kg*m/s sowie Fehlerfälle sind dokumentiert. Details:
[Mechanik-Konventionen](docs/CONVENTIONS.md) und
[CI-/Referenznachweise](docs/VALIDATION.md).

```rust
use azphysics::mechanics::{State, HarmonicOscillator, mechanics_work_len};
fn main() -> Result<(), azphysics::PhysicsError> {
    let model = HarmonicOscillator::new(2.0f64, 8.0, [0.0])?;
    let mut state = State::new([1.0], [0.0])?;
    let mut work = [0.0; mechanics_work_len(1)];
    model.integrate(&mut state, 0.0, 0.01, 1000, &mut work, |_, _| {})?;
    Ok(())
}
```

Ausführbare Beispiele: `cargo run --locked --example 01_oscillator` und
`cargo run --locked --example 02_integrator_comparison`.

Der zweite Lauf vergleicht Verlet (Ordnung 2) und implizites Gauss–Legendre-RK
(Ordnung 4) am selben Oszillator. Er zeigt dimensionslos skalierte Zustandsfehler,
maximale relative Energiefehler sowie Funktions-/Jacobi-Auswertungen. Gleiche
Schrittzahlen bedeuten keinen gleichen Rechenaufwand; eine Laufzeitmessung
ist dieser Vergleich nicht. Gauss–Legendre erhält die quadratische Energie bei
exakt gelösten Stufen; Newton-Toleranzen und Rundung begrenzen dies praktisch.
Negative Schrittweiten dienen der Rückwärtsintegration. Beim Oszillator muss
für Stabilität `|dt*omega|<2` gelten; kleine Energiefehler erfordern kleinere
Schritte. Die Modelle erzwingen keine automatische Schrittweitenanpassung.

## Räumliche Partikelsuche (optional alloc)

`spatial::SpatialHash<T, D>` besitzt 2D-/3D-Positionskopien mit stabilen IDs.
Einfügen, Update, Entfernung, atomarer Neuaufbau und deterministische Radius-/
Paarkandidaten sind verfügbar. Kandidaten sind konservative Übermengen; eine
exakte Distanz-/Kontaktprüfung folgt separat. Hash-Kollisionen verlieren keine Punkte.

Beispiel: `cargo run --locked --features alloc --example 03_spatial_hash`.
[Vertrag, Paketentscheidung, Tests und Grenzen](docs/SPATIAL.md).
Der Index benötigt `alloc`, bleibt aber mit libm/micromath `no_std`-fähig.

## XPBD-Distanzbedingungen (implementiert, ungeprüft)

`xpbd::DistanceConstraint`, `xpbd_step` und `max_constraint_error` bilden einen allokationsfreien
XPBD-Schritt mit inversen Massen, Fixpunkten und Compliance/dt². Lambda wird
innerhalb eines Schritts akkumuliert und vor jedem neuen Schritt zurückgesetzt.
Positionen, Geschwindigkeiten und inverse Massen werden als Slices übergeben.
Arbeitsbuffer kommen vom Aufrufer; Fehler erhalten Positionen/Geschwindigkeiten.

[Vertrag, Grenzen und offene Abnahme](docs/XPBD.md). Das Kettenbeispiel
`examples/04_xpbd_chain.rs` ist vorbereitet, wurde aber nicht ausgeführt.
**Keine Tests, Builds oder Fuzzing für dieses Paket ausgeführt.** Die Workflows
sind nur noch manuell auslösbar; automatische Push-/PR-CI ist deaktiviert.

## Features

| azphysics | Weitergereicht an azmath | Verwendung |
|---|---|---|
| `std` (Standard) | `std` | Native Mathematikfunktionen |
| `libm` | `libm` | Genaue no_std-Mathematik |
| `micromath` | `micromath` | Näherungen, nur lockere Smoke-Prüfungen |
| `alloc` | `alloc` | Dynamische Matrizen/Sparse-Typen, unabhängig vom Backend |

Bei mehreren Backends gilt die Reihenfolge von azmath: std, dann libm, dann
micromath. Ohne Backend bleibt das Gerüst baubar; für numerische f32/f64-
Anwendungen eines der drei Backends wählen. `azphysics::Real` exportiert das
Trait aus genau dieser azmath-Abhängigkeit.

```sh
cargo test --locked
cargo test --locked --no-default-features --features libm,alloc
cargo test --locked --no-default-features --features micromath,alloc
cargo build --locked --no-default-features --features libm,alloc --target thumbv7em-none-eabihf
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo doc --locked --no-deps --all-features
```

Das Embedded-Target benötigt `rustup target add thumbv7em-none-eabihf`.
Die nur manuell auslösbare CI-Konfiguration enthält Prüfungen für Backend-/alloc-Kombinationen, Builds ohne Backend, Embedded
mit libm und libm+alloc sowie Rust 1.85 mit std und libm+alloc. Die zwei
Integrationstests (mit alloc drei) prüfen die Git-Abhängigkeit und Feature-
Weiterleitung. Die Mechaniktests ergänzen physikalische Referenzlösungen und
Invarianten für std/libm.

## Gemeinsame Planung

[aztodo](https://github.com/Andie2302/aztodo) enthält
[Roadmap](https://github.com/Andie2302/aztodo/blob/master/ROADMAP.md),
[Aufgaben](https://github.com/Andie2302/aztodo/blob/master/BACKLOG.md) und
[Architektur](https://github.com/Andie2302/aztodo/blob/master/ARCHITECTURE.md).

SI-Konventionen und erste analytisch prüfbare Modelle sind implementiert.
Danach können XPBD, räumliche Nachbarschaftssuche,
Druckprojektion und APIC/FLIP folgen. Lineare Löser und allgemeine Integratoren
kommen aus azmath; Kräfte, Materialmodelle, Einheiten und physikalische
Randbedingungen gehören hierher.

Lizenz: MIT oder Apache-2.0.
