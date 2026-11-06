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
azmath = { git = "https://github.com/Andie2302/azmath", rev = "3c6068a5f218aaba940a85a73606c9f61c8ebdc9", default-features = false }
```

Geprüfte Revision: [`3c6068a5`](https://github.com/Andie2302/azmath/commit/3c6068a5f218aaba940a85a73606c9f61c8ebdc9).
Die zugehörige Verlet-Dokumentations-PR wurde nach 18 erfolgreichen CI-Prüfungen gemergt. `Cargo.lock` hält
zusätzlich den aufgelösten Git-Commit und die weiteren Abhängigkeiten fest.
Ein neuer Commit auf `master` ändert diese Abhängigkeit **nicht automatisch**.
Auch `cargo update` zieht bei festem `rev` keinen anderen azmath-Commit heran.
Eine Branch-Umbenennung beeinflusst diese Commit-Bindung ebenfalls nicht.

Ein bewusstes Update besteht aus drei Schritten:

1. Neue geprüfte SHA in `Cargo.toml` unter `rev` eintragen.
2. Mit `cargo check` das Lockfile aktualisieren und die CI-Prüfmatrix ausführen.
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
[Mechanik-Konventionen](docs/CONVENTIONS.md).

```rust
use azphysics::mechanics::{State, HarmonicOscillator, mechanics_work_len};
let model = HarmonicOscillator::new(2.0f64, 8.0, [0.0])?;
let mut state = State::new([1.0], [0.0])?;
let mut work = [0.0; mechanics_work_len(1)];
model.integrate(&mut state, 0.0, 0.01, 1000, &mut work, |_, _| {})?;
```

Ausführbares Beispiel: `cargo run --locked --example 01_oscillator`.
Negative Schrittweiten dienen der Rückwärtsintegration. Beim Oszillator muss
für Stabilität `|dt*omega|<2` gelten; kleine Energiefehler erfordern kleinere
Schritte. Die Modelle erzwingen keine automatische Schrittweitenanpassung.

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
Die CI prüft die Backend-/alloc-Kombinationen, Builds ohne Backend, Embedded
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
