# azphysics

Physikalische Modelle und Simulationen auf Grundlage von
[azmath](https://github.com/Andie2302/azmath).

**Stand 0.0.1: Projektgerüst.** Es sind noch keine Physikverfahren implementiert.
Die Struktur und Qualitätsregeln orientieren sich an azmath: Rust 2024,
`no_std`, optionales `alloc`, dokumentierte Voraussetzungen und passende
Referenztests. `std` ist standardmäßig aktiviert. Eine azmath-Abhängigkeit und
numerische Backend-Features werden beim ersten Modell ergänzt und geprüft.

Die gemeinsame Planung liegt in [aztodo](https://github.com/Andie2302/aztodo):
[Roadmap](https://github.com/Andie2302/aztodo/blob/master/ROADMAP.md),
[Aufgaben](https://github.com/Andie2302/aztodo/blob/master/BACKLOG.md) und
[Architektur](https://github.com/Andie2302/aztodo/blob/master/ARCHITECTURE.md).

Geplant sind zunächst physikalische Konventionen und analytisch prüfbare
Referenzmodelle. Danach können XPBD, räumliche Nachbarschaftssuche,
Druckprojektion und APIC/FLIP folgen. Lineare Löser und allgemeine Integratoren
werden aus azmath verwendet; Kräfte, Materialmodelle, Einheiten und
physikalische Randbedingungen gehören hierher.

```sh
cargo build
cargo build --no-default-features
cargo build --no-default-features --features alloc --target thumbv7em-none-eabihf
cargo clippy --all-targets --all-features -- -D warnings
cargo doc --no-deps --all-features
```

Das Embedded-Target muss zuvor mit `rustup target add thumbv7em-none-eabihf`
installiert werden. Diese Builds prüfen das Gerüst, noch keine numerischen
Verfahren. Die CI-Konfiguration enthält dieselben Prüfungen und einen
MSRV-Build mit Rust 1.85.

Lizenz: MIT oder Apache-2.0.
