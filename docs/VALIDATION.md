# Mechanik-Nachweise (2026-11-07)

Code-Stand: [azphysics PR #4](https://github.com/Andie2302/azphysics/pull/4),
Commit `1856e1ff55b8651d7f664bc5f99ce540f16eab76`, gemergt als
`1528bea62aba8a856b5cb9b39fae507d0d18e878`.

Die Tests wurden gemäß Projekt-Arbeitsweise **nicht lokal ausgeführt**.
Die GitHub-CI war vor dem Merge für Push und PR vollständig grün:
[PR-Lauf](https://github.com/Andie2302/azphysics/actions/runs/42314757082).

## Tests und Build-Prüfungen

- 20 Mechaniktests: zehn für f32 und zehn für f64, mit std/libm.
- Freier Flug: analytische Bahn, Geschwindigkeit, Impuls und kinetische Energie.
- Konstante Beschleunigung: analytische Bahn/Geschwindigkeit, Impulsbilanz und
  Gesamtenergie mit Potential `-m*a·q`.
- Oszillator: Masse, Federkonstante, verschobenes Gleichgewicht, zwei Koordinaten,
  analytische Lösung, Ordnung 2 bei Halbierung der Schrittweite, Zeitumkehr und
  begrenzter Energiefehler im stabilen Bereich.
- Fehlerfälle: unzulässige Massen/Federn/Zustände, null Dimensionen, nicht
  darstellbare Modellskalen/Uhren, zu kleiner Buffer, Rechenüberlauf und
  Erhalt des letzten vollständigen Zustands.
- Je nach alloc zusätzlich zwei oder drei bestehende Git-/Backend-Integrationstests:
  insgesamt **22 ohne alloc**, **23 mit alloc** bei den genauen Backends.

Die 14 CI-Jobs pro Lauf prüfen sieben Laufzeitkonfigurationen, zwei reine
Compile-Konfigurationen ohne Backend, zwei Embedded-Builds (libm mit/ohne alloc),
zwei Rust-1.85-Builds sowie Formatierung/Clippy/Rustdoc. Die präzisen Physiktests
werden bei reinem micromath übersprungen; dort laufen nur die bisherigen
Backend-/Integrations-Smoke-Tests. Build-Erfolg ist keine Genauigkeitsfreigabe.

## Historische mathematische Grundlage (PR #4)

Cargo.toml und Cargo.lock binden azmath an
`3c6068a5f218aaba940a85a73606c9f61c8ebdc9`.
[azmath PR #6](https://github.com/Andie2302/azmath/pull/6) präzisiert die
Verlet-Voraussetzungen und wurde nach 18 erfolgreichen CI-Jobs gemergt.
Die numerische Implementierung entspricht dem zuvor geprüften Verlet-Stand.

Die Mechanikmodelle delegieren an azmath-Verlet. Für f32/f64 gelten die
Konventionen und Grenzen in [CONVENTIONS.md](CONVENTIONS.md): insbesondere
`|dt*omega|<2` beim Oszillator, hinreichend kleine Schritte für geringe Fehler,
keine automatische Schrittweitenanpassung und keine globale Energiegarantie
für beliebige Kräfte oder Zeitabhängigkeiten. Die Referenzfälle sind ein
gezielter Nachweis, keine vollständige Validierung sämtlicher Parameterwerte.

## Gauss–Legendre-Vergleich (AZM-031)

Aktuelle Git-Abhängigkeit: `75c5dc04c51ddf02568f15f4ab3c46f7ba4c2c7c`,
[azmath PR #7](https://github.com/Andie2302/azmath/pull/7). Vor dem Merge waren
alle 18 CI-Jobs erfolgreich; die 25 neuen Tests umfassen f32/f64-Ordnung 4,
Zeitumkehr, quadratische Invarianten, nicht separierbare Hamilton-Systeme,
nichtlineare/nichtautonome Gleichungen und Fehlerpfade.

`tests/integrator_comparison.rs` ergänzt vier Tests (je zwei pro Präzision):
Verlet-Ordnung 2 gegen GL-Ordnung 4 bei derselben analytischen Referenz sowie
maximale Energieabweichung über 800 feste Schritte mit `omega*dt=0.5`.
Die quadratische Oszillatorenergie ist ein spezieller Invariantentest, keine
Garantie exakter Energieerhaltung für allgemeine nichtlineare Hamilton-Funktionen.
Mit den vorhandenen Tests sind es 26 ohne alloc bzw. 27 mit alloc bei std/libm.

`examples/02_integrator_comparison.rs` zeigt für 100/200/400 Schritte bis t=10 s
Zustands- und Energiefehler sowie die Anzahl der rechten Seiten/Jacobi-Matrizen.
Die CI führt dieses Beispiel mit std/libm und allen Features aus.
Die Ausgabe ist kein Geschwindigkeitsbenchmark. Tests wurden nicht lokal
ausgeführt; CI-Nachweis und Merge-Commit werden nach Abschluss in aztodo erfasst.
