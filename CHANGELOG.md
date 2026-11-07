# Changelog

Version 0.0.1, experimentell und noch nicht auf crates.io veröffentlicht.

## Unveröffentlicht

- XPBD-Fortsetzung ungeprüft: numerische Bereichsfehler werden nicht mehr als
  Koinzidenz behandelt, positive Compliance/dt² darf nicht still unterlaufen;
  Vorhersage, Lambda und Korrekturen werden vor der Zustandsübernahme auf
  Endlichkeit geprüft. Vorhandene Slice-API bleibt erhalten.
- Automatische Push-/PR-CI auf ausdrücklichen Projektwunsch deaktiviert;
  Workflows nur manuell auslösbar. In dieser Sitzung keine Tests/Builds/Fuzzing.

- AZP-040: Modul `xpbd` (ohne alloc) mit `DistanceConstraint`, `xpbd_step` und `max_constraint_error`: XPBD-Distanzbedingungen mit inversen Massen, Compliance, fixierten Teilchen und atomaren Schritten; Tests `tests/xpbd.rs`, Beschreibung in `docs/XPBD.md`. Neue `PhysicsError`-Varianten `InvalidConstraint`, `InvalidInverseMass`, `InvalidIterationCount`, `BufferMismatch`.

- GEO-001/AZP-030: optionales alloc-Modul `spatial`, Spatial Hash für 2D/3D
  mit stabilen IDs, dynamischen Updates und atomarem Neuaufbau.
- Deterministische konservative Radius-/Paarkandidaten, begrenzte Zellabfragen
  mit Scan-Fallback und explizite Eingabe-/Bereichsfehler.
- 38 f32/f64-/Brute-Force-Tests sowie ein CI-ausgeführtes Beispiel.


- SI-/Koordinaten-/Zeitkonventionen und geprüfter Position-Geschwindigkeit-State.
- Freie Punktmasse, konstante Beschleunigung und harmonischer Oszillator mit
  Masse, Federkonstante und Gleichgewichtslage.
- Wiederverwendung von azmath-Verlet, Energie-/Impulsdiagnostik und PhysicsError.
- Analytische f32/f64-Tests, Schrittweitenkonvergenz, Zeitumkehr, Energiefehler
  und Fehlerzustände; Beispiel `01_oscillator`.
- Bewusstes azmath-rev-Update auf den Verlet-Stand mit präzisierter Dokumentation.
