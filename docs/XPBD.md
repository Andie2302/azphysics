# XPBD-Distanzbedingungen (AZP-040)

**Implementiert, ungeprüft.** In dieser Sitzung wurden keine Tests, Builds,
Prüfwerkzeuge, Beispiele, Fuzzing-Läufe oder GitHub-CI ausgeführt. Vorhandener
Code und Tests aus dem parallel gepushten Stand bleiben erhalten. Historische
CI-Nachweise anderer Pakete sind keine Freigabe dieser Änderungen.

Automatische Push-/PR-Workflows sind deaktiviert. Ausschließlich direkt auf
master arbeiten, keine Branches oder Pull Requests; neue Commits mit `[skip ci]`.

Modul `azphysics::xpbd`, ohne `alloc`, generisch über `f32`/`f64` und die Dimension `D`.
Der Aufrufer stellt alle Felder bereit (Positionen, Geschwindigkeiten, inverse Massen,
Bedingungen, zwei Arbeitsfelder); es wird nichts allokiert.

## Zeitschritt `xpbd_step`

1. **Vorhersage** beweglicher Teilchen: `v' = v + dt·g`, `x* = x + dt·v'`; fixierte Teilchen
   (`w = 0`) bleiben stehen.
2. **Iterationen** (Gauss–Seidel in Listenreihenfolge): für `C = |xᵢ − xⱼ| − L` und
   `α̃ = Compliance/dt²`: `Δλ = (−C − α̃λ)/(wᵢ + wⱼ + α̃)`, `λ += Δλ`, `xᵢ += wᵢΔλ·n`,
   `xⱼ −= wⱼΔλ·n`, `n = (xᵢ − xⱼ)/|xᵢ − xⱼ|`.
3. **Geschwindigkeit** aus der Positionsänderung `v = (x_neu − x_alt)/dt`.

## Verträge

| Punkt | Festlegung |
|---|---|
| Einheiten | m, s, m/s², `w` in 1/kg, Compliance in m/N (0 = starr), `λ` in kg·m |
| Lebensdauer von `λ` | pro Zeitschritt; Reset am Schrittanfang, **kein** Warmstart |
| Kraft einer Bedingung | `λ·∇C/dt²` mit `∇ᵢC = n` |
| Fixierte Teilchen | bewegen sich nie, Geschwindigkeit 0, Schwerkraft wirkt nicht |
| Koinzidente Punkte | Richtung undefiniert: Bedingung wird in dieser Iteration übersprungen |
| Zwei fixierte Teilchen | immer übersprungen, Lambda bleibt null; geometrisches Residuum bleibt sichtbar |
| Fehler | atomar: bei jedem Fehler bleibt der Zustand unverändert |
| Zeitschritt | `dt` endlich und `> 0`; negative Schritte sind hier nicht definiert |
| Nicht enthalten | Kontakte, räumliche Suche, Dämpfung, Reibung |

## Vorbereitete Abnahme (vorhandene Tests in `tests/xpbd.rs`, nicht ausgeführt)

- Zwei Teilchen, starr: exakte Projektion in einer Iteration; ungleiche Massen
  (Verschiebung ∝ `w`, Impulserhaltung `Σ v/w = 0`); Compliance mit exaktem Restfehler
  `C·α̃/(Σw + α̃)` und streng wachsendem Residuum bei wachsender Compliance.
- Fixierte und koinzidente Teilchen ohne Division durch null.
- Freier Fall entspricht der geschlossenen Form des Prädiktors.
- **Lineare Feder = implizites Euler:** für Compliance `1/k` ist ein XPBD-Schritt gleich
  `u' = (u + h·v)/(1 + h²k/m)`, unabhängig von der Iterationszahl; der Fehler gegenüber der exakten
  Lösung halbiert sich bei halbem `dt` (Ordnung 1).
- Kette mit Schwerkraft und fixiertem Ende: endlich, reproduzierbar, `x`-Komponenten exakt 0,
  begrenztes Residuum; mehr Iterationen erhöhen das Residuum nicht.
- Eingabefehler (`dt`, Iterationen, inverse Massen, Bedingungen, Feldlängen, `D = 0`, NaN/Inf,
  Überlauf) mit unverändertem Zustand.

Die Tests wurden in dieser Sitzung nicht ausgeführt. Die beschriebenen Erwartungen
sind analytische Kriterien, keine hier gemessenen Ergebnisse. Auch GitHub-CI wird
nicht gestartet; eine spätere Ausführung braucht einen ausdrücklichen neuen Auftrag.


## Numerische Fehlerbehandlung und Grenzen

Die vorhandene Slice-API (`xpbd_step`, `DistanceConstraint`, `XpbdStats`,
`max_constraint_error`) bleibt erhalten. Nichtendliche Differenzen oder Längen
werden als Fehler gemeldet und nicht als koinzidente Punkte verschluckt.
`alpha_tilde` wird als `(compliance/dt)/dt` gebildet; nicht darstellbare positive
Werte liefern `PhysicsError::NumericalRange`. Nichtendliche Nenner, Lambda-Werte,
Positionskorrekturen oder rekonstruierte Geschwindigkeiten führen zum Fehler,
ohne öffentliche Positionen/Geschwindigkeiten zu übernehmen. Arbeitsbuffer sind
bei Fehlern unspezifiziert. Ohne Bedingungen laufen null Projektionsdurchläufe.

`max_constraint_error` ist ein geometrisches Residuum in m. Bei Compliance ist
eine Dehnung erlaubt; `Ok` bedeutet einen endlichen Schritt und keine automatisch
kleinen Residuen. Die Constraint-Reihenfolge beeinflusst endliche Iterationsläufe.
Der Prädiktor und die Projektion versprechen keine Energieerhaltung oder Zeitumkehr.
Die Gleichheit mit implizitem Euler gilt für eine isolierte lineare Feder bei
fester radialer Richtung ohne Durchqueren des Fixpunkts, nicht für beliebige
gekoppelte/nichtlineare Systeme. Micromath unterstützt nur f32 und benötigt bei
einer späteren Freigabe eigene Genauigkeitsbudgets.

`examples/04_xpbd_chain.rs` ist vorbereitet und wurde nicht ausgeführt.
Kontakte, Reibung, CCD, Warmstarts und räumliche Kopplung gehören nicht zu AZP-040.
