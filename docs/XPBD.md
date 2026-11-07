# XPBD-Distanzbedingungen (AZP-040)

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
| Zwei fixierte Teilchen, starr | `Σw + α̃ = 0`: übersprungen |
| Fehler | atomar: bei jedem Fehler bleibt der Zustand unverändert |
| Zeitschritt | `dt` endlich und `> 0`; negative Schritte sind hier nicht definiert |
| Nicht enthalten | Kontakte, räumliche Suche, Dämpfung, Reibung |

## Nachweise (Tests in `tests/xpbd.rs`)

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

Die Tests wurden in der Arbeitssitzung nicht lokal ausgeführt; die Erwartungswerte sind
analytisch hergeleitet. Maßgeblich ist die CI.
