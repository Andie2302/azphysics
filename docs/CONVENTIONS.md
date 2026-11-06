# Mechanik-Konventionen (AZP-010)

Stand: 2026-11-07. Diese erste Ausbaustufe beschreibt Punktmassen und autonome,
konservative Modelle. Implementierung: `azphysics::mechanics`.

## Größen und Koordinaten

| Größe | Einheit / Darstellung |
|---|---|
| Position q, Gleichgewichtslage | m; `[T; D]` |
| Geschwindigkeit v | m/s; `[T; D]` |
| Masse m | kg; endlicher, strikt positiver Skalar |
| Beschleunigung a | m/s²; expliziter Vektor |
| Federkonstante k | N/m = kg/s²; endlich und nichtnegativ |
| Zeit, Schrittweite | s; Skalar |
| Impuls p=m*v | kg*m/s; berechneter Vektor |
| Mechanische Energie | J; berechneter Skalar |

`D >= 1`. Üblich sind 1D (x), 2D (x,y) und ein rechtshändiges 3D-System (x,y,z).
Es gibt keine automatisch vorgegebene Gravitationsrichtung; z.B. kann ein
Beschleunigungsvektor `[0,0,-9.81]` gewählt werden. Größere D dienen entkoppelten
Koordinaten, nicht einer Behauptung über zusätzlichen physikalischen Raum.

Die SI-Bedeutung ist eine dokumentierte Konvention, keine Compile-Time-
Dimensionsanalyse. Es gibt vorerst keine Abhängigkeit von einer Einheitenbibliothek.

## Zustand und Parameter

`State<T,D>` hält **Position und Geschwindigkeit**, nicht Position und Impuls.
Der Konstruktor prüft endliche Werte; Felder sind privat, Getter liefern
unveränderliche Arrays. Für einen anderen Anfangszustand wird ein neuer
geprüfter State erzeugt. Impuls und kinetische Energie prüfen auch die
übergebene Masse. Modelle speichern ihre geprüfte Masse und weitere Parameter.

- `FreeParticle`: keine Kräfte; Gesamtenergie K.
- `ConstantAcceleration`: a ist bereits Beschleunigung, keine Kraft. Potential
  `V=-m*a·q` mit Nullpunkt im Koordinatenursprung; Gesamtenergie K+V.
- `HarmonicOscillator`: isotrop, `V=k*|q-q_eq|²/2`, `a=-(k/m)*(q-q_eq)`.
  `k=0` ist ein erlaubter freier Grenzfall. Nicht darstellbares positives k/m
  wird als NumericalRange abgewiesen.

Impuls ist beim freien Modell konstant. Unter äußerer Kraft gilt eine
Impulsbilanz; Impulserhaltung wird für Beschleunigungsfeld und Federkraft
nicht pauschal zugesagt. Energie-Nullpunkte der Modelle sind wie oben definiert.

## Zeitintegration und Fehler

Alle Modelle verwenden `azmath::ode::verlet`. Es gibt keinen zweiten
Integrator in azphysics. Aufruf: `integrate(state,t0,dt,steps,work,observer)`.
`mechanics_work_len(D)` liefert den nötigen Arbeitsbuffer. State, Modelle und
Kernintegration brauchen keinen Allokator.

- t0 muss endlich sein; dt muss endlich, ungleich null und als Halbschritt
  darstellbar sein; steps muss positiv sein.
- Ein negatives dt bedeutet ausdrücklich **Rückwärtsintegration**. Es ist kein
  Fehler, weil Zeitumkehr zu den Eigenschaften des Integrators gehört.
- Endzeit und Zeitfortschritt müssen darstellbar sein. Das Floating-Point-Raster
  wird von azmath aus den Zeitgrenzen gebildet; dt ist die nominelle Schrittweite.
- Der Observer erhält Anfangszustand und jeden vollständig abgeschlossenen
  Schritt. Sein State ist eine temporäre Kopie mit unveränderlichem Zugriff.
- Ungültige Uhren/Puffer ändern den Anfangszustand nicht. Bei einem Fehler im
  späteren Integrationsschritt bleibt der letzte vollständige Zustand erhalten.
- PhysicsError unterscheidet Parameter-, Zustands-/Rechenwert- und Zeitfehler;
  Integrationsfehler erhalten den zugrunde liegenden OdeError.

Für den harmonischen Oszillator gilt `|dt*omega|<2` als lineare Stabilitätsgrenze.
Für kleine Fehler sind deutlich kleinere Schritte nötig. Große Schritte werden
nicht automatisch verkleinert; ein numerisch endlicher Zustand ist kein
Stabilitäts- oder Genauigkeitsnachweis. Für glatte autonome Hamilton-Systeme
und geeignete Schritte wird das typische beschränkte Energiefehlerverhalten
geprüft. Für beliebige Kräfte oder Zeitabhängigkeiten gilt diese Aussage nicht.

Numerische Referenzprüfungen erfolgen für std/libm. Micromath bleibt ein
Näherungsbackend ohne präzise Physik-Freigabe. Endliche Eingaben garantieren
nicht die Darstellbarkeit aller Zwischenwerte; erkannte Überläufe und nicht
auflösbare Modellskalen werden diagnostiziert.
