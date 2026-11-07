# Räumliche Partikelsuche (GEO-001 / AZP-030)

## Entscheidung und Bestand

azmath enthält allgemeine Vektoren sowie `geometry` für Quaternionen, SO(3),
SE(3), starre Transformationen und Euler-Winkel. Diese Bausteine bleiben dort.
Für den ersten Partikelindex werden keine neuen AABB-/Ray-/Mesh-Primitiven benötigt.

Zwei API-Varianten wurden für GEO-001 verglichen:

| Variante | API-Skizze | Folgen |
|---|---|---|
| Allgemeiner mathematischer Zellindex in azmath | `GridIndex::build(cells, ids)`; Position/Zellzuordnung extern | Gut für generische Ganzzahl-Zellabfragen, würde zunächst ohne weiteren Nutzer Positionen, Zellgrenzen und Update-Vertrag auf zwei Crates verteilen |
| Dynamischer Partikelindex in azphysics | `SpatialHash<T,D>::insert(Particle { id, position })`, `update`, `query_candidates` | Besitzt die für Simulationsschritte benötigten Positionskopien/IDs und deren Update-Vertrag gemeinsam; nur Real wird aus azmath verwendet |

Gewählt ist die zweite Variante: `azphysics::spatial`, optionales `alloc`,
weiterhin `no_std`, f32/f64 und zunächst genau 2D/3D. Es entsteht keine Abhängigkeit
von azmath zu azphysics und kein viertes Repository. Allgemeine AABB-Primitiven
werden erst beim getrennten BVH-Paket AZP-031 anhand tatsächlicher Anforderungen
zugeordnet. BVH ist keine Voraussetzung für diesen Index oder XPBD-Distanzbedingungen.

## Vertrag

- Kartesische Positionskopien, stabile externe `usize`-IDs. Eine ID kommt nur
  einmal vor; nach Entfernung darf der Aufrufer sie neu vergeben. Generationen
  alter Handles werden nicht automatisch verfolgt.
- Zellgröße strikt positiv/endlich; Positionen endlich und
  `|position/cell_size| <= 2^50`. Diese Grenze gilt für gespeicherte Punkte und
  Abfragezentren. Division/Floor-Zuordnung erfolgt kontrolliert in f64, auch bei
  f32 und micromath; Float→Integer-Sättigung wird nicht als gültige Zelle übernommen.
- Koordinaten/Radien/Zellgröße in derselben Längeneinheit, in Physikbeispielen m.
  Radius endlich und >= 0. Negative Koordinaten verwenden Floor, nicht Abschneiden.
- `query_candidates(center, radius)` liefert aufsteigend sortierte eindeutige IDs
  einer konservativen Zellbox. Kugelnachbarn einschließlich Grenzberührung müssen
  enthalten sein. Eine zusätzliche Zelle an jeder Grenze schützt gegen Rundung;
  dadurch entstehen absichtlich zusätzliche Kandidaten, auch bei Radius null.
- `candidate_pairs(radius)` liefert lexikographisch sortierte eindeutige Paare
  `(a,b)` mit `a < b`, ohne Selbstpaare. Das ist eine Broad Phase, keine exakte
  Nachbarschafts- oder Kontaktprüfung. Für Teilchenradien zunächst einen gemeinsamen
  konservativen Suchabstand wählen; individuelle Radien werden nicht gespeichert.
- Hash-Kollisionen vergleichen vollständige Zellkoordinaten; IDs/Ergebnisse hängen
  nicht von Einfügefolge oder Bucketzahl ab. Große Query-Boxen scannen die Partikel
  statt unbeschränkt viele leere Zellen zu durchlaufen.
- Updates/Entfernung wirken sofort; `rebuild` ersetzt den Index atomar erst nach
  erfolgreicher Validierung sämtlicher Einträge. Bei `SpatialError` bleibt der
  bisherige Index erhalten. Allokationsversagen folgt dem normalen Rust-Verhalten.

## Grenzen und Kosten

Der Hash verwendet separat verkettete Buckets und eine BTreeMap für IDs/Positionen.
Bucketzahl ist eine Zweierpotenz und wächst bei mittlerer Belegung über zwei
Partikeln; ID-Zugriff O(log n), Speicher O(n + Buckets). Queries enumerieren
höchstens 4096 Zellen, sonst wird der Bestand konservativ gescannt. Hash-Bucket-
Kollisionen kosten zusätzliche ID-/Zellprüfungen. Ergebnislisten werden sortiert
und allokiert. Es gibt keine garantierte konstante Abfragezeit; Cluster oder
große Radien können quadratisch viele Paare verursachen. Keine Parallelisierung,
periodischen Ränder, CCD, AABB-/Mesh-/Ray-Abfragen oder no-alloc-Variante in AZP-030.

## Nachweise

`tests/spatial_hash.rs` enthält 38 Tests: je neun für f32/f64 × 2D/3D sowie
Konstruktorfehler und drei Szenengrößen. Feste Seeds, Brute-Force-Radius-/Paarreferenzen,
Zell-/Floatgrenzen, negative Koordinaten, Koinzidenz, Hash-Kollisionen, Updates,
Entfernung, atomarer Neuaufbau, stabile Reihenfolge, NaN/Inf und Bereichsgrenzen.
Der Test läuft mit alloc und jedem Backend, einschließlich micromath; es werden
keine näherungsweisen Transzendentalfunktionen für die Zellzuordnung verwendet.

`cargo run --locked --features alloc --example 03_spatial_hash` zeigt bei
64/256/1024 Partikeln Bucketzahl, alle möglichen Paare, Kandidaten und echte
Nachbarn; jeder echte Nachbar wird gegen die Kandidaten geprüft. Die CI führt
dieses Beispiel mit alloc und allen Backends aus. Die Anzahl ist ein struktureller
Vergleich, kein Laufzeitbenchmark oder allgemeines Skalierungsversprechen.

Tests werden gemäß Projektwunsch nicht lokal ausgeführt. CI- und Merge-Nachweise
werden nach erfolgreichem Lauf in aztodo/STATUS.md festgehalten.
