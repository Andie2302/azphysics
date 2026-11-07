# Roadmap

Die zentrale Planung liegt in
[aztodo](https://github.com/Andie2302/aztodo/blob/master/ROADMAP.md).

- AZP-001: Ursprüngliches Gerüst 0.0.1 vorhanden; die ersten Modelle stehen unter AZP-020.
- AZP-002: Direkte azmath-Git-Abhängigkeit mit festem rev und Backend-/alloc-Weiterleitung geprüft.
- AZP-010: SI-Konventionen, geprüfter Zustand und Fehlerdiagnosen implementiert.
- AZP-020: Freier Flug, konstante Beschleunigung und harmonischer Oszillator mit azmath-Verlet implementiert.
- AZM-031: Gauss–Legendre-RK4 aus azmath mit Verlet-Vergleich (Beispiel und f32/f64-Tests).
- GEO-001: Partikelindex in azphysics::spatial; mathematische Rotation/Transformation bleibt in azmath (docs/SPATIAL.md).
- AZP-030: Spatial Hash mit alloc für 2D/3D, stabile IDs, Updates und deterministische Kandidaten.
- AZP-031: BVH separat; AZP-040/041: XPBD-Distanzbedingungen und Kontakte.
- AZP-051/050: Druckprojektion und APIC/FLIP.
- AZP-060/070: Materialmodelle und spätere weitere Physikgebiete.

Verbindliche Voraussetzungen und Abnahmekriterien stehen im
[Backlog](https://github.com/Andie2302/aztodo/blob/master/BACKLOG.md).
