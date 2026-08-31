# physics – Physik-Bibliothek

C-Bibliothek für physikalische Berechnungen: Himmelsmechanik,
Kepler-Orbit, Weber-Gravitation, Rayleigh-Kriterium und Konstanten.

## Umfang

- Physikalische Konstanten (G, c, AU, Planetenmassen, Orbitalparameter)
- Himmelsmechanik: Kepler-Radius, Barycenter, Umlaufzeit
- Weber-Gravitation: Periodzeit, Position, Winkelgeschwindigkeit
- Optik: Rayleigh-Kriterium, Sensor-Pixel-Größe
- Hilfsfunktionen: `physics_modulo`, `physics_deg_to_rad`, `physics_rad_to_deg`

## Nutzung

```c
#include <physics/physics.h>

ld earth_period = physics_weber_periodtime(&earth_body, PHYSICS_SUN_MASS);
ld rayleigh = physics_rayleigh_criteria_deg(550e-9, 0.1);
```

## Abhängigkeiten

- `string`, `logging`, `vector`
- System: GSL (GNU Scientific Library)

## Build

```bash
cmake -S . -B build
cmake --build build
```

In ein Projekt einbinden: `add_subdirectory(../../libraries/physics …)`,
Einbindung des Headers über den Include-Pfad `<physics/physics.h>`.

## Wiki

Dokumentation: [physics – Physik-Bibliothek](https://czybor.i234.me/wiki/sw-module/physics/) (Quartz-Wiki)
