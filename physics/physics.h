#pragma once

#include <stdint.h>
#include <vector/vector.h>

#define PHYSICS_C 299792458.0L
#define PHYSICS_G 6.6743015e-11L
#define PHYSICS_C_SQUARE (299792458.0L * 299792458.0L)
#define PHYSICS_AU 149597870700.0L
#define PHYSICS_CALCIUM_FILTER_WAVELENGTH_M 393.3e-9L
#define PHYSICS_GREEN_FILTER_M 540.0e-9L

// https://www.imcce.fr/

#define PHYSICS_SUN_MASS 1.98841e30L
#define PHYSICS_MERCURY_MASS 3.301e23L
#define PHYSICS_VENUS_MASS 4.8673e24L
#define PHYSICS_EARTH_MASS 5.9722e24L
#define PHYSICS_MOON_MASS 7.346e22L
#define PHYSICS_MARS_MASS 6.417e23L
#define PHYSICS_MARS_PHOBOS_MASS 1.072e16L
#define PHYSICS_MARS_DEIMOS_MASS 1.8e15L
#define PHYSICS_MARS_SYSTEM_MASS (PHYSICS_MARS_MASS + PHYSICS_MARS_PHOBOS_MASS + PHYSICS_MARS_DEIMOS_MASS)
#define PHYSICS_JUPITER_MASS 1.89813e27L
#define PHYSICS_JUPITER_GANYMED_MASS 1.4819e23L
#define PHYSICS_JUPITER_KALLISTO_MASS 1.0759e23L
#define PHYSICS_JUPITER_EUROPA_MASS 4.800e22L
#define PHYSICS_JUPITER_IO_MASS 8.93193797e22L
#define PHYSICS_JUPITER_SYSTEM_MASS (PHYSICS_JUPITER_MASS + PHYSICS_JUPITER_GANYMED_MASS + \
    PHYSICS_JUPITER_KALLISTO_MASS + PHYSICS_JUPITER_EUROPA_MASS + PHYSICS_JUPITER_IO_MASS)
#define PHYSICS_SATURN_MASS 5.683e26L
#define PHYSICS_SATURN_TITAN_MASS 1.345e23L
#define PHYSICS_SATURN_RHEA_MASS 2.31e21L
#define PHYSICS_SATURN_DIONE_MASS 1.10e21L
#define PHYSICS_SATURN_TETHYS_MASS 6.18e20L
#define PHYSICS_SATURN_IAPETUS_MASS 1.81e21L
#define PHYSICS_SATURN_SYSTEM_MASS (PHYSICS_SATURN_MASS + PHYSICS_SATURN_TITAN_MASS + PHYSICS_SATURN_RHEA_MASS + \
    PHYSICS_SATURN_DIONE_MASS + PHYSICS_SATURN_TETHYS_MASS + PHYSICS_SATURN_IAPETUS_MASS)
#define PHYSICS_URANUS_MASS 8.681e25L
#define PHYSICS_NEPTUNE_MASS 1.024e26L

#define PHYSICS_SUN_A 0.0L
#define PHYSICS_MERCURY_A 0.3870983098L
#define PHYSICS_VENUS_A 0.7233298200L
#define PHYSICS_EARTH_A 1.0000010178L
#define PHYSICS_MOON_A (384400.0e3L / PHYSICS_AU)
#define PHYSICS_MARS_A 1.5236793419L
#define PHYSICS_JUPITER_A 5.2026032092L
#define PHYSICS_SATURN_A 9.5549091915L
#define PHYSICS_URANUS_A 19.2184460618L
#define PHYSICS_NEPTUNE_A 30.1103868694L

#define PHYSICS_SUN_ECCENTRICITY 0.0L
#define PHYSICS_MERCURY_ECCENTRICITY 0.2056317526L
#define PHYSICS_VENUS_ECCENTRICITY 0.0067719164L
#define PHYSICS_EARTH_ECCENTRICITY 0.0167086342L
#define PHYSICS_MOON_ECCENTRICITY 0.0549L
#define PHYSICS_MARS_ECCENTRICITY 0.0934006477L
#define PHYSICS_JUPITER_ECCENTRICITY 0.0484979255L
#define PHYSICS_SATURN_ECCENTRICITY 0.0555481426L
#define PHYSICS_URANUS_ECCENTRICITY 0.0463812221L
#define PHYSICS_NEPTUNE_ECCENTRICITY 0.009455747L

#define PHYSICS_SUN_I 0.0L
#define PHYSICS_MERCURY_I 7.00498625L
#define PHYSICS_VENUS_I 3.39466189L
#define PHYSICS_EARTH_I 0.0L
#define PHYSICS_MOON_I 0.0L
#define PHYSICS_MARS_I 1.84972648L
#define PHYSICS_JUPITER_I 1.30326698L
#define PHYSICS_SATURN_I 2.48887878L
#define PHYSICS_URANUS_I 0.77319689L
#define PHYSICS_NEPTUNE_I 1.76995259L

// Aufsteigender Knoten
#define PHYSICS_SUN_NODE 0.0L
#define PHYSICS_MERCURY_NODE 48.33089304L
#define PHYSICS_VENUS_NODE 76.67992019L
#define PHYSICS_EARTH_NODE 0.0L
#define PHYSICS_MOON_NODE 0.0L
#define PHYSICS_MARS_NODE 49.55809321L
#define PHYSICS_JUPITER_NODE 100.46440702L
#define PHYSICS_SATURN_NODE 113.66550252L
#define PHYSICS_URANUS_NODE 74.00595701L
#define PHYSICS_NEPTUNE_NODE 131.78405702L

// Perihellänge
#define PHYSICS_SUN_PL 0.0L
#define PHYSICS_MERCURY_PL 77.45611904L
#define PHYSICS_VENUS_PL 131.56370300L
#define PHYSICS_EARTH_PL 102.93734808L
#define PHYSICS_MOON_PL 102.93734808L
#define PHYSICS_MARS_PL 336.06023395L
#define PHYSICS_JUPITER_PL 14.33120687L
#define PHYSICS_SATURN_PL 93.05723748L
#define PHYSICS_URANUS_PL 173.00529106L
#define PHYSICS_NEPTUNE_PL 48.12027554L

#define PHYSICS_SUN_W 0.0L
#define PHYSICS_MERCURY_W 29.125226L
#define PHYSICS_VENUS_W 54.88378281L
#define PHYSICS_EARTH_W 100.46645683L
#define PHYSICS_MOON_W 100.46645683L
#define PHYSICS_MARS_W 286.50214074L
#define PHYSICS_JUPITER_W 273.86679985L
#define PHYSICS_SATURN_W 339.39173496L
#define PHYSICS_URANUS_W 98.99933405L
#define PHYSICS_NEPTUNE_W 276.33621852L

// Mittlere Anomalie M0 zur Epoche J2000.0 (in Radiant)
// Quelle: IMCCE VSOP87 (https://www.imcce.fr/content/medias/recherche/equipes/asd/vsop87/vsop87.html)
#define PHYSICS_MERCURY_M0   3.050L          // VSOP87: 174.793° (Merkur)
#define PHYSICS_VENUS_M0     0.880L          // VSOP87: 50.416°  (Venus)
#define PHYSICS_EARTH_M0     6.240L          // VSOP87: 357.517° (Erde)
#define PHYSICS_MARS_M0      0.338L          // VSOP87: 19.373°  (Mars)
#define PHYSICS_JUPITER_M0   4.598L          // VSOP87: 263.451° (Jupiter)
#define PHYSICS_SATURN_M0    2.559L          // VSOP87: 146.636° (Saturn)
#define PHYSICS_URANUS_M0    4.482L          // VSOP87: 256.828° (Uranus)
#define PHYSICS_NEPTUNE_M0   3.196L          // VSOP87: 183.120° (Neptun)

typedef struct celestial_body
{
    struct vector_3d r_m;
    struct vector_3d r_bary_m;
    struct vector_3d v_m_s;
    struct vector_3d v_bary_m_s;
    struct vector_3d w_rad_s;
    ld h;
    ld K_1;
    ld K_2;
    ld A;
    ld mass_kg;
    ld e;
    ld e_square;
    ld a_m;
    ld T_s;
    ld perihel_m;
    const char* name;
    uint32_t index;
} *celestial_body_t;

ld physics_pi();
ld physics_frac(cld x);
ld physics_modulo(cld a, cld b);
ld physics_deg_to_rad(cld angle_deg);
ld physics_rad_to_deg(cld angle_rad);
ld physics_seconds_per_day();
ld physics_barycenter_AU(cld distance_AU, cld mass_center_kg, cld mass_satellite_kg);
ld physics_rayleigh_criteria_deg(cld wavelength_light_m, cld objective_aperture_m);
ld physics_image_sensor_object_size_m(cld object_size_deg, cld focal_length_m);
ld physics_needed_image_sensor_pixel_size_m(cld wavelength_light_m, cld objective_aperture_m, cld focal_length_m);
ld physics_kepler_radius(cld a_m, cld eccentricity, cld phi_rad);
ld physics_kinetic_energy(cld mass_kg, const vector_3d_t v);
ld physics_kinetic_energy_body(const celestial_body_t body);

ld physics_weber_potential_energy(const celestial_body_t body, cld mass_center_kg, cld phi_rad);
ld physics_weber_h(const celestial_body_t body, cld mass_center_kg);
ld physics_weber_alpha(const celestial_body_t body, cld mass_center_kg);
ld physics_weber_k_1st_order(const celestial_body_t body, cld mass_center_kg);
ld physics_weber_k_2nd_order(const celestial_body_t body, cld mass_center_kg);
ld physics_weber_periodtime(const celestial_body_t body, cld mass_center_kg);
ld physics_weber_deltaphi(const celestial_body_t body, cld mass_center_kg, cld t_step_s, cld phi_0_rad);
ld physics_deltaphi_per_revolution_1(const celestial_body_t body, cld mass_center_kg);
ld physics_deltaphi_per_revolution_2(const celestial_body_t body, cld mass_center_kg);
struct vector_3d physics_weber_position(const celestial_body_t body, cld mass_center_kg, cld phi_rad);
struct vector_3d physics_weber_angular_speed(const celestial_body_t body, cld mass_center_kg, cld phi_rad);
