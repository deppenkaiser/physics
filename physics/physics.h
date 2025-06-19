#pragma once

#include <stdint.h>
#include <vector/vector.h>

#define PHYSICS_C 299792458.0
#define PHYSICS_G 6.6743015e-11
#define PHYSICS_C_SQUARE (299792458.0 * 299792458.0)
#define PHYSICS_AU 149597870700.0
#define PHYSICS_CALCIUM_FILTER_WAVELENGTH_M 393.3e-9
#define PHYSICS_GREEN_FILTER_M 540.0e-9

// https://www.imcce.fr/

#define PHYSICS_SUN_MASS 1.98841e30
#define PHYSICS_MERCURY_MASS 3.301e23
#define PHYSICS_VENUS_MASS 4.8673e24
#define PHYSICS_EARTH_MASS 5.9722e24
#define PHYSICS_MOON_MASS 7.346e22
#define PHYSICS_MARS_MASS 6.417e23
#define PHYSICS_MARS_PHOBOS_MASS 1.072e16
#define PHYSICS_MARS_DEIMOS_MASS 1.8e15
#define PHYSICS_MARS_SYSTEM_MASS (PHYSICS_MARS_MASS + PHYSICS_MARS_PHOBOS_MASS + PHYSICS_MARS_DEIMOS_MASS)
#define PHYSICS_JUPITER_MASS 1.89813e27
#define PHYSICS_JUPITER_GANYMED_MASS 1.4819e23
#define PHYSICS_JUPITER_KALLISTO_MASS 1.0759e23
#define PHYSICS_JUPITER_EUROPA_MASS 4.800e22
#define PHYSICS_JUPITER_IO_MASS 8.93193797e22
#define PHYSICS_JUPITER_SYSTEM_MASS (PHYSICS_JUPITER_MASS + PHYSICS_JUPITER_GANYMED_MASS + \
    PHYSICS_JUPITER_KALLISTO_MASS + PHYSICS_JUPITER_EUROPA_MASS + PHYSICS_JUPITER_IO_MASS)
#define PHYSICS_SATURN_MASS 5.683e26
#define PHYSICS_SATURN_TITAN_MASS 1.345e23
#define PHYSICS_SATURN_RHEA_MASS 2.31e21
#define PHYSICS_SATURN_DIONE_MASS 1.10e21
#define PHYSICS_SATURN_TETHYS_MASS 6.18e20
#define PHYSICS_SATURN_IAPETUS_MASS 1.81e21
#define PHYSICS_SATURN_SYSTEM_MASS (PHYSICS_SATURN_MASS + PHYSICS_SATURN_TITAN_MASS + PHYSICS_SATURN_RHEA_MASS + \
    PHYSICS_SATURN_DIONE_MASS + PHYSICS_SATURN_TETHYS_MASS + PHYSICS_SATURN_IAPETUS_MASS)
#define PHYSICS_URANUS_MASS 8.681e25
#define PHYSICS_NEPTUNE_MASS 1.024e26

#define PHYSICS_SUN_A 0.0
#define PHYSICS_MERCURY_A 0.3870983098
#define PHYSICS_VENUS_A 0.7233298200
#define PHYSICS_EARTH_A 1.0000010178
#define PHYSICS_MOON_A (384400.0e3 / PHYSICS_AU)
#define PHYSICS_MARS_A 1.5236793419
#define PHYSICS_JUPITER_A 5.2026032092
#define PHYSICS_SATURN_A 9.5549091915
#define PHYSICS_URANUS_A 19.2184460618
#define PHYSICS_NEPTUNE_A 30.1103868694

#define PHYSICS_SUN_ECCENTRICITY 0.0
#define PHYSICS_MERCURY_ECCENTRICITY 0.2056317526
#define PHYSICS_VENUS_ECCENTRICITY 0.0067719164
#define PHYSICS_EARTH_ECCENTRICITY 0.0167086342
#define PHYSICS_MOON_ECCENTRICITY 0.0549
#define PHYSICS_MARS_ECCENTRICITY 0.0934006477
#define PHYSICS_JUPITER_ECCENTRICITY 0.0484979255
#define PHYSICS_SATURN_ECCENTRICITY 0.0555481426
#define PHYSICS_URANUS_ECCENTRICITY 0.0463812221
#define PHYSICS_NEPTUNE_ECCENTRICITY 0.009455747

#define PHYSICS_SUN_I 0.0
#define PHYSICS_MERCURY_I 7.00498625
#define PHYSICS_VENUS_I 3.39466189
#define PHYSICS_EARTH_I 0.0
#define PHYSICS_MOON_I 0.0
#define PHYSICS_MARS_I 1.84972648
#define PHYSICS_JUPITER_I 1.30326698
#define PHYSICS_SATURN_I 2.48887878
#define PHYSICS_URANUS_I 0.77319689
#define PHYSICS_NEPTUNE_I 1.76995259

// Aufsteigender Knoten
#define PHYSICS_SUN_NODE 0.0
#define PHYSICS_MERCURY_NODE 48.33089304
#define PHYSICS_VENUS_NODE 76.67992019
#define PHYSICS_EARTH_NODE 0.0
#define PHYSICS_MOON_NODE 0.0
#define PHYSICS_MARS_NODE 49.55809321
#define PHYSICS_JUPITER_NODE 100.46440702
#define PHYSICS_SATURN_NODE 113.66550252
#define PHYSICS_URANUS_NODE 74.00595701
#define PHYSICS_NEPTUNE_NODE 131.78405702

// Perihellänge
#define PHYSICS_SUN_PL 0.0
#define PHYSICS_MERCURY_PL 77.45611904
#define PHYSICS_VENUS_PL 131.56370300
#define PHYSICS_EARTH_PL 102.93734808
#define PHYSICS_MOON_PL 102.93734808
#define PHYSICS_MARS_PL 336.06023395
#define PHYSICS_JUPITER_PL 14.33120687
#define PHYSICS_SATURN_PL 93.05723748
#define PHYSICS_URANUS_PL 173.00529106
#define PHYSICS_NEPTUNE_PL 48.12027554

#define PHYSICS_SUN_W 0.0
#define PHYSICS_MERCURY_W 29.125226
#define PHYSICS_VENUS_W 54.88378281
#define PHYSICS_EARTH_W 100.46645683
#define PHYSICS_MOON_W 100.46645683
#define PHYSICS_MARS_W 286.50214074
#define PHYSICS_JUPITER_W 273.86679985
#define PHYSICS_SATURN_W 339.39173496
#define PHYSICS_URANUS_W 98.99933405
#define PHYSICS_NEPTUNE_W 276.33621852

// Mittlere Anomalie M0 zur Epoche J2000.0 (in Radiant)
// Quelle: IMCCE VSOP87 (https://www.imcce.fr/content/medias/recherche/equipes/asd/vsop87/vsop87.html)
#define PHYSICS_MERCURY_M0   3.050          // VSOP87: 174.793° (Merkur)
#define PHYSICS_VENUS_M0     0.880          // VSOP87: 50.416°  (Venus)
#define PHYSICS_EARTH_M0     6.240          // VSOP87: 357.517° (Erde)
#define PHYSICS_MARS_M0      0.338          // VSOP87: 19.373°  (Mars)
#define PHYSICS_JUPITER_M0   4.598          // VSOP87: 263.451° (Jupiter)
#define PHYSICS_SATURN_M0    2.559          // VSOP87: 146.636° (Saturn)
#define PHYSICS_URANUS_M0    4.482          // VSOP87: 256.828° (Uranus)
#define PHYSICS_NEPTUNE_M0   3.196          // VSOP87: 183.120° (Neptun)

typedef struct celestial_body
{
    struct vector_3d r_m;
    struct vector_3d v_m_s;
    struct vector_3d w_rad_s;
    double mass_kg;
    double e;
    double a_m;
    const char* name;
} *celestial_body_t;

typedef const double cd;

double physics_pi();
double physics_frac(cd x);
double physics_modulo(cd a, cd b);
double physics_deg_to_rad(cd angle_deg);
double physics_rad_to_deg(cd angle_rad);
double physics_seconds_per_day();
double physics_barycenter_AU(cd distance_AU, cd mass_center_kg, cd mass_satellite_kg);
double physics_rayleigh_criteria_deg(cd wavelength_light_m, cd objective_aperture_m);
double physics_image_sensor_object_size_m(cd object_size_deg, cd focal_length_m);
double physics_needed_image_sensor_pixel_size_m(cd wavelength_light_m, cd objective_aperture_m, cd focal_length_m);
double physics_kepler_radius(cd a_m, cd eccentricity, cd phi_rad);
double physics_kinetic_energy(cd mass_kg, const vector_3d_t v);
double physics_kinetic_energy_body(const celestial_body_t body);

double physics_weber_potential_energy(const celestial_body_t body, cd mass_center_kg, cd phi_rad);
double physics_weber_specific_angular_momentum(const celestial_body_t body, cd mass_center_kg);
struct vector_3d physics_weber_position(const celestial_body_t body, cd mass_center_kg, cd phi_rad);
struct vector_3d physics_weber_position_perturbed(const celestial_body_t body, const celestial_body_t body_i, cd mass_center_kg);
struct vector_3d physics_weber_velocity(const celestial_body_t body, cd mass_center_kg, cd phi_rad);
struct vector_3d physics_weber_velocity_perturbed(const celestial_body_t body, const celestial_body_t body_i, cd mass_center_kg);
struct vector_3d physics_weber_angular_speed(const celestial_body_t body, cd mass_center_kg, cd phi_rad);
struct vector_3d physics_weber_angular_speed_perturbed(const celestial_body_t body, const celestial_body_t body_i, cd mass_center_kg);
