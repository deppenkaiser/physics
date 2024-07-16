#pragma once

#include <stdint.h>
#include <vector/vector.h>

#define PHYSICS_C 299792458.0
#define PHYSICS_G 6.6743015e-11
#define PHYSICS_C_SQUARE (299792458.0 * 299792458.0)
#define PHYSICS_AU 149597870700.0

typedef enum
{
    SUN = 0,
    SUN_BARYCENTER = 1,
    MERCURY = 2,
    VENUS = 3,
    EARTH = 4,
    MOON = 5,
    EARTH_MOON_BARYCENTER = 6,
    MARS = 7,
    JUPITER = 8,
    SATURN = 9,
    URANUS = 10,
    NEPTUNE = 11,
    PHYSICS_BODY_COUNT = 12
} physics_body_id_t;

const char* physics_body_name(physics_body_id_t id);
double physics_body_mass_kg(physics_body_id_t id);
double physics_body_a_AU(physics_body_id_t id);
double physics_body_eccentricity(physics_body_id_t id);
double physics_body_i_rad(physics_body_id_t id);
double physics_body_node_rad(physics_body_id_t id);
double physics_body_w_rad(physics_body_id_t id);
double physics_body_perihel_AU(physics_body_id_t id);
double physics_body_aphel_AU(physics_body_id_t id);
double physics_czybor_angular_speed_rad_per_s(physics_body_id_t planet_id, physics_body_id_t center_id, double phi_rad);
double physics_czybor_mean_angular_speed_rad_per_s(physics_body_id_t planet_id, physics_body_id_t center_id);
double physics_pi();
double physics_frac(double x);
double physics_modulo(double a, double b);
double physics_deg_to_rad();
double physics_rad_to_deg();
double physics_seconds_per_day();
double physics_barycenter_AU(double distance_AU, double mass_center_kg, double mass_satellite_kg);
double physics_orbital_period_s(physics_body_id_t planet_id, physics_body_id_t center_id, double phi_1_rad, double phi_2_rad, uint32_t step_count);
double physics_optimized_orbital_period_s(physics_body_id_t planet_id, physics_body_id_t center_id, double phi_1_rad, double phi_2_rad);
double physics_full_orbital_period_days(physics_body_id_t planet_id, physics_body_id_t center_id);
struct vector_3d physics_kepler_r_AU(physics_body_id_t id, double angle_rad);
