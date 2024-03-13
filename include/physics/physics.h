#pragma once

#include <vector/vector.h>

#define PHYSICS_C 299792458.0
#define PHYSICS_G 6.6743015e-11
#define PHYSICS_C_SQUARE (299792458.0 * 299792458.0)
#define PHYSICS_AU 149597870700.0
#define PHYSICS_SUN_MASS 1.9884e30
#define PHYSICS_MERCURY_MASS 3.301e23
#define PHYSICS_VENUS_MASS 4.8673e24
#define PHYSICS_EARTH_MASS 5.9722e24
#define PHYSICS_MOON_MASS 7.346e22
#define PHYSICS_MARS_MASS 6.417e23
#define PHYSICS_MOON_APOAPSIS 405500.0e3

typedef enum
{
    SUN = 0,
    SUN_BARYCENTER = 1,
    MERCURY = 2,
    VENUS = 3,
    EARTH = 4,
    MOON = 5,
    EARTH_MOON_BARYCENTER = 6,
    MARS = 7
} physics_body_id_t;

const char* physics_body_name(physics_body_id_t id);
double physics_body_mass(physics_body_id_t id);
double physics_pi();
double physics_seconds_per_day();
double physics_seconds_per_siderial_day();
double physics_moon_angular_speed();
double physics_newton_angular_speed(vector_3d_t r, double mass_kg);
struct vector_3d physics_weber_gravitation_force(double m1_kg, double m2_kg, vector_3d_t r, vector_3d_t v, vector_3d_t a);
struct vector_3d physics_tangential_speed(vector_3d_t w, vector_3d_t r);
struct vector_3d physics_linear_motion(vector_3d_t v, double t_s);
double physics_barycenter(double distance_m, double mass_center_kg, double mass_satellite_kg);
