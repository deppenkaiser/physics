#pragma once

#include <vector/vector.h>

#define PHYSICS_C 299792458.0
#define PHYSICS_C_SQUARE (299792458.0 * 299792458.0)
#define PHYSICS_AU 149597870700.0
#define PHYSICS_SUN_MASS 1.9884e30
#define PHYSICS_MERCURY_MASS 3.301e23
#define PHYSICS_VENUS_MASS 4.8673e24
#define PHYSICS_EARTH_MASS 5.9722e24
#define PHYSICS_MOON_MASS 7.346e22
#define PHYSICS_MOON_ANGULAR_SPEED (2 * physics_pi() / (physics_seconds_per_siderial_day() * 27.322))
#define PHYSICS_EARTH_MOON_DISTANCE 384400.0e3

double physics_pi();
double physics_seconds_per_day();
double physics_seconds_per_siderial_day();
struct vector_3d physics_weber_gravitation_force(double m1_kg, double m2_kg, vector_3d_t r, vector_3d_t v, vector_3d_t a);
struct vector_3d physics_tangential_speed(vector_3d_t w, vector_3d_t r);
struct vector_3d physics_linear_motion(vector_3d_t v, double t_s);
double physics_barycenter(double distance_m, double mass_center_kg, double mass_satellite_kg);
