#include "physics/physics.h"
#include <vector/vector.h>
#include <math.h>
#include <stdint.h>

#define PHYSICS_SUN_MASS 1.9884e30
#define PHYSICS_MERCURY_MASS 3.301e23
#define PHYSICS_VENUS_MASS 4.8673e24
#define PHYSICS_EARTH_MASS 5.9722e24
#define PHYSICS_MOON_MASS 7.346e22
#define PHYSICS_MARS_MASS 6.417e23
#define PHYSICS_JUPITER_MASS 1.89813e27
#define PHYSICS_SATURN_MASS 5.683e26
#define PHYSICS_URANUS_MASS 8.681e25
#define PHYSICS_NEPTUNE_MASS 1.024e26

#define PHYSICS_SUN_P 0.0
#define PHYSICS_MERCURY_P 0.3876
#define PHYSICS_VENUS_P 0.7233
#define PHYSICS_EARTH_P 1.0
#define PHYSICS_MOON_P (384400.0e3 / PHYSICS_AU)
#define PHYSICS_MARS_P 1.524
#define PHYSICS_JUPITER_P 5.204
#define PHYSICS_SATURN_P 9.582
#define PHYSICS_URANUS_P 19.201
#define PHYSICS_NEPTUNE_P 30.047

#define PHYSICS_SUN_ECCENTRICITY 0.0
#define PHYSICS_MERCURY_ECCENTRICITY 0.2056
#define PHYSICS_VENUS_ECCENTRICITY 0.0067
#define PHYSICS_EARTH_ECCENTRICITY 0.0167
#define PHYSICS_MOON_ECCENTRICITY 0.0549
#define PHYSICS_MARS_ECCENTRICITY 0.0934
#define PHYSICS_JUPITER_ECCENTRICITY 0.0489
#define PHYSICS_SATURN_ECCENTRICITY 0.0542
#define PHYSICS_URANUS_ECCENTRICITY 0.0472
#define PHYSICS_NEPTUNE_ECCENTRICITY 0.0113

double _physics_body_p_AU[] =
    {
        PHYSICS_SUN_P,
        0,
        PHYSICS_MERCURY_P,
        PHYSICS_VENUS_P,
        PHYSICS_EARTH_P,
        PHYSICS_MOON_P,
        0,
        PHYSICS_MARS_P,
        PHYSICS_JUPITER_P,
        PHYSICS_SATURN_P,
        PHYSICS_URANUS_P,
        PHYSICS_NEPTUNE_P
    };

double _physics_body_eccentricity[] =
    {
        PHYSICS_SUN_ECCENTRICITY,
        0,
        PHYSICS_MERCURY_ECCENTRICITY,
        PHYSICS_VENUS_ECCENTRICITY,
        PHYSICS_EARTH_ECCENTRICITY,
        PHYSICS_MOON_ECCENTRICITY,
        0,
        PHYSICS_MARS_ECCENTRICITY,
        PHYSICS_JUPITER_ECCENTRICITY,
        PHYSICS_SATURN_ECCENTRICITY,
        PHYSICS_URANUS_ECCENTRICITY,
        PHYSICS_NEPTUNE_ECCENTRICITY
    };

double _physics_body_masses_kg[] =
    {
        PHYSICS_SUN_MASS,
        PHYSICS_SUN_MASS + PHYSICS_MERCURY_MASS + PHYSICS_VENUS_MASS + PHYSICS_EARTH_MASS + PHYSICS_MOON_MASS +
            PHYSICS_MARS_MASS + PHYSICS_JUPITER_MASS + PHYSICS_SATURN_MASS + PHYSICS_URANUS_MASS + PHYSICS_NEPTUNE_MASS,
        PHYSICS_MERCURY_MASS,
        PHYSICS_VENUS_MASS,
        PHYSICS_EARTH_MASS,
        PHYSICS_MOON_MASS,
        PHYSICS_EARTH_MASS + PHYSICS_MOON_MASS,
        PHYSICS_MARS_MASS,
        PHYSICS_JUPITER_MASS,
        PHYSICS_SATURN_MASS,
        PHYSICS_URANUS_MASS,
        PHYSICS_NEPTUNE_MASS
    };

const char* _physics_body_names[] =
    {
        "sun",
        "sun_system_barycenter",
        "mercury",
        "venus",
        "earth",
        "moon",
        "earth_moon_barycenter",
        "mars",
        "jupiter",
        "saturn",
        "uranus",
        "neptune"
    };

void physics_initialize()
{
}

double physics_body_perihel_AU(physics_body_id_t id)
{
    return 0.0;
}

double physics_body_aphel_AU(physics_body_id_t id)
{
    return 0.0;
}

const char* physics_body_name(physics_body_id_t id)
{
    return _physics_body_names[id];
}

double physics_body_mass_kg(physics_body_id_t id)
{
    return _physics_body_masses_kg[id];
}

double physics_body_p_AU(physics_body_id_t id)
{
    return _physics_body_p_AU[id];
}

double physics_body_eccentricity(physics_body_id_t id)
{
    return _physics_body_eccentricity[id];
}

double physics_body_angle_rad(double t_s, double epoche_rad, double p_AU)
{
    return sqrt(PHYSICS_G * _physics_body_masses_kg[SUN_BARYCENTER] / pow(p_AU * PHYSICS_AU, 3.0)) * t_s + epoche_rad;
}

double physics_pi()
{
    return acos(-1.0);
}

double physics_seconds_per_day()
{
    return 24.0 * 3600.0;
}

double physics_barycenter_AU(double distance_AU, double mass_center_kg, double mass_satellite_kg)
{
    return distance_AU * mass_satellite_kg / (mass_center_kg + mass_satellite_kg);
}

struct vector_3d physics_newton_angular_speed(vector_3d_t r_AU, double center_mass_kg)
{
    struct vector_3d w = {0};
    return w;
}

struct vector_3d physics_weber_gravitation_force(double m1_kg, double m2_kg, vector_3d_t r, vector_3d_t v, vector_3d_t a)
{
    struct vector_3d F = {0};
    return F;
}

struct vector_3d physics_tangential_speed(vector_3d_t w, vector_3d_t r)
{
    return vector_cross(w, r);
}

struct vector_3d physics_linear_motion(vector_3d_t v, double t_s)
{
    return vector_multiply_scalar(v, t_s);
}

struct vector_3d physics_linear_angle(vector_3d_t w, double t_s)
{
    return vector_multiply_scalar(w, t_s);
}

struct vector_3d physics_kepler_position_AU(physics_body_id_t id, double angle_rad)
{
    double R_AU = physics_body_p_AU(id) / (physics_body_eccentricity(id) * cos(angle_rad) + 1.0);
    struct vector_3d r_AU = {cos(angle_rad) * R_AU, sin(angle_rad) * R_AU, 0.0};
    return r_AU;
}
