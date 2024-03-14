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

double _physics_perihels[] =
{
    0.0, // sun
    0.0, // sun_system_barycenter
    0.308, // mercury
    0.718, // venus
    0.0, // earth
    363300.0e3 / PHYSICS_AU, // moon
    0.983, // earth-moon-barycenter
    1.382,
    4.95,
    9.063,
    18.295,
    29.707
};

double _physics_aphels[] =
{
    0.0, // sun
    0.0, // sun_system_barycenter
    0.467, // mercury
    0.728, // venus
    0.0, // earth
    405500.0e3 / PHYSICS_AU, // moon
    1.017, // earth-moon-barycenter
    1.666,
    5.458,
    10.101,
    20.107,
    30.387
};

double _physics_body_masses[] =
    {
        PHYSICS_SUN_MASS,
        PHYSICS_SUN_MASS + PHYSICS_MERCURY_MASS + PHYSICS_VENUS_MASS + PHYSICS_EARTH_MASS + PHYSICS_MOON_MASS + PHYSICS_MARS_MASS,
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

double _physics_initial_radii_AU[] =
    {
        0.0, // sun
        0.0, // sun_system_barycenter
        0.0,
        0.0,
        0.0, // earth
        0.0, // moon
        0.0, // earth-moon-barycenter
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    };

void physics_initialize()
{
    for (uint32_t i = 0; i < PHYSICS_BODY_COUNT; ++i)
    {
        switch (i)
        {
            case EARTH:
                _physics_initial_radii_AU[i] = physics_barycenter_AU(physics_body_aphel_AU(MOON), physics_body_mass(EARTH), physics_body_mass(MOON));
                break;

            case MOON:
                _physics_initial_radii_AU[i] = -(physics_body_aphel_AU(MOON) - _physics_initial_radii_AU[EARTH]);
                break;

            default:
                _physics_initial_radii_AU[i] = physics_body_aphel_AU(i);
                break;
        }
    }
}

double physics_body_perihel_AU(physics_body_id_t id)
{
    return _physics_perihels[id];
}

double physics_body_aphel_AU(physics_body_id_t id)
{
    return _physics_aphels[id];
}

double physics_initial_radii_AU(physics_body_id_t id)
{
    return _physics_initial_radii_AU[id];
}

const char* physics_body_name(physics_body_id_t id)
{
    return _physics_body_names[id];
}

double physics_body_mass(physics_body_id_t id)
{
    return _physics_body_masses[id];
}

double physics_pi()
{
    return acos(-1.0);
}

double physics_seconds_per_day()
{
    return 24.0 * 3600.0;
}

double physics_seconds_per_siderial_day()
{
    return 86164.0989;
}

double physics_barycenter_AU(double distance_AU, double mass_center_kg, double mass_satellite_kg)
{
    return distance_AU * mass_satellite_kg / (mass_center_kg + mass_satellite_kg);
}

struct vector_3d physics_newton_angular_speed(vector_3d_t r_AU, double center_mass_kg)
{
    struct vector_3d w = {0.0, 0.0, sqrt(center_mass_kg * PHYSICS_G / pow(vector_norm(r_AU) * PHYSICS_AU, 3.0))};
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
