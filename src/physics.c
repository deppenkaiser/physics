#include "physics/physics.h"
#include <vector/vector.h>
#include <math.h>

double _physics_body_masses[] =
    {
        PHYSICS_SUN_MASS,
        PHYSICS_SUN_MASS + PHYSICS_MERCURY_MASS + PHYSICS_VENUS_MASS + PHYSICS_EARTH_MASS + PHYSICS_MOON_MASS + PHYSICS_MARS_MASS,
        PHYSICS_MERCURY_MASS,
        PHYSICS_VENUS_MASS,
        PHYSICS_EARTH_MASS,
        PHYSICS_MOON_MASS,
        PHYSICS_EARTH_MASS + PHYSICS_MOON_MASS,
        PHYSICS_MARS_MASS
    };

const char* _physics_body_names[] =
    {
        "sun",
        "sun_barycenter",
        "mercury",
        "venus",
        "earth",
        "moon",
        "earth_moon_barycenter",
        "mars"
    };

const char* physics_body_name(physics_body_id_t id)
{
    return _physics_body_names[id];
}

double physics_moon_angular_speed()
{
    return 2.0 * physics_pi() / (physics_seconds_per_siderial_day() * 27.322);
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

double physics_barycenter(double distance_m, double mass_center_kg, double mass_satellite_kg)
{
    return distance_m * mass_satellite_kg / (mass_center_kg + mass_satellite_kg);
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
