#include "physics/physics.h"
#include <vector/vector.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

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

#define PHYSICS_SUN_A 0.0
#define PHYSICS_MERCURY_A 0.38709893
#define PHYSICS_VENUS_A 0.72333199
#define PHYSICS_EARTH_A 1.00000011
#define PHYSICS_MOON_A (384400.0e3 / PHYSICS_AU)
#define PHYSICS_MARS_A 1.52366231
#define PHYSICS_JUPITER_A 5.20336301
#define PHYSICS_SATURN_A 9.53707032
#define PHYSICS_URANUS_A 19.19126393
#define PHYSICS_NEPTUNE_A 30.06896348

#define PHYSICS_SUN_ECCENTRICITY 0.0
#define PHYSICS_MERCURY_ECCENTRICITY 0.20563069
#define PHYSICS_VENUS_ECCENTRICITY 0.00677323
#define PHYSICS_EARTH_ECCENTRICITY 0.01671022
#define PHYSICS_MOON_ECCENTRICITY 0.0549
#define PHYSICS_MARS_ECCENTRICITY 0.09341233
#define PHYSICS_JUPITER_ECCENTRICITY 0.04839266
#define PHYSICS_SATURN_ECCENTRICITY 0.05415060
#define PHYSICS_URANUS_ECCENTRICITY 0.04716771
#define PHYSICS_NEPTUNE_ECCENTRICITY 0.00858587

#define PHYSICS_SUN_I 0.0
#define PHYSICS_MERCURY_I 7.00487
#define PHYSICS_VENUS_I 3.39471
#define PHYSICS_EARTH_I 0.00005
#define PHYSICS_MOON_I 0.0549
#define PHYSICS_MARS_I 1.85061
#define PHYSICS_JUPITER_I 1.30530
#define PHYSICS_SATURN_I 2.48446
#define PHYSICS_URANUS_I 0.76986
#define PHYSICS_NEPTUNE_I 1.76917

// Aufsteigender Knoten
#define PHYSICS_SUN_NODE 0.0
#define PHYSICS_MERCURY_NODE 48.33167
#define PHYSICS_VENUS_NODE 76.68069
#define PHYSICS_EARTH_NODE -11.26064
#define PHYSICS_MOON_NODE -11.26064
#define PHYSICS_MARS_NODE 49.57854
#define PHYSICS_JUPITER_NODE 100.55615
#define PHYSICS_SATURN_NODE 113.71504
#define PHYSICS_URANUS_NODE 74.22988
#define PHYSICS_NEPTUNE_NODE 131.72169

// Perihellänge
#define PHYSICS_SUN_PL 0.0
#define PHYSICS_MERCURY_PL 77.45645
#define PHYSICS_VENUS_PL 131.53298
#define PHYSICS_EARTH_PL 102.94719
#define PHYSICS_MOON_PL 102.94719
#define PHYSICS_MARS_PL 336.04084
#define PHYSICS_JUPITER_PL 14.75385
#define PHYSICS_SATURN_PL 92.43194
#define PHYSICS_URANUS_PL 170.96424
#define PHYSICS_NEPTUNE_PL 44.97135

#define PHYSICS_SUN_W 0.0
#define PHYSICS_MERCURY_W (PHYSICS_MERCURY_PL - PHYSICS_MERCURY_NODE)
#define PHYSICS_VENUS_W (PHYSICS_VENUS_PL - PHYSICS_VENUS_NODE)
#define PHYSICS_EARTH_W (PHYSICS_EARTH_PL - PHYSICS_EARTH_NODE)
#define PHYSICS_MOON_W (PHYSICS_EARTH_PL - PHYSICS_EARTH_NODE)
#define PHYSICS_MARS_W (PHYSICS_MARS_PL - PHYSICS_MARS_NODE)
#define PHYSICS_JUPITER_W (PHYSICS_JUPITER_PL - PHYSICS_JUPITER_NODE)
#define PHYSICS_SATURN_W (PHYSICS_SATURN_PL - PHYSICS_SATURN_NODE)
#define PHYSICS_URANUS_W (PHYSICS_URANUS_PL - PHYSICS_URANUS_NODE)
#define PHYSICS_NEPTUNE_W (PHYSICS_NEPTUNE_PL - PHYSICS_NEPTUNE_NODE)

double _physics_body_a_AU[] =
    {
        PHYSICS_SUN_A,
        0,
        PHYSICS_MERCURY_A,
        PHYSICS_VENUS_A,
        PHYSICS_EARTH_A,
        PHYSICS_MOON_A,
        PHYSICS_EARTH_A,
        PHYSICS_MARS_A,
        PHYSICS_JUPITER_A,
        PHYSICS_SATURN_A,
        PHYSICS_URANUS_A,
        PHYSICS_NEPTUNE_A
    };

double _physics_body_eccentricity[] =
    {
        PHYSICS_SUN_ECCENTRICITY,
        0,
        PHYSICS_MERCURY_ECCENTRICITY,
        PHYSICS_VENUS_ECCENTRICITY,
        PHYSICS_EARTH_ECCENTRICITY,
        PHYSICS_MOON_ECCENTRICITY,
        PHYSICS_EARTH_ECCENTRICITY,
        PHYSICS_MARS_ECCENTRICITY,
        PHYSICS_JUPITER_ECCENTRICITY,
        PHYSICS_SATURN_ECCENTRICITY,
        PHYSICS_URANUS_ECCENTRICITY,
        PHYSICS_NEPTUNE_ECCENTRICITY
    };

double _physics_body_i_deg[] =
    {
        PHYSICS_SUN_I,
        0,
        PHYSICS_MERCURY_I,
        PHYSICS_VENUS_I,
        PHYSICS_EARTH_I,
        PHYSICS_EARTH_I,
        PHYSICS_EARTH_I,
        PHYSICS_MARS_I,
        PHYSICS_JUPITER_I,
        PHYSICS_SATURN_I,
        PHYSICS_URANUS_I,
        PHYSICS_NEPTUNE_I
    };

double _physics_body_node_deg[] =
    {
        PHYSICS_SUN_NODE,
        0,
        PHYSICS_MERCURY_NODE,
        PHYSICS_VENUS_NODE,
        PHYSICS_EARTH_NODE,
        PHYSICS_EARTH_NODE,
        PHYSICS_EARTH_NODE,
        PHYSICS_MARS_NODE,
        PHYSICS_JUPITER_NODE,
        PHYSICS_SATURN_NODE,
        PHYSICS_URANUS_NODE,
        PHYSICS_NEPTUNE_NODE
    };

double _physics_body_w_deg[] =
    {
        PHYSICS_SUN_W,
        0,
        PHYSICS_MERCURY_W,
        PHYSICS_VENUS_W,
        PHYSICS_EARTH_W,
        PHYSICS_EARTH_W,
        PHYSICS_EARTH_W,
        PHYSICS_MARS_W,
        PHYSICS_JUPITER_W,
        PHYSICS_SATURN_W,
        PHYSICS_URANUS_W,
        PHYSICS_NEPTUNE_W
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
    struct vector_3d r = physics_kepler_r_AU(id, 0.0);
    return vector_norm(&r);
}

double physics_body_aphel_AU(physics_body_id_t id)
{
    struct vector_3d r = physics_kepler_r_AU(id, physics_pi());
    return vector_norm(&r);
}

const char* physics_body_name(physics_body_id_t id)
{
    return _physics_body_names[id];
}

double physics_body_mass_kg(physics_body_id_t id)
{
    return _physics_body_masses_kg[id];
}

double physics_body_a_AU(physics_body_id_t id)
{
    return _physics_body_a_AU[id];
}

double physics_body_eccentricity(physics_body_id_t id)
{
    return _physics_body_eccentricity[id];
}

double physics_body_i_rad(physics_body_id_t id)
{
    return physics_deg_to_rad(_physics_body_i_deg[id]);
}

double physics_body_node_rad(physics_body_id_t id)
{
    return physics_deg_to_rad(_physics_body_node_deg[id]);
}

double physics_body_w_rad(physics_body_id_t id)
{
    return physics_deg_to_rad(_physics_body_w_deg[id]);
}

double physics_body_angle_rad(double t_s, double epoche_rad, double p_AU)
{
    return sqrt(PHYSICS_G * _physics_body_masses_kg[SUN_BARYCENTER] / pow(p_AU * PHYSICS_AU, 3.0)) * t_s + epoche_rad;
}

double physics_pi()
{
    return acos(-1.0);
}

double physics_deg_to_rad(double angle_deg)
{
    return angle_deg * physics_pi() / 180.0;
}

double physics_rad_to_deg(double angle_rad)
{
    return angle_rad * 180.0 / physics_pi();
}

double physics_seconds_per_day()
{
    return 24.0 * 3600.0;
}

double physics_barycenter_AU(double distance_AU, double mass_center_kg, double mass_satellite_kg)
{
    return distance_AU * mass_satellite_kg / (mass_center_kg + mass_satellite_kg);
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

struct vector_3d physics_kepler_r_AU(physics_body_id_t id, double angle_rad)
{
    double e = physics_body_eccentricity(id);
    double R_AU = physics_body_a_AU(id) * (1.0 - e * e) / (e * cos(angle_rad) + 1.0);
    struct vector_3d r_AU = {cos(angle_rad) * R_AU, sin(angle_rad) * R_AU, 0.0};
    return r_AU;
}

struct vector_3d physics_czybor_w_rad_per_sec(physics_body_id_t id, double angle_rad)
{
    double W_rad_per_sec = (sqrt(1.0 - cos(angle_rad)) * sqrt(1.0 + cos(angle_rad)) *
        sqrt(PHYSICS_G * physics_body_mass_kg(SUN_BARYCENTER) * pow(physics_body_eccentricity(id) * cos(angle_rad) + 1.0, 3.0) /
            pow(physics_body_a_AU(id) * PHYSICS_AU, 3.0))) / fabs(sin(angle_rad));
    struct vector_3d w_rad_per_sec = {0.0, 0.0, W_rad_per_sec};
    return w_rad_per_sec;
}
