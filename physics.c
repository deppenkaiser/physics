#include "physics.h"

#include <logging/logging.h>
#include <string/string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <float.h>

// https://www.imcce.fr/

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

double physics_frac(double x)
{
    return x - floor(x);
}

double physics_modulo(double a, double b)
{
    return b * physics_frac(a / b);
}

double physics_czybor_angular_speed_rad_per_s(physics_body_id_t planet_id, physics_body_id_t center_id, double phi_rad)
{
    double e = physics_body_eccentricity(planet_id);
    double M_kg = physics_body_mass_kg(center_id);
    double a_m = physics_body_a_AU(planet_id) * PHYSICS_AU;
    double angular_speed_rad_per_s = 0.0;

    phi_rad = physics_modulo(phi_rad, 2.0 * physics_pi());
    double cos_phi = cos(phi_rad);

    if (phi_rad != 0.0)
    {
        angular_speed_rad_per_s = sqrt(1.0 - cos_phi) * sqrt(1.0 + cos_phi) *
            sqrt(-(PHYSICS_G * M_kg * pow(e * cos_phi + 1.0, 3.0)) /
                (pow(a_m, 3.0) * pow(e - 1.0, 3.0) * pow(e + 1.0, 3.0))) / fabs(sin(phi_rad));

        bool is_zero = angular_speed_rad_per_s == 0.0;

        if (is_zero)
        {
            double angular_speed_top_rad_per_sec = physics_czybor_angular_speed_rad_per_s(planet_id, center_id, phi_rad + FLT_EPSILON);
            double angular_speed_bottom_rad_per_sec = physics_czybor_angular_speed_rad_per_s(planet_id, center_id, phi_rad - FLT_EPSILON);
            angular_speed_rad_per_s = (angular_speed_top_rad_per_sec + angular_speed_bottom_rad_per_sec) / 2.0;
        }
    }
    else
    {
        angular_speed_rad_per_s = sqrt(2.0) * sqrt(-(PHYSICS_G * M_kg / (pow(a_m, 3.0) * pow(e - 1.0, 3.0))));
    }

    return angular_speed_rad_per_s;
}

double physics_czybor_mean_angular_speed_rad_per_s(physics_body_id_t planet_id, physics_body_id_t center_id)
{
    double e = physics_body_eccentricity(planet_id);
    double M_kg = physics_body_mass_kg(center_id);
    double a_m = physics_body_a_AU(planet_id) * PHYSICS_AU;
    return sqrt(-(PHYSICS_G * M_kg) / (pow(a_m, 3.0) * pow(e - 1.0, 3.0) * pow(e + 1.0, 3.0)));
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

double physics_orbital_period_s(physics_body_id_t planet_id, physics_body_id_t center_id,
    double phi_1_rad, double phi_2_rad, uint32_t step_count)
{
    double dx_rad = (phi_2_rad - phi_1_rad) / step_count;
    double sum_s = 0.0;
    for (uint32_t i = 0; i < step_count; ++i)
    {
        double phi_rad = phi_1_rad + i * dx_rad;
        sum_s += 1.0 / physics_czybor_angular_speed_rad_per_s(planet_id, center_id, phi_rad);
    }
    return sum_s * dx_rad;
}

double physics_optimized_orbital_period_s(physics_body_id_t planet_id, physics_body_id_t center_id,
    double phi_1_rad, double phi_2_rad)
{
    double last_period_s = 0.0;
    double period_s = 0.0;
    uint32_t i = 0;

    while (phi_2_rad > phi_1_rad)
    {
        uint32_t step_count = (uint32_t) pow(10.0, (double) i);
        
        period_s = physics_orbital_period_s(planet_id, center_id, phi_1_rad, phi_2_rad, step_count);
    
        if (fabs(period_s - last_period_s) < (period_s * 0.00001))
        {
            string_t message = {0};
            sprintf(message, "orbital period: %s/%s %f s in 10^%d steps; iteration instability: %f %%", physics_body_name(center_id),
                physics_body_name(planet_id), last_period_s, i, fabs(100.0 - last_period_s / period_s * 100.0));
            logging_log_message(message, true);
            break;
        }

        last_period_s = period_s;
        ++i;
    }

    return period_s;
}

double physics_full_orbital_period_days(physics_body_id_t planet_id, physics_body_id_t center_id)
{
    return physics_optimized_orbital_period_s(planet_id, center_id, 0.0, 2.0 * physics_pi()) / physics_seconds_per_day();
}

double physics_rayleigh_criteria_deg(double wavelength_light_m, double objective_aperture_m)
{
    return physics_rad_to_deg(asin(1.22 * wavelength_light_m / objective_aperture_m));
}

double physics_image_sensor_object_size_m(double object_size_deg, double focal_length_m)
{
    return 2.0 * focal_length_m * tan(physics_deg_to_rad(object_size_deg) / 2.0);
}

struct vector_3d physics_kepler_r_AU(physics_body_id_t id, double angle_rad)
{
    double e = physics_body_eccentricity(id);
    double a_AU = physics_body_a_AU(id);
    double R_AU = a_AU * (1.0 - e * e) / (e * cos(angle_rad) + 1.0);
    struct vector_3d r_AU = {cos(angle_rad) * R_AU, sin(angle_rad) * R_AU, 0.0};
    return r_AU;
}
