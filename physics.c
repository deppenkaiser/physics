#include "physics/physics.h"

#include <logging/logging.h>
#include <string/string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <float.h>

double physics_frac(double x)
{
    return x - floor(x);
}

double physics_modulo(double a, double b)
{
    return b * physics_frac(a / b);
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

double physics_rayleigh_criteria_deg(double wavelength_light_m, double objective_aperture_m)
{
    return physics_rad_to_deg(asin(1.22 * wavelength_light_m / objective_aperture_m));
}

double physics_image_sensor_object_size_m(double object_size_deg, double focal_length_m)
{
    return 2.0 * focal_length_m * tan(physics_deg_to_rad(object_size_deg) / 2.0);
}

double physics_needed_image_sensor_pixel_size_m(double wavelength_light_m, double objective_aperture_m, double focal_length_m)
{
    return physics_image_sensor_object_size_m(physics_rayleigh_criteria_deg(wavelength_light_m, objective_aperture_m), focal_length_m) / 2.0;
}

double physics_specific_angular_momentum(double mass_center_kg, double a_m, double eccentricity)
{
    return sqrt(PHYSICS_G * mass_center_kg * a_m * (1.0 - pow(eccentricity, 2.0)));
}

double physics_angular_speed(double angular_moment, double mass_center_kg, double eccentricity, double phi_rad)
{
    double sin_phi = sin(phi_rad);
    double cos_phi = cos(phi_rad);
    double e = 1.0 + eccentricity * cos_phi;
    double factor = 1.0 + pow(eccentricity, 2.0) / 2.0 + eccentricity * phi_rad * sin_phi;
    return pow(angular_moment, 3.0) / (pow(PHYSICS_G, 2.0) * pow(mass_center_kg, 2.0) * pow(e, 2.0)) +
        3.0 * angular_moment / (PHYSICS_C_SQUARE * e) * factor;
}

double physics_kepler_radius(double a_m, double eccentricity, double phi_rad)
{
    return a_m * (1 - pow(eccentricity, 2.0)) / (1.0 + eccentricity * cos(phi_rad));
}

struct vector_3d compute_weber_force(const celestial_body_t body1, const celestial_body_t body2) 
{
    struct vector_3d force = {0};

    // 1. Relativvektoren berechnen
    struct vector_3d r_rel = vector_sub(&body2->position, &body1->position);
    struct vector_3d v_rel = vector_sub(&body2->velocity, &body1->velocity);
    struct vector_3d a_rel = vector_sub(&body2->acceleration, &body1->acceleration);

    // 2. Abstand und Einheitsvektor
    double r = vector_norm(&r_rel);
    struct vector_3d r_rel_unit = vector_divide_scalar(&r_rel, r);

    // 3. Weber-Kraft
    double v_square = vector_dot(&v_rel, &v_rel);
    double ra = vector_dot(&r_rel, &a_rel);
    force = vector_multiply_scalar(&r_rel_unit, -PHYSICS_G * body1->mass * body2->mass / pow(r, 2.0) *
        (1.0 - v_square / PHYSICS_C_SQUARE + ra / (2.0 * PHYSICS_C_SQUARE)));

    return force;
}
