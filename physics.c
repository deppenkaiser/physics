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
    // 1. Relativvektoren berechnen
    struct vector_3d r_vec = vector_sub(&body2->position, &body1->position);
    struct vector_3d v_rel = vector_sub(&body2->velocity, &body1->velocity);
    struct vector_3d a_rel = vector_sub(&body2->acceleration, &body1->acceleration);

    // 2. Abstand und Einheitsvektor
    double r = vector_norm(&r_vec);
    struct vector_3d r_unit = vector_divide_scalar(&r_vec, r);

    // 3. Radialkomponenten von v und a
    double v_r = vector_dot(&v_rel, &r_unit);
    
    // Tangentialgeschwindigkeit (v_tan = v - (v·r)r)
    struct vector_3d v_tan =
    {
        .x = v_rel.x - v_r * r_unit.x,
        .y = v_rel.y - v_r * r_unit.y,
        .z = v_rel.z - v_r * r_unit.z
    };
    double v_tan_sq = vector_dot(&v_tan, &v_tan);
    
    // Radialbeschleunigung (a_r = a·r + v_tan²/r)
    double a_r = vector_dot(&a_rel, &r_unit) + v_tan_sq / r;

    // 4. Weber-Kraft berechnen (β = 0.5 für Gravitation)
    double F_magnitude = -(PHYSICS_G * body1->mass * body2->mass) / (r * r);
    double weber_factor = 1.0 - (v_r*v_r)/(PHYSICS_C_SQUARE) + (r*a_r)/(2.0 * PHYSICS_C_SQUARE);

    // 5. Kraftvektor
    struct vector_3d force =
    {
        .x = F_magnitude * weber_factor * r_unit.x,
        .y = F_magnitude * weber_factor * r_unit.y,
        .z = F_magnitude * weber_factor * r_unit.z
    };

    return force;
}
