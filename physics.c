#include "physics/physics.h"

#include <logging/logging.h>
#include <string/string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <float.h>

double physics_frac(cd x)
{
    return x - floor(x);
}

double physics_modulo(cd a, cd b)
{
    return b * physics_frac(a / b);
}

double physics_pi()
{
    return acos(-1.0);
}

double physics_deg_to_rad(cd angle_deg)
{
    return angle_deg * physics_pi() / 180.0;
}

double physics_rad_to_deg(cd angle_rad)
{
    return angle_rad * 180.0 / physics_pi();
}

double physics_seconds_per_day()
{
    return 24.0 * 3600.0;
}

double physics_barycenter_AU(cd distance_AU, cd mass_center_kg, cd mass_satellite_kg)
{
    return distance_AU * mass_satellite_kg / (mass_center_kg + mass_satellite_kg);
}

double physics_rayleigh_criteria_deg(cd wavelength_light_m, cd objective_aperture_m)
{
    return physics_rad_to_deg(asin(1.22 * wavelength_light_m / objective_aperture_m));
}

double physics_image_sensor_object_size_m(cd object_size_deg, cd focal_length_m)
{
    return 2.0 * focal_length_m * tan(physics_deg_to_rad(object_size_deg) / 2.0);
}

double physics_needed_image_sensor_pixel_size_m(cd wavelength_light_m, cd objective_aperture_m, cd focal_length_m)
{
    return physics_image_sensor_object_size_m(physics_rayleigh_criteria_deg(wavelength_light_m, objective_aperture_m), focal_length_m) / 2.0;
}

double physics_kepler_radius(cd a_m, cd eccentricity, cd phi_rad)
{
    return a_m * (1 - pow(eccentricity, 2.0)) / (1.0 + eccentricity * cos(phi_rad));
}

double physics_kinetic_energy(cd mass_kg, const vector_3d_t v)
{
    return 0.5 * mass_kg * vector_dot(v, v);
}

double physics_kinetic_energy_body(const celestial_body_t body)
{
    return physics_kinetic_energy(body->mass_kg, &body->velocity);
}

double _physics_weber_k(const celestial_body_t body, cd mass_center_kg)
{
    double A = 6.0 * PHYSICS_G * mass_center_kg;
    double B = PHYSICS_C_SQUARE * body->a_m * (1.0 - pow(body->e, 2.0));
    double C = sqrt(1.0 - A / B);
    return C;
}

double physics_weber_potential_energy(const celestial_body_t body, cd mass_center_kg, cd phi_rad)
{
    double h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    double K = _physics_weber_k(body, mass_center_kg);
    double v_radial = h * body->e * K * sin(K * phi_rad) / (1.0 - pow(body->e, 2.0));
    return -PHYSICS_G * mass_center_kg * body->mass_kg / vector_norm(&body->position) * (1.0 - pow(v_radial, 2.0) / (2.0 * PHYSICS_C_SQUARE));
}

struct vector_3d physics_weber_force(const celestial_body_t body1, const celestial_body_t body2)
{
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
    
    return vector_multiply_scalar(&r_rel_unit, -PHYSICS_G * body1->mass_kg * body2->mass_kg / pow(r, 2.0) *
        (1.0 - v_square / PHYSICS_C_SQUARE + ra / (2.0 * PHYSICS_C_SQUARE)));
}

double physics_weber_specific_angular_momentum(const celestial_body_t body, cd mass_center_kg)
{
    return sqrt(PHYSICS_G * mass_center_kg * body->a_m * (1.0 - pow(body->e, 2.0)));
}

double physics_weber_angular_speed(cd angular_moment, cd mass_center_kg, cd eccentricity, cd phi_rad)
{
    double sin_phi = sin(phi_rad);
    double cos_phi = cos(phi_rad);
    double e = 1.0 + eccentricity * cos_phi;
    double factor = 1.0 + pow(eccentricity, 2.0) / 2.0 + eccentricity * phi_rad * sin_phi;
    return pow(angular_moment, 3.0) / (pow(PHYSICS_G, 2.0) * pow(mass_center_kg, 2.0) * pow(e, 2.0)) +
        3.0 * angular_moment / (PHYSICS_C_SQUARE * e) * factor;
}

struct vector_3d physics_weber_position(const celestial_body_t body, cd mass_center_kg, cd phi_rad)
{
    struct vector_3d position = (struct vector_3d)
    {
        .x = cos(phi_rad),
        .y = sin(phi_rad),
        .z = 0.0
    };

    double h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    double K = _physics_weber_k(body, mass_center_kg);
    double A = body->a_m * (1.0 - pow(body->e, 2.0)) / (1.0 + body->e * cos(K * phi_rad));
    double B = 3.0 * pow(PHYSICS_G, 2.0) * pow(mass_center_kg, 2.0) / (PHYSICS_C_SQUARE * pow(h, 4.0));
    double C = A * (1.0 + B * (1.0 + pow(body->e, 2.0) / 2.0 + body->e * phi_rad * sin(K * phi_rad)));
    
    return vector_multiply_scalar(&position, C);
}

struct vector_3d physics_weber_velocity(const celestial_body_t body, cd mass_center_kg, cd phi_rad)
{
    struct vector_3d velocity_1 = (struct vector_3d)
    {
        .x = cos(phi_rad),
        .y = sin(phi_rad),
        .z = 0.0
    };

    struct vector_3d velocity_2 = (struct vector_3d)
    {
        .x = -sin(phi_rad),
        .y = cos(phi_rad),
        .z = 0.0
    };

    double K = _physics_weber_k(body, mass_center_kg);
    double A = sqrt(PHYSICS_G * mass_center_kg / (body->a_m * (1 - pow(body->e, 2.0))));
    double B = body->e * K * sin(K * phi_rad) / (1.0 + body->e * cos(K * phi_rad));
    double C = 1.0 + body->e * cos(K * phi_rad);

    velocity_1 = vector_multiply_scalar(&velocity_1, A * B);
    velocity_2 = vector_multiply_scalar(&velocity_2, A * C);
    
    return vector_add(&velocity_1, &velocity_2);
}

struct vector_3d physics_weber_acceleration(const celestial_body_t body, cd mass_center_kg, cd phi_rad)
{
    struct vector_3d acceleration = {0};
    return acceleration;
}
