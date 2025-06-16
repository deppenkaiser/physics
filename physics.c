#include "physics/physics.h"

#include <logging/logging.h>
#include <string/string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <float.h>

double _physics_weber_k(const celestial_body_t body, cd mass_center_kg);

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
    return physics_kinetic_energy(body->mass_kg, &body->v_m_s);
}

double physics_weber_potential_energy(const celestial_body_t body, cd mass_center_kg, cd phi_rad)
{
    double h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    double K = _physics_weber_k(body, mass_center_kg);
    double v_radial = h * body->e * K * sin(K * phi_rad) / (1.0 - pow(body->e, 2.0));
    return -PHYSICS_G * mass_center_kg * body->mass_kg / vector_norm(&body->r_m) * (1.0 - pow(v_radial, 2.0) / (2.0 * PHYSICS_C_SQUARE));
}

double _physics_weber_k(const celestial_body_t body, cd mass_center_kg)
{
    double A = 6.0 * PHYSICS_G * mass_center_kg;
    double B = PHYSICS_C_SQUARE * body->a_m * (1.0 - pow(body->e, 2.0));
    double C = sqrt(1.0 - A / B);
    return C;
}

double physics_weber_specific_angular_momentum(const celestial_body_t body, cd mass_center_kg)
{
    return sqrt(PHYSICS_G * mass_center_kg * body->a_m * (1.0 - pow(body->e, 2.0)));
}

struct vector_3d physics_weber_angular_speed(const celestial_body_t body, cd mass_center_kg, cd phi_rad)
{
    struct vector_3d w = {0};
    double h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    struct vector_3d r = physics_weber_position(body, mass_center_kg, phi_rad);
    w.z = h / vector_dot(&r, &r);
    return w;
}

struct vector_3d physics_weber_angular_speed_perturbed(const celestial_body_t body, const celestial_body_t body_i, cd mass_center_kg)
{
    struct vector_3d dr = physics_weber_position_perturbed(body, body_i, mass_center_kg);
    struct vector_3d dv = physics_weber_velocity_perturbed(body, body_i, mass_center_kg);
    struct vector_3d A = vector_cross(&body->r_m, &dv);
    struct vector_3d B = vector_cross(&dr, &body->v_m_s);
    struct vector_3d w = vector_add(&A, &B);
    return vector_divide_scalar(&w, vector_dot(&body->r_m, &body->r_m));
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

struct vector_3d physics_weber_position_perturbed(const celestial_body_t body, const celestial_body_t body_i, cd mass_center_kg)
{
    struct vector_3d r = vector_sub(&body->r_m, &body_i->r_m);
    double A = body_i->mass_kg / mass_center_kg * pow(body_i->a_m, 2.0) / vector_dot(&r, &r);
    return vector_multiply_scalar(&r, A);
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

struct vector_3d physics_weber_velocity_perturbed(const celestial_body_t body, const celestial_body_t body_i, cd mass_center_kg)
{
    struct vector_3d vr = vector_sub(&body->r_m, &body_i->r_m);
    struct vector_3d z = (struct vector_3d)
    {
        .x = 0.0,
        .y = 0.0,
        .z = 1.0
    };
    vr = vector_cross(&vr, &z);
    double r = vector_norm(&vr);
    double h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    double A = PHYSICS_G * body_i->mass_kg / (h * pow(r, 3.0));
    return vector_multiply_scalar(&vr, A);
}
