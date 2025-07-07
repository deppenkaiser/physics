#include "physics/physics.h"

#include <logging/logging.h>
#include <string/string.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <float.h>

ld _physics_weber_k(const celestial_body_t body, cld mass_center_kg);
ld _physics_weber_specific_angular_momentum(const celestial_body_t body, cld mass_center_kg);

ld physics_frac(cld x)
{
    return x - floorl(x);
}

ld physics_modulo(cld a, cld b)
{
    return b * physics_frac(a / b);
}

ld physics_pi()
{
    return acos(-1.0);
}

ld physics_deg_to_rad(cld angle_deg)
{
    return angle_deg * physics_pi() / 180.0;
}

ld physics_rad_to_deg(cld angle_rad)
{
    return angle_rad * 180.0 / physics_pi();
}

ld physics_seconds_per_day()
{
    return 24.0L * 3600.0L;
}

ld physics_barycenter_AU(cld distance_AU, cld mass_center_kg, cld mass_satellite_kg)
{
    return distance_AU * mass_satellite_kg / (mass_center_kg + mass_satellite_kg);
}

ld physics_rayleigh_criteria_deg(cld wavelength_light_m, cld objective_aperture_m)
{
    return physics_rad_to_deg(asinl(1.22L * wavelength_light_m / objective_aperture_m));
}

ld physics_image_sensor_object_size_m(cld object_size_deg, cld focal_length_m)
{
    return 2.0L * focal_length_m * tanl(physics_deg_to_rad(object_size_deg) / 2.0L);
}

ld physics_needed_image_sensor_pixel_size_m(cld wavelength_light_m, cld objective_aperture_m, cld focal_length_m)
{
    return physics_image_sensor_object_size_m(physics_rayleigh_criteria_deg(wavelength_light_m, objective_aperture_m), focal_length_m) / 2.0L;
}

ld physics_kepler_radius(cld a_m, cld eccentricity, cld phi_rad)
{
    return a_m * (1.0L - powl(eccentricity, 2.0L)) / (1.0L + eccentricity * cosl(phi_rad));
}

ld physics_kinetic_energy(cld mass_kg, const vector_3d_t v)
{
    return 0.5L * mass_kg * vector_dot(v, v);
}

ld physics_kinetic_energy_body(const celestial_body_t body)
{
    return physics_kinetic_energy(body->mass_kg, &body->v_m_s);
}

ld physics_weber_potential_energy(const celestial_body_t body, cld mass_center_kg, cld phi_rad)
{
    cld v_radial = body->h * body->e * body->K_2 * sin(body->K_2 * phi_rad) / (1.0 - body->e_square);
    return -PHYSICS_G * mass_center_kg * body->mass_kg / vector_norm(&body->r_m) * (1.0L - powl(v_radial, 2.0L) / (2.0L * PHYSICS_C_SQUARE));
}

ld physics_weber_h(const celestial_body_t body, cld mass_center_kg)
{
    return sqrtl(PHYSICS_G * mass_center_kg * body->a_m * (1.0L - body->e_square));
}

ld physics_weber_alpha(const celestial_body_t body, cld mass_center_kg)
{
    return 3.0L * powl(PHYSICS_G, 2.0L) * powl(mass_center_kg, 2.0L) * body->e / (8.0L * powl(body->h, 4.0L) * powl(PHYSICS_C, 4.0L));
}

ld physics_weber_k_1st_order(const celestial_body_t body, cld mass_center_kg)
{
    cld A = 6.0 * PHYSICS_G * mass_center_kg;
    cld B = PHYSICS_C_SQUARE * body->a_m * (1.0 - body->e_square);
    cld C = sqrt(1.0 - A / B);
    return C;
}

ld physics_weber_k_2nd_order(const celestial_body_t body, cld mass_center_kg)
{
    cld A = 6.0 * PHYSICS_G * mass_center_kg;
    cld B = PHYSICS_C_SQUARE * body->a_m * (1.0 - body->e_square);
    cld C = 27.0L * powl(PHYSICS_G, 2.0L) * powl(mass_center_kg, 2.0L);
    cld D = 2.0L * powl(PHYSICS_C, 4.0L) * powl(body->a_m, 2.0L) * powl(1.0 - body->e_square, 2.0L);
    cld E = sqrt(1.0 - A / B + C / D);
    return E;
}

struct vector_3d physics_weber_position(const celestial_body_t body, cld mass_center_kg, cld phi_rad)
{
    struct vector_3d position = (struct vector_3d)
    {
        .x = cosl(phi_rad),
        .y = sinl(phi_rad),
        .z = 0.0
    };

    cld B = body->a_m * (1.0L - body->e_square) / (1.0L + body->e * cosl(body->K_2 * phi_rad + body->A * powl(phi_rad, 2.0L)));
    return vector_multiply_scalar(&position, B);
}

struct vector_3d physics_weber_angular_speed(const celestial_body_t body, cld mass_center_kg, cld phi_rad)
{
    struct vector_3d w = {0};
    cld B = body->h * powl(1.0L + body->e * cosl(body->K_2 * phi_rad + body->A * powl(phi_rad, 2.0L)), 2.0L);
    w.z = B / (powl(body->a_m, 2.0L) * powl(1.0L - body->e_square, 2.0L));
    return w;
}

ld physics_weber_periodtime(const celestial_body_t body, cld mass_center_kg)
{
    cld A = 2.0L * physics_pi() * powl(body->a_m, 3.0L / 2.0L) / (sqrtl(PHYSICS_G * mass_center_kg));
    cld B = 3.0L * PHYSICS_G * mass_center_kg / (2.0L * PHYSICS_C_SQUARE * body->a_m * (1.0L - body->e_square));
    cld C = 45.0L * powl(PHYSICS_G, 2.0L) * powl(mass_center_kg, 2.0L) /
        (8.0L * powl(PHYSICS_C, 4.0L) * powl(body->a_m, 2.0L) * powl(1.0L - body->e_square, 2.0L));
    return A * (1.0L + B + C * (1.0L - body->e_square / 3.0L));
}

ld physics_weber_deltaphi(const celestial_body_t body, cld mass_center_kg, cld t_step_s, cld phi_0_rad)
{
    cld phi = physics_weber_angular_speed(body, mass_center_kg, phi_0_rad).z * t_step_s;
    return phi;
}

ld physics_deltaphi_per_revolution_1(const celestial_body_t body, cld mass_center_kg)
{
    return 2.0L * physics_pi() * (1.0L / body->K_1 - 1.0L);
}

ld physics_deltaphi_per_revolution_2(const celestial_body_t body, cld mass_center_kg)
{
    return 2.0L * physics_pi() * (1.0L / body->K_2 - 1.0L) - 4.0L * powl(physics_pi(), 2.0L) * body->A / powl(body->K_2, 3.0L);
}
