#include "physics/physics.h"

#include <logging/logging.h>
#include <string/string.h>
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
    cld h = _physics_weber_specific_angular_momentum(body, mass_center_kg);
    cld K = _physics_weber_k(body, mass_center_kg);
    cld v_radial = h * body->e * K * sin(K * phi_rad) / (1.0 - body->e_square);
    return -PHYSICS_G * mass_center_kg * body->mass_kg / vector_norm(&body->r_m) * (1.0L - powl(v_radial, 2.0L) / (2.0L * PHYSICS_C_SQUARE));
}

ld _physics_weber_specific_angular_momentum(const celestial_body_t body, cld mass_center_kg)
{
    return sqrtl(PHYSICS_G * mass_center_kg * body->a_m * (1.0L - body->e_square));
}

ld _physics_weber_k(const celestial_body_t body, cld mass_center_kg)
{
    cld A = 6.0 * PHYSICS_G * mass_center_kg;
    cld B = PHYSICS_C_SQUARE * body->a_m * (1.0 - body->e_square);
    cld C = sqrt(1.0 - A / B);
    return C;
}

struct vector_3d physics_weber_position(const celestial_body_t body, cld mass_center_kg, cld phi_rad)
{
    struct vector_3d position = (struct vector_3d)
    {
        .x = cosl(phi_rad),
        .y = sinl(phi_rad),
        .z = 0.0
    };

    cld h = _physics_weber_specific_angular_momentum(body, mass_center_kg);
    cld K = _physics_weber_k(body, mass_center_kg);
    cld A = body->a_m * (1.0L - body->e_square) / (1.0L + body->e * cosl(K * phi_rad));
    cld B = 3.0L * powl(PHYSICS_G, 2.0) * powl(mass_center_kg, 2.0L) / (PHYSICS_C_SQUARE * powl(h, 4.0L));
    cld C = A * (1.0L + B * (1.0L + body->e_square / 2.0L + body->e * phi_rad * sinl(K * phi_rad)));
    
    return vector_multiply_scalar(&position, C);
}

ld physics_weber_periodtime(const celestial_body_t body, cld mass_center_kg)
{
    cld A = 2.0L * physics_pi() * sqrtl(powl(body->a_m, 3.0L) / (PHYSICS_G * mass_center_kg));
    cld B = 3.0L * PHYSICS_G * mass_center_kg / (2.0L * PHYSICS_C_SQUARE * body->a_m * (1.0L - body->e_square));
    return A * (1.0L + B);
}

ld physics_weber_deltaphi(const celestial_body_t body, cld mass_center_kg, cld T_step_s)
{
    cld A = vector_norm(&body->w_rad_s) * T_step_s;
    cld h = _physics_weber_specific_angular_momentum(body, mass_center_kg);
    cld r = vector_norm(&body->r_m);
    cld B = 3.0L * PHYSICS_G * mass_center_kg * h * T_step_s / (PHYSICS_C_SQUARE * powl(r, 3.0L));
    return A + B;
}

struct vector_3d physics_weber_angular_speed(const celestial_body_t body, cld mass_center_kg, cld phi_rad)
{
    struct vector_3d w = {0};
    cld h = _physics_weber_specific_angular_momentum(body, mass_center_kg);
    struct vector_3d r = physics_weber_position(body, mass_center_kg, phi_rad);
    w.z = h / vector_dot(&r, &r);
    return w;
}
