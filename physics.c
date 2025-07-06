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
    cld h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    cld K = _physics_weber_k(body, mass_center_kg);
    cld v_radial = h * body->e * K * sin(K * phi_rad) / (1.0 - body->e_square);
    return -PHYSICS_G * mass_center_kg * body->mass_kg / vector_norm(&body->r_m) * (1.0L - powl(v_radial, 2.0L) / (2.0L * PHYSICS_C_SQUARE));
}

ld physics_weber_specific_angular_momentum(const celestial_body_t body, cld mass_center_kg)
{
    return sqrtl(PHYSICS_G * mass_center_kg * body->a_m * (1.0L - body->e_square));
}

ld _physics_weber_alpha(const celestial_body_t body, cld mass_center_kg)
{
    cld h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    return 3.0L * powl(PHYSICS_G, 2.0L) * powl(mass_center_kg, 2.0L) * body->e / (8.0L * powl(h, 4.0L));
}

ld _physics_weber_k(const celestial_body_t body, cld mass_center_kg)
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

    cld h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    cld K = _physics_weber_k(body, mass_center_kg);
    cld A = _physics_weber_alpha(body, mass_center_kg);
    cld B = body->a_m * (1.0L - body->e_square) / (1.0L + body->e * cosl(K * phi_rad + A * powl(phi_rad, 2.0L) / powl(PHYSICS_C, 4.0L)));
    return vector_multiply_scalar(&position, B);
}

struct vector_3d physics_weber_angular_speed(const celestial_body_t body, cld mass_center_kg, cld phi_rad)
{
    struct vector_3d w = {0};
    cld h = physics_weber_specific_angular_momentum(body, mass_center_kg);
    cld K = _physics_weber_k(body, mass_center_kg);
    cld A = _physics_weber_alpha(body, mass_center_kg);
    cld B = h * powl(1.0L + body->e * cosl(K * phi_rad + A * powl(phi_rad, 2.0L)), 2.0L);
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

typedef struct rk4_params
{
    celestial_body_t body;
    cld mass_center_kg;
} *rk4_params_t;

int func(double t, const double phi[], double dphidt[], void *p)
{
    rk4_params_t params = (rk4_params_t) p;
    dphidt[0] = physics_weber_angular_speed(params->body, params->mass_center_kg, phi[0]).z;
    return GSL_SUCCESS;
}

ld physics_weber_deltaphi(const celestial_body_t body, cld mass_center_kg, cld t_0_s, cld t_1_s, cld phi_0_rad)
{
    struct rk4_params p = (struct rk4_params)
    {
        .body = body,
        .mass_center_kg = mass_center_kg
    };

    gsl_odeiv2_system sys = {func, NULL, 1, &p}; // 1D-System (phi)
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0); // RK4 mit tolerierbarem Fehler
    double phi = phi_0_rad; // Startwinkel (z. B. Perihel bei phi=0)
    double t = t_0_s, t_end = t_1_s; // Integrationszeitraum
    gsl_odeiv2_driver_apply(driver, &t, t_end, &phi);
    gsl_odeiv2_driver_free(driver);

    return phi - phi_0_rad;
}
