#include "physics.h"
#include "constants.h"
#include <cmath>

namespace Physics {

// Compute radial and angular accelerations from Schwarzschild geodesic equations
// Geodesic equations in geometric units (RS = 1)
void computeAccelerations(float r, float v_r, float v_phi, float E_val,
                          float& a_r, float& a_phi) {
    // f = 1 - 1/r (geometric units: RS = 1)
    const float f{1.0f - 1.0f / r};
    const float dt_dlambda{E_val / f};

    // a_r = dv_r/dλ (radial acceleration from geodesic equation)
    // Exact Schwarzschild null geodesic equation (geometric units: RS = 1):
    // dv_r/dλ = -(1/(2r²))f(dt/dλ)² + (1/(2r²f))(v_r)² + (r - 1)(v_φ)²
    a_r = -(1.0f / (2.0f * r * r)) * f * (dt_dlambda * dt_dlambda)
          + (1.0f / (2.0f * r * r * f)) * (v_r * v_r)
          + (r - 1.0f) * (v_phi * v_phi);

    // a_phi = dv_phi/dλ (angular acceleration from geodesic equation)
    // dv_φ/dλ = -2v_r·v_φ/r
    a_phi = -2.0f * v_r * v_phi / r;
}

// Velocity Verlet integration step (symplectic integrator)
// Uses "kick-drift-kick" scheme for better energy conservation
void verletStep(Ray& ray, float dlambda) {
    const float half_dt{dlambda / 2.0f};

    // Compute initial accelerations
    float a_r, a_phi;
    computeAccelerations(ray.r, ray.v_r, ray.v_phi, ray.E, a_r, a_phi);

    // Half-kick: update velocities by half step
    ray.v_r += a_r * half_dt;
    ray.v_phi += a_phi * half_dt;

    // Drift: update positions using half-step velocities
    ray.r += ray.v_r * dlambda;
    ray.phi += ray.v_phi * dlambda;

    // Compute new accelerations at updated position
    computeAccelerations(ray.r, ray.v_r, ray.v_phi, ray.E, a_r, a_phi);

    // Half-kick: update velocities by another half step
    ray.v_r += a_r * half_dt;
    ray.v_phi += a_phi * half_dt;
}

// Adaptive timestep: smaller steps near black hole for accuracy
// Returns smaller step sizes near the black hole (r close to 1) for accuracy
// Returns larger step sizes far away (r > 10) for performance
float adaptiveStep(float r) {
    // Scale linearly based on (r - 1.0f) / 10.0f, clamped to [0, 1]
    float factor{std::min((r - 1.0f) / 10.0f, 1.0f)};
    factor = std::max(factor, 0.0f);  // Clamp to non-negative
    return Simulation::MIN_STEP + factor * (Simulation::MAX_STEP - Simulation::MIN_STEP);
}

// Integrate ray by a fixed "budget" using adaptive substeps
// This maintains constant visual speed while using adaptive accuracy
void integrateWithBudget(Ray& ray, float budget, float maxDistance, int currentFrame) {
    if (!ray.isActive(currentFrame)) {
        return;
    }

    if (ray.isCaptured() || ray.hasEscaped(maxDistance)) {
        return;
    }

    float accumulated{0.0f};
    while (accumulated < budget) {
        float step{adaptiveStep(ray.r)};
        step = std::min(step, budget - accumulated);  // Don't overshoot budget

        verletStep(ray, step);
        accumulated += step;

        if (ray.isCaptured() || ray.hasEscaped(maxDistance)) {
            break;
        }
    }
    ray.recordPosition();
    ray.updateDeflection();
}

} // namespace Physics