#pragma once

#include "ray.h"

// Physics simulation namespace (geometric units: c = G = 1, RS = 1)
namespace Physics {
    // Compute radial and angular accelerations from Schwarzschild geodesic equations
    // a_r and a_phi are filled with the computed accelerations
    void computeAccelerations(float r, float v_r, float v_phi, float E_val,
                              float& a_r, float& a_phi);

    // Velocity Verlet integration step (symplectic integrator)
    // Uses "kick-drift-kick" scheme for better energy conservation
    void verletStep(Ray& ray, float dlambda);

    // Adaptive timestep based on distance from black hole
    // Returns smaller step sizes near the black hole (r close to 1) for accuracy
    // Returns larger step sizes far away (r > 10) for performance
    float adaptiveStep(float r);

    // Integrate ray by a fixed "budget" using adaptive substeps
    // This maintains constant visual speed while using adaptive accuracy
    void integrateWithBudget(Ray& ray, float budget, float maxDistance, int currentFrame);
}

