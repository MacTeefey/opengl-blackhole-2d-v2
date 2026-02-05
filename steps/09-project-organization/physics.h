#pragma once

#include "ray.h"

#error "TODO: Create Physics namespace with function declarations - remove this line when done"

// TODO: Create Physics namespace (geometric units: c = G = 1, RS = 1)

// TODO: Declare computeAccelerations function
// - Takes float r, v_r, v_phi, E_val, and output refs a_r, a_phi
// - Computes radial and angular accelerations from Schwarzschild geodesic equations

// TODO: Declare verletStep function
// - Takes Ray& and float dlambda
// - Velocity Verlet integration using kick-drift-kick scheme

// TODO: Declare adaptiveStep function
// - Takes float r
// - Returns float step size (smaller near black hole, larger far away)

// TODO: Declare integrateWithBudget function
// - Takes Ray&, float budget, float maxDistance, int currentFrame
// - Integrates ray using adaptive substeps within fixed budget
