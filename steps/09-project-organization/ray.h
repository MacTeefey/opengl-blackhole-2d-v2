#pragma once

#include <glm/glm.hpp>
#include <vector>

#error "TODO: Create ray.h with RayScenario enum and Ray struct - remove this line when done"

// TODO: Create RayScenario enum class with PARALLEL, POINT_SOURCE, ORBITING

// TODO: Create Ray struct with:
// Member variables (all float for geometric units):
// - r, phi (polar coordinates)
// - v_r, v_phi (velocities)
// - E, L (conserved quantities)
// - trail (std::vector<glm::vec2>)
// - initialVelocityAngle, deflection
// - scenario (RayScenario), startFrame (int)

// Methods to declare:
// - Constructor: Ray(float x, float y, float vx, float vy, RayScenario scenario, int startFrame = 0)
// - bool isCaptured() const
// - bool hasEscaped(float maxDistance) const
// - void recordPosition()
// - void updateDeflection()
// - bool isActive(int currentFrame) const
