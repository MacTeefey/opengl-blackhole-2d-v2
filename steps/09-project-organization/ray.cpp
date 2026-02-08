#include "ray.h"
#include "constants.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Ray::Ray(float x, float y, float vx, float vy, RayScenario scenario, int startFrame)
    : scenario{scenario}, startFrame{startFrame}, deflection{0.0f} {
    // Convert to polar coordinates
    r = std::sqrt(x * x + y * y);
    phi = std::atan2(y, x);

    // Store initial velocity direction angle (not position angle)
    initialVelocityAngle = std::atan2(vy, vx);

    // Convert velocity to polar
    const float cos_phi{std::cos(phi)};
    const float sin_phi{std::sin(phi)};
    v_r = vx * cos_phi + vy * sin_phi;
    v_phi = (-vx * sin_phi + vy * cos_phi) / r;

    // Calculate conserved quantities (geometric units: RS = 1)
    L = r * r * v_phi;
    const float f{1.0f - 1.0f / r};
    const float dt_dlambda{std::sqrt((v_r * v_r) / (f * f) + (r * r * v_phi * v_phi) / f)};
    E = f * dt_dlambda;

    trail.push_back(glm::vec2(x, y));
}

bool Ray::isCaptured() const {
    return r <= 1.01f;  // Captured if within 1% of event horizon (RS = 1)
}

bool Ray::hasEscaped(float maxDistance) const {
    return r > maxDistance;
}

void Ray::recordPosition() {
    const float x{r * std::cos(phi)};
    const float y{r * std::sin(phi)};
    trail.push_back(glm::vec2(x, y));
}

void Ray::updateDeflection() {
    // Convert current polar velocity back to Cartesian to get velocity direction
    const float vx_current{v_r * std::cos(phi) - r * v_phi * std::sin(phi)};
    const float vy_current{v_r * std::sin(phi) + r * v_phi * std::cos(phi)};

    // Calculate current velocity direction angle
    const float currentVelocityAngle{std::atan2(vy_current, vx_current)};

    // Deflection = change in velocity direction
    deflection = std::abs(currentVelocityAngle - initialVelocityAngle);

    // Handle angle wrapping (deflections > π map back down)
    if (deflection > static_cast<float>(M_PI)) {
        deflection = 2.0f * static_cast<float>(M_PI) - deflection;
    }
}

bool Ray::isActive(int currentFrame) const {
    return currentFrame >= startFrame;
}