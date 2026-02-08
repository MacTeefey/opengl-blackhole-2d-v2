#pragma once

#include <glm/glm.hpp>
#include <vector>

// Ray scenario types for different visual effects
enum class RayScenario {
    PARALLEL,
    POINT_SOURCE,
    ORBITING
};

// Ray representation in Schwarzschild coordinates (geometric units)
struct Ray {
    // Polar coordinates
    float r;      // Radial distance from black hole
    float phi;    // Angular position
    float v_r;    // Radial velocity (dr/dλ)
    float v_phi;  // Angular velocity (dφ/dλ)

    // Conserved quantities (fundamental to general relativity)
    float E;      // Conserved energy
    float L;      // Conserved angular momentum

    // Visualization data
    std::vector<glm::vec2> trail;
    float initialVelocityAngle;
    float deflection;

    // Scenario tracking
    RayScenario scenario;
    int startFrame;

    // Constructor: Initialize ray from Cartesian position and velocity
    Ray(float x, float y, float vx, float vy,
        RayScenario scenario = RayScenario::PARALLEL,
        int startFrame = 0);

    // Check if ray has been captured by black hole
    bool isCaptured() const;

    // Check if ray has escaped to infinity
    bool hasEscaped(float maxDistance) const;

    // Record current position to trail
    void recordPosition();

    // Update deflection angle (for color coding)
    void updateDeflection();

    // Check if ray should be active at current frame
    bool isActive(int currentFrame) const;
};