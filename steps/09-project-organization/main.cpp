#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Project headers - these must be implemented for the project to compile!
// TODO: Create constants.h with BlackHole, Visual, and Simulation namespaces
#include "constants.h"
// TODO: Create ray.h with RayScenario enum and Ray struct
#include "ray.h"
// TODO: Create physics.h with Physics namespace functions
#include "physics.h"
// TODO: Create rendering.h with Rendering namespace and RenderEngine
#include "rendering.h"

// Your task: Refactor this working single-file program into organized modules
// The project won't compile until you implement all the header files above!

const int WIDTH{800};
const int HEIGHT{600};

// Geometric units: c = G = 1, Schwarzschild radius = 1
const float RS{1.0f};

const int NUM_STARS{200};

// Forward declarations (needed for Ray::integrate method)
struct Ray;
void verletStep(Ray& ray, float dlambda);

void drawCircle(float x, float y, float radius, int segments) {
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y);

    for (int i{0}; i <= segments; ++i) {
        const float angle{(static_cast<float>(i) * 2.0f * static_cast<float>(M_PI)) / static_cast<float>(segments)};
        const float px{x + radius * std::cos(angle)};
        const float py{y + radius * std::sin(angle)};
        glVertex2f(px, py);
    }

    glEnd();
}

void drawCircleOutline(float x, float y, float radius, int segments) {
    glBegin(GL_LINE_LOOP);

    for (int i{0}; i < segments; ++i) {
        const float angle{(static_cast<float>(i) * 2.0f * static_cast<float>(M_PI)) / static_cast<float>(segments)};
        const float ox{x + radius * std::cos(angle)};
        const float oy{y + radius * std::sin(angle)};
        glVertex2f(ox, oy);
    }

    glEnd();
}

void drawDashedCircle(float x, float y, float radius, int segments) {
    for (int i{0}; i < segments; i += 2) {
        const float theta1{(static_cast<float>(i) * 2.0f * static_cast<float>(M_PI)) / static_cast<float>(segments)};
        const float theta2{(static_cast<float>(i + 1) * 2.0f * static_cast<float>(M_PI)) / static_cast<float>(segments)};

        glBegin(GL_LINES);
        glVertex2f(x + radius * std::cos(theta1), y + radius * std::sin(theta1));
        glVertex2f(x + radius * std::cos(theta2), y + radius * std::sin(theta2));
        glEnd();
    }
}

// TODO: Move to ray.h
enum class RayScenario {
    PARALLEL,
    POINT_SOURCE,
    ORBITING
};

// TODO: Move to ray.h (declaration) and ray.cpp (implementation)
struct Ray {
    // Polar coordinates
    float r;      // Radial distance from black hole
    float phi;    // Angular position
    float v_r;    // Radial velocity (dr/dλ)
    float v_phi;  // Angular velocity (dφ/dλ)

    // Conserved quantities (fundamental to general relativity)
    float E;      // Conserved energy
    float L;      // Conserved angular momentum

    std::vector<glm::vec2> trail;

    float initialVelocityAngle;
    float deflection{0.0f};

    RayScenario scenario;
    int startFrame;

    Ray(float x, float y, float vx, float vy, RayScenario scenario, int startFrame = 0)
        : scenario{scenario}, startFrame{startFrame} {
        // Step 1: Convert Cartesian position to polar coordinates
        r = std::sqrt(x * x + y * y);
        phi = std::atan2(y, x);

        // Step 2: Convert Cartesian velocity to polar coordinates
        // Using transformation:
        // v_r = vx*cos(φ) + vy*sin(φ)
        // v_φ = (-vx*sin(φ) + vy*cos(φ)) / r
        float cos_phi = std::cos(phi);
        float sin_phi = std::sin(phi);
        v_r = vx * cos_phi + vy * sin_phi;
        v_phi = (-vx * sin_phi + vy * cos_phi) / r;

        // Step 3: Calculate conserved quantities
        // These remain constant along the geodesic (ray path)

        // Angular momentum L = r² * v_φ
        L = r * r * v_phi;

        // Metric coefficient f = 1 - 1/r (simplified from f = 1 - rs/r since RS = 1)
        float f{1.0f - 1.0f / r};

        // For null geodesics (light rays), the normalization condition gives:
        // dt/dλ = √[(v_r)²/f² + (r²·v_φ)²/f]
        float dt_dlambda{std::sqrt((v_r * v_r) / (f * f) + (r * r * v_phi * v_phi) / f)};

        // Energy E = f * dt/dλ
        E = f * dt_dlambda;

        trail.push_back(glm::vec2(x, y));
    }

    bool isCaptured() const {
        return r <= RS * 1.01f;
    }

    bool hasEscaped(float maxDistance) const {
        return r > maxDistance;
    }

    void recordPosition() {
        const float x{r * std::cos(phi)};
        const float y{r * std::sin(phi)};
        trail.push_back(glm::vec2(x, y));
    }

    bool isActive(int currentFrame) const {
        return currentFrame >= startFrame;
    }

    void updateDeflection() {
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
};

// TODO: Move to physics.h (declaration) and physics.cpp (implementation)
// Compute radial and angular accelerations from Schwarzschild geodesic equations
// Geodesic equations in geometric units (RS = 1)
void computeAccelerations(float r, float v_r, float v_phi, float E_val,
                          float& a_r, float& a_phi) {
    // f = 1 - 1/r (simplified from f = 1 - rs/r since RS = 1)
    float f{1.0f - 1.0f / r};
    float dt_dlambda{E_val / f};

    // a_r = dv_r/dλ (radial acceleration from geodesic equation)
    // Exact Schwarzschild null geodesic equation (with RS = 1):
    // dv_r/dλ = -(1/(2r²))f(dt/dλ)² + (1/(2r²f))(v_r)² + (r - 1)(v_φ)²
    a_r = -(1.0f / (2.0f * r * r)) * f * (dt_dlambda * dt_dlambda)
          + (1.0f / (2.0f * r * r * f)) * (v_r * v_r)
          + (r - 1.0f) * (v_phi * v_phi);

    // a_phi = dv_phi/dλ (angular acceleration from geodesic equation)
    // dv_φ/dλ = -2v_r·v_φ/r
    a_phi = -2.0f * v_r * v_phi / r;
}

// TODO: Move to physics.h (declaration) and physics.cpp (implementation)
// Velocity Verlet integration step (symplectic integrator)
// Uses "kick-drift-kick" scheme for better energy conservation
void verletStep(Ray& ray, float dlambda) {
    float half_dt{dlambda / 2.0f};

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

// TODO: Move to physics.h (declaration) and physics.cpp (implementation)
// Adaptive timestep: smaller steps near black hole for accuracy
// Returns smaller step sizes near the black hole (r close to 1) for accuracy
// Returns larger step sizes far away (r > 10) for performance
float adaptiveStep(float r) {
    const float minStep{0.01f};
    const float maxStep{0.1f};
    // Scale linearly based on (r - 1.0f) / 10.0f, clamped to [0, 1]
    float factor{std::min((r - 1.0f) / 10.0f, 1.0f)};
    factor = std::max(factor, 0.0f);  // Clamp to non-negative
    return minStep + factor * (maxStep - minStep);
}

// TODO: Move to physics.h (declaration) and physics.cpp (implementation)
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

// TODO: Move to rendering.h (declaration) and rendering.cpp (implementation)
std::vector<glm::vec2> generateStars(int count, float viewWidth, float viewHeight) {
    std::vector<glm::vec2> stars;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distX(-viewWidth, viewWidth);
    std::uniform_real_distribution<float> distY(-viewHeight, viewHeight);

    for (int i{0}; i < count; ++i) {
        stars.push_back(glm::vec2(distX(gen), distY(gen)));
    }
    return stars;
}

// TODO: Move to rendering.h (declaration) and rendering.cpp (implementation)
void drawStars(const std::vector<glm::vec2>& stars) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> brightness(0.5f, 1.0f);

    glPointSize(2.0f);
    glBegin(GL_POINTS);
    for (const auto& star : stars) {
        float b = brightness(gen);
        glColor3f(b, b, b);  // Grayscale brightness
        glVertex2f(star.x, star.y);
    }
    glEnd();
}

// TODO: Move to rendering.h (declaration) and rendering.cpp (implementation)
// Draw ray trails with fading
void drawRays(const std::vector<Ray>& rays, int currentFrame) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.5f);

    for (const auto& ray : rays) {
        if (!ray.isActive(currentFrame)) {
            continue;
        }

        if (ray.trail.size() < 2) {
            continue;
        }

        // Color based on scenario type
        float r_color, g_color, b_color;
        const float maxDeflection{static_cast<float>(M_PI)};
        const float t{std::min(1.0f, ray.deflection / maxDeflection)};

        switch (ray.scenario) {
            case RayScenario::POINT_SOURCE:
                // Green to yellow gradient based on deflection
                r_color = t;           // 0 → 1
                g_color = 1.0f;        // Always green
                b_color = 0.0f;        // No blue
                break;
            case RayScenario::ORBITING:
                // Bright magenta (constant color)
                r_color = 1.0f;
                g_color = 0.0f;
                b_color = 1.0f;
                break;
            case RayScenario::PARALLEL:
            default:
                // Cyan (minimal deflection) → Purple (moderate) → Red (extreme deflection)
                r_color = t;                    // 0 → 1
                g_color = 0.5f * (1.0f - t);    // 0.5 → 0
                b_color = 1.0f - t;             // 1 → 0
                break;
        }

        // Draw trail with fading
        glBegin(GL_LINE_STRIP);
        for (size_t i{0}; i < ray.trail.size(); ++i) {
            const float alpha{0.2f + 0.8f * static_cast<float>(i) / static_cast<float>(ray.trail.size() - 1)};

            glColor4f(r_color, g_color, b_color, alpha);
            glVertex2f(ray.trail[i].x, ray.trail[i].y);
        }
        glEnd();
    }

    // Draw current ray positions with color by scenario
    glPointSize(3.0f);
    glBegin(GL_POINTS);
    for (const auto& ray : rays) {
        if (!ray.isActive(currentFrame)) {
            continue;
        }

        if (!ray.trail.empty() && !ray.isCaptured()) {
            switch (ray.scenario) {
                case RayScenario::POINT_SOURCE:
                    glColor3f(1.0f, 1.0f, 0.0f);  // Yellow
                    break;
                case RayScenario::ORBITING:
                    glColor3f(1.0f, 0.0f, 1.0f);  // Magenta
                    break;
                case RayScenario::PARALLEL:
                default:
                    glColor3f(1.0f, 1.0f, 0.0f);  // Yellow
                    break;
            }
            glVertex2f(ray.trail.back().x, ray.trail.back().y);
        }
    }
    glEnd();

    glDisable(GL_BLEND);
}

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    GLFWwindow* window{glfwCreateWindow(WIDTH, HEIGHT, "2D Black Hole Simulator", nullptr, nullptr)};
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // Enable VSync to limit frame rate to monitor refresh rate

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return -1;
    }

    // Dark space background with better contrast for black hole visibility
    glClearColor(0.08f, 0.08f, 0.12f, 1.0f);

    int fbWidth, fbHeight;
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    glViewport(0, 0, fbWidth, fbHeight);

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    std::cout << "Schwarzschild radius: " << RS << " (geometric units)\n";

    // Viewport in geometric units (units of RS)
    const float viewWidth{8.0f};    // 8 Schwarzschild radii
    const float viewHeight{6.0f};   // 6 Schwarzschild radii
    const float maxDistance{30.0f}; // Escape when beyond 30 RS

    std::vector<glm::vec2> stars{generateStars(NUM_STARS, viewWidth, viewHeight)};

    std::vector<Ray> rays;

    // SCENARIO 1 - Special orbiting ray near photon sphere (1 ray, starts at frame 0)
    // Position just below the critical impact parameter (2.598 RS)
    const float orbitStartX{-0.9f * viewWidth};
    const float orbitStartY{2.577934f};  // Just below critical 2.598 RS
    rays.emplace_back(orbitStartX, orbitStartY, 1.0f, 0.0f, RayScenario::ORBITING, 0);
    std::cout << "Orbiting ray: y = " << orbitStartY << " RS (starts frame 0)\n";

    // SCENARIO 2 - Point source from top-left corner (25 rays, starts at frame 700)
    const float sourceX{-0.95f * viewWidth};
    const float sourceY{0.85f * viewHeight};
    const int numPointRays{25};
    const float spreadAngle{60.0f * static_cast<float>(M_PI) / 180.0f};  // 60 degrees in radians

    // Calculate angle from source to black hole center
    const float angleToCenter{std::atan2(-sourceY, -sourceX)};

    for (int i{0}; i < numPointRays; ++i) {
        // Spread rays evenly across the cone
        const float angleOffset{-spreadAngle / 2.0f + (static_cast<float>(i) / static_cast<float>(numPointRays - 1)) * spreadAngle};
        const float rayAngle{angleToCenter + angleOffset};

        // Velocity direction at speed 1.0 (speed of light)
        const float vx{std::cos(rayAngle)};
        const float vy{std::sin(rayAngle)};

        rays.emplace_back(sourceX, sourceY, vx, vy, RayScenario::POINT_SOURCE, 700);
    }
    std::cout << "Point source: " << numPointRays << " rays from (" << sourceX << ", " << sourceY << ") (starts frame 700)\n";

    // SCENARIO 3 - Parallel rays from left (70 rays, starts at frame 1200)
    const int numParallelRays{70};
    const float parallelStartX{-viewWidth};

    for (int i{0}; i < numParallelRays; ++i) {
        // Evenly distribute from -viewHeight to +viewHeight
        const float startY{-viewHeight + (static_cast<float>(i) / static_cast<float>(numParallelRays - 1)) * 2.0f * viewHeight};

        rays.emplace_back(parallelStartX, startY, 1.0f, 0.0f, RayScenario::PARALLEL, 1200);
    }
    std::cout << "Parallel rays: " << numParallelRays << " rays (starts frame 1200)\n";

    // Summary
    std::cout << "\nTotal rays: " << rays.size() << "\n";
    std::cout << "  - Orbiting: 1\n";
    std::cout << "  - Point source: " << numPointRays << "\n";
    std::cout << "  - Parallel: " << numParallelRays << "\n";

    // Fixed "budget" per frame = constant visual speed
    // Adaptive substeps within budget = accuracy near black hole
    const float budgetPerFrame{0.08f};

    // Frame counter for staggered ray activation
    int frame{0};

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        // Setup viewport in geometric units
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-viewWidth, viewWidth, -viewHeight, viewHeight, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Background stars (drawn first)
        drawStars(stars);

        // Event horizon at r = RS = 1
        glColor3f(0.0f, 0.0f, 0.0f);
        drawCircle(0.0f, 0.0f, RS, 100);

        // Photon sphere outline at r = 1.5 RS (dashed for visual distinction)
        glColor3f(0.0f, 0.8f, 0.8f);
        glLineWidth(2.0f);
        drawDashedCircle(0.0f, 0.0f, 1.5f * RS, 100);

        // Draw point source marker (bright green dot)
        glPointSize(8.0f);
        glColor3f(0.0f, 1.0f, 0.0f);
        glBegin(GL_POINTS);
        glVertex2f(sourceX, sourceY);
        glEnd();

        // Integrate with budget-based adaptive stepping
        for (auto& ray : rays) {
            integrateWithBudget(ray, budgetPerFrame, maxDistance, frame);
        }

        // Color-coded ray trails
        drawRays(rays, frame);

        glfwSwapBuffers(window);
        glfwPollEvents();

        ++frame;
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
