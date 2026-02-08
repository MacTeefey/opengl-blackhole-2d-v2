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

enum class RayScenario {
    PARALLEL,      // Parallel rays from left side (distant starlight)
    POINT_SOURCE,  // Rays emanating from a single point (gravitational lensing)
    ORBITING       // Special orbiting ray (photon sphere)
};

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

    Ray(double x, double y, double vx, double vy, RayScenario scenario, int startFrame = 0)
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

    bool isActive(int currentFrame) const {
        return currentFrame >= startFrame;
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

// Adds isActive check before isCaptured/hasEscaped checks
// Integrates ray by a fixed distance using adaptive substeps
// This maintains constant visual speed while using adaptive accuracy
void integrateWithBudget(Ray& ray, float dlambda, float maxDistance, int currentFrame) {
    if (!ray.isActive(currentFrame) || ray.isCaptured() || ray.hasEscaped(maxDistance)) {
        return;
    }

    float remaining{dlambda};
    while (remaining > 0.0f) {
        float step{adaptiveStep(ray.r)};
        step = std::min(step, remaining);  // Don't exceed remaining distance

        verletStep(ray, step);
        remaining -= step;

        if (ray.isCaptured() || ray.hasEscaped(maxDistance)) {
            break;
        }
    }
    ray.recordPosition();
    ray.updateDeflection();
}

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

// Draw ray trails with fading
void drawRays(const std::vector<Ray>& rays, int currentFrame) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.5f);

    for (const auto& ray : rays) {
        if (!ray.isActive(currentFrame) || ray.trail.size() < 2) {
            continue;
        }

        // Calculate color based on velocity deflection angle (0 to π radians)
        const float maxDeflection{static_cast<float>(M_PI)};
        const float t{std::min(1.0f, ray.deflection / maxDeflection)};

        // Color based on ray scenario
        float r, g, b;

        if (ray.scenario == RayScenario::POINT_SOURCE) {
            // Green to Yellow gradient based on deflection
            const float maxDeflection{static_cast<float>(M_PI)};
            const float t{std::min(1.0f, static_cast<float>(ray.deflection) / maxDeflection)};
            r = 0.5f + 0.5f * t;  // 0.5 to 1.0 (increases with deflection)
            g = 1.0f;             // Constant green
            b = 0.0f;
        } else if (ray.scenario == RayScenario::ORBITING) {
            // Solid bright magenta (no gradient)
            r = 1.0f;
            g = 0.2f;
            b = 1.0f;
        } else {
            // PARALLEL: Cyan to Purple to Red gradient
            const float maxDeflection{static_cast<float>(M_PI)};
            const float t{std::min(1.0f, static_cast<float>(ray.deflection) / maxDeflection)};
            r = t;                    // 0 to 1
            g = 0.5f * (1.0f - t);    // 0.5 to 0
            b = 1.0f - t;             // 1 to 0
        }

        // Draw trail with fading
        glBegin(GL_LINE_STRIP);
        for (size_t i{0}; i < ray.trail.size(); ++i) {
            const float alpha{0.2f + 0.8f * static_cast<float>(i) / static_cast<float>(ray.trail.size() - 1)};

            glColor4f(r, g, b, alpha);
            glVertex2f(ray.trail[i].x, ray.trail[i].y);
        }
        glEnd();
    }

    // Draw current ray positions with scenario-specific colors
    glPointSize(3.0f);
    glBegin(GL_POINTS);
    for (const auto& ray : rays) {
        if (ray.isActive(currentFrame) && !ray.trail.empty() && !ray.isCaptured()) {
            // Color by scenario
            if (ray.scenario == RayScenario::POINT_SOURCE) {
                glColor3f(0.5f, 1.0f, 0.0f);  // Lime green
            } else if (ray.scenario == RayScenario::ORBITING) {
                glColor3f(1.0f, 0.2f, 1.0f);  // Magenta
            } else {
                glColor3f(1.0f, 1.0f, 0.0f);  // Yellow
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

    // TODO: Replace this simple ray generation with THREE scenarios:

    // SCENARIO 1: Special orbiting ray near photon sphere (1 ray)
    // Starts at frame 0 from left side of screen
    std::cout << "=== Scenario 1: Orbiting Ray ===\n";
    const float rs = 1.0f;
    const double orbitStartX{-0.9 * viewWidth}; // 90% left
    const double orbitStartY{2.577934 * rs}; // Tuned for multiple spirals
    const double orbitVx{1.0f}; // Horizontal toward black hole
    const double orbitVy{0.0};
    rays.emplace_back(orbitStartX, orbitStartY, orbitVx, orbitVy, RayScenario::ORBITING, 0);

    std::cout << "  Orbiting ray starts at x = " << orbitStartX / 1e9 << " Gm\n";
    std::cout << "  Impact parameter = " << orbitStartY / rs << " rs (creates multiple spirals!)\n";
    std::cout << "  Horizontal velocity = " << orbitVx / 1.0f << " c\n";
    std::cout << "  Starts at frame: 0\n";
    std::cout << "  Expected: Multiple dramatic spiral orbits before escape\n\n";

    // SCENARIO 2: Point source from top-left corner (25 rays)
    // Starts at frame 180
    std::cout << "=== Scenario 2: Point Source from Top-Left ===\n";

    const double sourceX{-0.95 * viewWidth};  // 95% left for dramatic lensing
    const double sourceY{0.85 * viewHeight};  // 85% up
    const int numPointRays{25};

    // Calculate direction from source to black hole
    const double toBlackHoleX{0.0 - sourceX};
    const double toBlackHoleY{0.0 - sourceY};
    const double baseAngle{std::atan2(toBlackHoleY, toBlackHoleX)};

    std::cout << "  Source position: (" << sourceX / 1e9 << ", " << sourceY / 1e9 << ") Gm\n";
    std::cout << "  Aiming toward black hole at angle: " << baseAngle * 180.0 / M_PI << " degrees\n";
    std::cout << "  Starts at frame: 180\n";

    for (int i{0}; i < numPointRays; ++i) {
        // Spread rays in 60-degree cone around direction to black hole
        const double coneSpread{M_PI / 3.0};  // 60 degrees total (π/3 radians)
        const double angleOffset{-coneSpread / 2.0 + coneSpread * static_cast<double>(i) / static_cast<double>(numPointRays - 1)};
        const double angle{baseAngle + angleOffset};

        // Calculate velocity components (all rays travel at speed c)
        const double vx{1.0f * std::cos(angle)};
        const double vy{1.0f * std::sin(angle)};

        rays.emplace_back(sourceX, sourceY, vx, vy, RayScenario::POINT_SOURCE, 180);

        // Log every 5th ray
        if (i % 5 == 0) {
            std::cout << "  Point ray " << i << ": angle = " << angle * 180.0 / M_PI
                      << " degrees (offset " << angleOffset * 180.0 / M_PI << ")\n";
        }
    }

    std::cout << "\n";

    // SCENARIO 3: Parallel rays from left (70 rays)
    // Starts at frame 360
    std::cout << "=== Scenario 3: Parallel Rays ===\n";

    const int numParallelRays{70};
    const double parallelStartX{-1e11};  // 100 Gm to the left
    const double parallelVx{1.0f};          // Horizontal, rightward
    const double parallelVy{0.0};        // No vertical component

    std::cout << "  Starts at frame: 360\n";
    
    for (int i{0}; i < numParallelRays; ++i) {
        // Calculate interpolation parameter [0, 1]
        const double t{static_cast<double>(i) / static_cast<double>(numParallelRays - 1)};

        // Map to full viewport height: -viewHeight to +viewHeight
        const double startY{-viewHeight + t * 2.0 * viewHeight};

        rays.emplace_back(parallelStartX, startY, parallelVx, parallelVy,   RayScenario::PARALLEL, 360);

        // Log every 10th ray
        if (i % 10 == 0) {
            const double impactParam{std::abs(startY)};
            std::cout << "  Parallel ray " << i << ": y = " << startY / 1e9
                      << " Gm, impact = " << impactParam / rs << " rs\n";
        }
    }

    std::cout << "\n";


    // TODO: Print summary showing total count and breakdown by scenario type

    const int totalRays{static_cast<int>(rays.size())};
    std::cout << "=== Total: " << totalRays << " rays ===\n";
    std::cout << "  1 orbiting (magenta) - starts frame 0\n";
    std::cout << "  " << numPointRays << " point source (green-yellow gradient) - starts frame 180\n";
    std::cout << "  " << numParallelRays << " parallel (cyan-red gradient) - starts frame 360\n\n";

    // Fixed distance per frame = constant visual speed
    // Adaptive substeps = accuracy near black hole
    const float distancePerFrame{0.08f};

    int frame{0};

    std::cout << "=== Simulation Starting ===\n";
    std::cout << "Watch for sequential activation:\n";
    std::cout << "  Frame 0: Magenta orbiting ray near photon sphere begins\n";
    std::cout << "  Frame 180: Green-yellow point source rays appear, aimed at black hole\n";
    std::cout << "  Frame 360: Cyan-red parallel rays fill the screen\n";
    std::cout << "  - Some rays captured by event horizon, others escape\n\n";


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

        // Mark the point source location (bright green)
        glPointSize(8.0f);
        glColor3f(0.5f, 1.0f, 0.0f);
        glBegin(GL_POINTS);
        glVertex2f(static_cast<float>(sourceX), static_cast<float>(sourceY));
        glEnd();


        // Integrate with adaptive stepping
        for (auto& ray : rays) {
            integrateWithBudget(ray, distancePerFrame, maxDistance, frame);
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
