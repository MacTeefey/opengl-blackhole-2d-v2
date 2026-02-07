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

    Ray(float x, float y, float vx, float vy) {
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


    void integrate(float dlambda, float maxDistance) {
        if (isCaptured() || hasEscaped(maxDistance)) {
            return;
        }
        verletStep(*this, dlambda);
        recordPosition();
        updateDeflection();
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

// This function computes radial and angular accelerations directly:
// Compute accelerations for geodesic motion
// Returns radial and angular accelerations from Schwarzschild geodesic equations
void computeAccelerations(float r, float v_r, float v_phi, float E_val,
                          float& a_r, float& a_phi) {
    float f{1.0f - 1.0f / r};
    float dt_dlambda{E_val / f};

    // Radial acceleration from Schwarzschild geodesic equation
    a_r = -(1.0f / (2.0f * r * r)) * f * (dt_dlambda * dt_dlambda)
          + (1.0f / (2.0f * r * r * f)) * (v_r * v_r)
          + (r - 1.0f) * (v_phi * v_phi);

    // Angular acceleration: dv_phi/dlambda = -2 * v_r * v_phi / r
    a_phi = -2.0f * v_r * v_phi / r;
}

// Symplectic Velocity Verlet (Stormer-Verlet) integration step
// Uses "kick-drift-kick" scheme that preserves phase space structure
void verletStep(Ray& ray, float dlambda) {
    float a_r, a_phi;

    // KICK: Half-step velocity update using initial accelerations
    computeAccelerations(ray.r, ray.v_r, ray.v_phi, ray.E, a_r, a_phi);
    float v_r_half{ray.v_r + 0.5f * dlambda * a_r};
    float v_phi_half{ray.v_phi + 0.5f * dlambda * a_phi};

    // DRIFT: Full-step position update using half-step velocities
    ray.r += dlambda * v_r_half;
    ray.phi += dlambda * v_phi_half;

    // KICK: Half-step velocity update using new accelerations
    computeAccelerations(ray.r, v_r_half, v_phi_half, ray.E, a_r, a_phi);
    ray.v_r = v_r_half + 0.5f * dlambda * a_r;
    ray.v_phi = v_phi_half + 0.5f * dlambda * a_phi;
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
void drawRays(const std::vector<Ray>& rays) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.5f);

    for (const auto& ray : rays) {
        if (ray.trail.size() < 2) {
            continue;
        }

        // Calculate color based on velocity deflection angle (0 to π radians)
        const float maxDeflection{static_cast<float>(M_PI)};
        const float t{std::min(1.0f, ray.deflection / maxDeflection)};

        // Cyan (minimal deflection) → Purple (moderate) → Red (extreme deflection)
        const float r{t};              // 0 → 1
        const float g{0.5f * (1.0f - t)};  // 0.5 → 0
        const float b{1.0f - t};       // 1 → 0

        // Draw trail with fading
        glBegin(GL_LINE_STRIP);
        for (size_t i{0}; i < ray.trail.size(); ++i) {
            const float alpha{0.2f + 0.8f * static_cast<float>(i) / static_cast<float>(ray.trail.size() - 1)};

            glColor4f(r, g, b, alpha);
            glVertex2f(ray.trail[i].x, ray.trail[i].y);
        }
        glEnd();
    }

    // Draw current ray positions as bright yellow points
    glPointSize(3.0f);
    glColor3f(1.0f, 1.0f, 0.0f);
    glBegin(GL_POINTS);
    for (const auto& ray : rays) {
        if (!ray.trail.empty() && !ray.isCaptured()) {
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

    // Ray launch parameters in geometric units
    // Position: simple values (units of RS)
    // Velocity: 1.0f (speed of light in geometric units)
    const int numRays{30};
    const float startX{-viewWidth};  // Left edge of viewport
    const float vx{1.0f};   // Speed of light = 1 in geometric units
    const float vy{0.0f};

    for (int i{0}; i < numRays; ++i) {
        const float startY{2.5f + i * 0.05f};

        rays.emplace_back(startX, startY, vx, vy);

        std::cout << "Ray " << i << ": y = " << startY << " RS\n";
    }

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

        // Integrate all rays (each ray checks if it should continue)
        // Smaller step size (0.05f) needed for geometric units
        for (auto& ray : rays) {
            ray.integrate(0.05f, maxDistance);
        }

        // Color-coded ray trails
        drawRays(rays);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
