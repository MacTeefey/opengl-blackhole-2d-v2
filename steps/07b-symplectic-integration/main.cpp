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
// TODO: Replace rk4Step forward declaration with verletStep
void rk4Step(Ray& ray, float dlambda);

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
        // TODO: Update to call verletStep instead of rk4Step
        rk4Step(*this, dlambda);
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

// TODO: Replace geodesicRHS with computeAccelerations function
// The new function should compute radial and angular accelerations directly:
// void computeAccelerations(float r, float v_r, float v_phi, float E_val,
//                           float& a_r, float& a_phi)
// Geodesic equations in geometric units (RS = 1)
// Compute geodesic right-hand side using exact Schwarzschild equations
void geodesicRHS(const Ray& ray, float rhs[4]) {
    float r{ray.r};
    float v_r{ray.v_r};
    float v_phi{ray.v_phi};
    float E{ray.E};

    // f = 1 - 1/r (simplified from f = 1 - rs/r since RS = 1)
    float f_rhs{1.0f - 1.0f / r};
    float dt_dlambda{E / f_rhs};

    // rhs[0] = dr/dλ (position derivative = velocity)
    rhs[0] = v_r;

    // rhs[1] = dφ/dλ (angular position derivative = angular velocity)
    rhs[1] = v_phi;

    // rhs[2] = dv_r/dλ (radial acceleration from geodesic equation)
    // Exact Schwarzschild null geodesic equation (with RS = 1):
    // dv_r/dλ = -(1/(2r²))f(dt/dλ)² + (1/(2r²f))(v_r)² + (r - 1)(v_φ)²
    rhs[2] = -(1.0f / (2.0f * r * r)) * f_rhs * (dt_dlambda * dt_dlambda)
             + (1.0f / (2.0f * r * r * f_rhs)) * (v_r * v_r)
             + (r - 1.0f) * (v_phi * v_phi);

    // rhs[3] = dv_phi/dλ (angular acceleration from geodesic equation)
    // dv_φ/dλ = -2v_r·v_φ/r
    rhs[3] = -2.0f * v_r * v_phi / r;
}

// TODO: Delete addState function (not needed for Verlet)
// Helper function: add two state vectors
void addState(const float a[4], const float b[4], float factor, float out[4]) {
    for (int i{0}; i < 4; ++i) {
        out[i] = a[i] + b[i] * factor;
    }
}

// TODO: Replace rk4Step with verletStep using Velocity Verlet algorithm
// Velocity Verlet uses a "kick-drift-kick" scheme:
// 1. Half-kick: update velocities by half step using current accelerations
// 2. Drift: update positions using the half-step velocities
// 3. Half-kick: update velocities by another half step using new accelerations
// RK4 integration step
void rk4Step(Ray& ray, float dlambda) {
    float y0[4]{ray.r, ray.phi, ray.v_r, ray.v_phi};
    float k1[4], k2[4], k3[4], k4[4], temp[4];

    // k1 = f(y0)
    geodesicRHS(ray, k1);

    // k2 = f(y0 + k1*dlambda/2)
    addState(y0, k1, dlambda / 2.0f, temp);
    Ray r2{ray};
    r2.r = temp[0];
    r2.phi = temp[1];
    r2.v_r = temp[2];
    r2.v_phi = temp[3];
    geodesicRHS(r2, k2);

    // k3 = f(y0 + k2*dlambda/2)
    addState(y0, k2, dlambda / 2.0f, temp);
    Ray r3{ray};
    r3.r = temp[0];
    r3.phi = temp[1];
    r3.v_r = temp[2];
    r3.v_phi = temp[3];
    geodesicRHS(r3, k3);

    // k4 = f(y0 + k3*dlambda)
    addState(y0, k3, dlambda, temp);
    Ray r4{ray};
    r4.r = temp[0];
    r4.phi = temp[1];
    r4.v_r = temp[2];
    r4.v_phi = temp[3];
    geodesicRHS(r4, k4);

    // Update: y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4) * dlambda/6
    ray.r     += (dlambda / 6.0f) * (k1[0] + 2.0f * k2[0] + 2.0f * k3[0] + k4[0]);
    ray.phi   += (dlambda / 6.0f) * (k1[1] + 2.0f * k2[1] + 2.0f * k3[1] + k4[1]);
    ray.v_r   += (dlambda / 6.0f) * (k1[2] + 2.0f * k2[2] + 2.0f * k3[2] + k4[2]);
    ray.v_phi += (dlambda / 6.0f) * (k1[3] + 2.0f * k2[3] + 2.0f * k3[3] + k4[3]);
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
