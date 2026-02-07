#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <cmath>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const int WIDTH{800};
const int HEIGHT{600};

// Physical constants (SI units)
const double c{299792458.0};      // Speed of light (m/s)
const double G{6.67430e-11};      // Gravitational constant (m³/kg·s²)

// Black hole parameters - Sagittarius A* at the center of our galaxy
const double BH_MASS{8.54e36};    // Mass in kg
const double rs{2.0 * G * BH_MASS / (c * c)};  // Schwarzschild radius

// Forward declarations (needed for Ray::integrate method)
struct Ray;
void rk4Step(Ray& ray, double dlambda);


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

struct Ray {
    // Polar coordinates
    double r;      // Radial distance from black hole
    double phi;    // Angular position
    double v_r;    // Radial velocity (dr/dλ)
    double v_phi;  // Angular velocity (dφ/dλ)

    // Conserved quantities (fundamental to general relativity)
    double E;      // Conserved energy
    double L;      // Conserved angular momentum

    std::vector<glm::vec2> trail;

    Ray(double x, double y, double vx, double vy) {
        // Step 1: Convert Cartesian position to polar coordinates
        r = std::sqrt(x * x + y * y);
        phi = std::atan2(y, x);

        // Step 2: Convert Cartesian velocity to polar coordinates
        // Using transformation:
        // v_r = vx*cos(φ) + vy*sin(φ)
        // v_φ = (-vx*sin(φ) + vy*cos(φ)) / r
        double cos_phi = std::cos(phi);
        double sin_phi = std::sin(phi);
        v_r = vx * cos_phi + vy * sin_phi;
        v_phi = (-vx * sin_phi + vy * cos_phi) / r;

        // Step 3: Calculate conserved quantities
        // These remain constant along the geodesic (ray path)

        // Angular momentum L = r² * v_φ
        L = r * r * v_phi;

        // Metric coefficient (f = 1 - rs/r)
        double f{1.0 - rs / r};

        // For null geodesics (light rays), the normalization condition gives:
        // dt/dλ = √[(v_r)²/f² + (r²·v_φ)²/f]
        double dt_dlambda{std::sqrt((v_r * v_r) / (f * f) + (r * r * v_phi * v_phi) / f)};

        // Energy E = f * dt/dλ
        E = f * dt_dlambda;

        trail.push_back(glm::vec2(static_cast<float>(x), static_cast<float>(y))); // Make sure it starts at launch position
    }

    // Helper function to check if a ray has crossed the event horizon
    bool isCaptured() const {
        return r <= rs * 1.01;
    }

    // Helper function to check if a ray has escaped the event horizon
    bool hasEscaped(double maxDistance) const {
        return r > maxDistance;
    }


    // Helper function to convert polar to Cartesian and store in trail
    void recordPosition() {
        const float x{static_cast<float>(r * std::cos(phi))};
        const float y{static_cast<float>(r * std::sin(phi))};
        trail.push_back(glm::vec2(x, y));
    }


    // Helper function that checks capture/escape, calls rk4Step, and records position
    void integrate(double dlambda, double maxDistance) {
        if (isCaptured() || hasEscaped(maxDistance)) {
            return;
        }
        rk4Step(*this, dlambda);
        recordPosition();
    }

};

// Compute geodesic right-hand side using exact Schwarzschild equations
void geodesicRHS(const Ray& ray, double rhs[4]) {
    double r{ray.r};
    double v_r{ray.v_r};
    double v_phi{ray.v_phi};
    double E{ray.E};

    double f_rhs{1.0 - rs / r};
    double dt_dlambda{E / f_rhs};

    // rhs[0] = dr/dλ (position derivative = velocity)
    rhs[0] = v_r;

    // rhs[1] = dφ/dλ (angular position derivative = angular velocity)
    rhs[1] = v_phi;

    // rhs[2] = dv_r/dλ (radial acceleration from geodesic equation)
    // Exact Schwarzschild null geodesic equation:
    // dv_r/dλ = -(rs/(2r²))f(dt/dλ)² + (rs/(2r²f))(v_r)² + (r - rs)(v_φ)²
    rhs[2] = -(rs / (2.0 * r * r)) * f_rhs * (dt_dlambda * dt_dlambda)
             + (rs / (2.0 * r * r * f_rhs)) * (v_r * v_r)
             + (r - rs) * (v_phi * v_phi);

    // rhs[3] = dv_phi/dλ (angular acceleration from geodesic equation)
    // dv_φ/dλ = -2v_r·v_φ/r
    rhs[3] = -2.0 * v_r * v_phi / r;
}

// Helper function: add two state vectors
void addState(const double a[4], const double b[4], double factor, double out[4]) {
    for (int i{0}; i < 4; ++i) {
        out[i] = a[i] + b[i] * factor;
    }
}

// RK4 integration step
void rk4Step(Ray& ray, double dlambda) {
    double y0[4]{ray.r, ray.phi, ray.v_r, ray.v_phi};
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    // k1 = f(y0)
    geodesicRHS(ray, k1);

    // k2 = f(y0 + k1*dlambda/2)
    addState(y0, k1, dlambda / 2.0, temp);
    Ray r2{ray};
    r2.r = temp[0];
    r2.phi = temp[1];
    r2.v_r = temp[2];
    r2.v_phi = temp[3];
    geodesicRHS(r2, k2);

    // k3 = f(y0 + k2*dlambda/2)
    addState(y0, k2, dlambda / 2.0, temp);
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
    ray.r     += (dlambda / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    ray.phi   += (dlambda / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    ray.v_r   += (dlambda / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
    ray.v_phi += (dlambda / 6.0) * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]);
}

// Helper function to render all ray trails with fading
void drawRays(const std::vector<Ray>& rays) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.5f);

    // Draws ray trails with fading effect
    for (const auto& ray : rays) {
        if (ray.trail.size() < 2) { // Skips rays that cant draw lines
            continue;
        }
        glBegin(GL_LINE_STRIP); // Connect all points in sequence
        for (size_t i{0}; i < ray.trail.size(); ++i) {
            const float alpha{0.2f + 0.8f * static_cast<float>(i) / static_cast<float>(ray.trail.size() - 1)}; // gradient from 20% opacity to 100% opacity
            glColor4f(1.0f, 1.0f, 1.0f, alpha);
            glVertex2f(ray.trail[i].x, ray.trail[i].y);
        }
        glEnd();
    }

    // Draws current ray positions as bright points
    glPointSize(3.0f);
    glColor3f(1.0f, 1.0f, 0.0f); // Bright yellow color
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

    glClearColor(0.1f, 0.1f, 0.15f, 1.0f);

    int fbWidth, fbHeight;
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    glViewport(0, 0, fbWidth, fbHeight);

    // Create test ray - horizontal ray approaching from left
    double initialX{-1e11};
    double initialY{-3.27606302719999999e10};
    double velocityX{c};
    double velocityY{0.0};

    Ray testRay{initialX, initialY, velocityX, velocityY};

    // Viewport dimensions (physical coordinates)
    const double viewWidth{1e11};   // 100 billion meters
    const double viewHeight{7.5e10}; // 75 billion meters

    // maxDistance derived from viewport width (2.0 * viewWidth)
    const double maxDistance{2.0 * viewWidth};

    std::vector<Ray> rays; // A vector to store multiple Ray objects
    
    // Ray launch parameters:
    const int numRays{10};
    const double startX{-viewWidth};
    const double vx{c};
    const double vy{0.0};


    // Launch rays in a loop:
    for (int i{0}; i < numRays; ++i) {
        const double startY{(2.5 + i * 0.1) * rs};
        rays.emplace_back(startX, startY, vx, vy);
    }

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        const float aspect{static_cast<float>(WIDTH) / static_cast<float>(HEIGHT)};

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-viewWidth, viewWidth, -viewHeight, viewHeight, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Event horizon
        glColor3f(0.0f, 0.0f, 0.0f);
        const float eventRadius{static_cast<float>(rs)};
        drawCircle(0.0f, 0.0f, eventRadius, 100);
        
        // Photon sphere outline
        glColor3f(0.0f, 0.8f, 0.8f);
        glLineWidth(2.0f);
        const float photonRadius{static_cast<float>(1.5 * rs)};
        drawCircleOutline(0.0f, 0.0f, photonRadius, 100);
        
        // Integrates all rays by calling ray.integrate(1.0, maxDistance) for each ray
        for (auto& ray : rays) {
            ray.integrate(1.0, maxDistance);
        }

        drawRays(rays);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
