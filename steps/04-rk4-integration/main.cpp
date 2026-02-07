#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <cmath>
#include <iostream>

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

        std::cout << "\nRay initialized:\n";
        std::cout << "  Position: r = " << (r / rs) << " rs\n";
        std::cout << "  Energy: E = " << E << "\n";
        std::cout << "  Angular momentum: L = " << L << "\n\n";
    }

    void displayConservedQuantities() const {
        std::cout << "Conserved quantities:\n";
        std::cout << "  E = " << E << " (energy)\n";
        std::cout << "  L = " << L << " (angular momentum)\n";

        // Show that these define the ray's behavior
        double impactParameter{L / E};
        std::cout << "  Impact parameter b = " << impactParameter / rs << " rs\n";
    }
};

// This computes the right-hand side of the geodesic differential equations
// Parameters:
//   ray - the current state
//   rhs - output array for [dr/dλ, dφ/dλ, dv_r/dλ, dv_phi/dλ]
void geodesicRHS(const Ray& ray, double rhs[4]) {
    double r{ray.r};
    double v_r{ray.v_r};
    double v_phi{ray.v_phi};
    double E{ray.E};

    double f_rhs{1.0 - rs / r}; // Schwartzschild metric coefficient
    double dt_dlambda{E / f_rhs}; // change in time over affine parameter

    // rhs[0] = dr/dλ (position derivative = velocity)
    rhs[0] = v_r;
    // rhs[1] = dφ/dλ (angular position derivative = angular velocity)
    rhs[1] = v_phi;

    // rhs[2] = dv_r/dλ (radial acceleration from geodesic equation)
    // Exact Schwarzschild null geodesic equation:
    // dv_r/dλ = -(rs/(2r²))f(dt/dλ)² + (rs/(2r²f))(v_r)² + (r - rs)(v_φ)²
    rhs[2] = -(rs / (2.0 * r * r)) * f_rhs * (dt_dlambda * dt_dlambda) // Gravitational attraction (negative)
            + (rs / (2.0 * r * r * f_rhs)) * (v_r * v_r) // Radial Velocity Correction (Positive)
            + (r - rs) * (v_phi * v_phi); // Centrifugal Effect (Positive)

    // rhs[3] = dv_phi/dλ (angular acceleration from geodesic equation)
    // dv_φ/dλ = -2v_r·v_φ/r
    rhs[3] = -2.0 * v_r * v_phi / r;

}

// This adds two state vectors: out = a + b * factor
// Parameters:
//   a[4] - first state vector
//   b[4] - second state vector
//   factor - scaling factor for b
//   out[4] - result
void addState(const double a[4], const double b[4], double factor, double out[4]) {
    for (int i{0}; i < 4; ++i) {
        out[i] = a[i] + b[i] * factor;
    }
}

// This performs one RK4 integration step
// Parameters:
//   ray - the ray to update (modified in place)
//   dlambda - step size along affine parameter
void rk4Step(Ray& ray, double dlambda) {
    // Current state
    double y0[4]{ray.r, ray.phi, ray.v_r, ray.v_phi};
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    // k1: Slope at current point
    geodesicRHS(ray, k1);

    // k2: Slope at midpoint using k1
    addState(y0, k1, dlambda / 2.0, temp);
    Ray r2{ray};
    r2.r = temp[0];
    r2.phi = temp[1];
    r2.v_r = temp[2];
    r2.v_phi = temp[3];
    geodesicRHS(r2, k2);

    // k3: Slope at midpoint using k2 (better estimate which is key to RK4)
    addState(y0, k2, dlambda / 2.0, temp);
    Ray r3{ray};
    r3.r = temp[0];
    r3.phi = temp[1];
    r3.v_r = temp[2];
    r3.v_phi = temp[3];
    geodesicRHS(r3, k3);

    // k4: Slope at endpoint using k3
    addState(y0, k3, dlambda, temp);
    Ray r4{ray};
    r4.r = temp[0];
    r4.phi = temp[1];
    r4.v_r = temp[2];
    r4.v_phi = temp[3];
    geodesicRHS(r4, k4);

    // Final weighted average: (k1 + 2*k2 + 2*k3 + k4) / 6
    ray.r     += (dlambda / 6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
    ray.phi   += (dlambda / 6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
    ray.v_r   += (dlambda / 6.0) * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]);
    ray.v_phi += (dlambda / 6.0) * (k1[3] + 2.0*k2[3] + 2.0*k3[3] + k4[3]);

}

// This function tests the RK4 integration and verifies conserved quantities
// Parameters:
//   ray - the ray to integrate (passed by value, so original is not modified)
//   num_steps - number of integration steps to perform
//   dlambda - step size for each integration step
void debugRK4Integration(Ray ray, int num_steps, double dlambda) {
    std::cout << "\n=== RK4 Integration Test ===\n";
    std::cout << "Step size (dlambda): " << dlambda << "\n";
    std::cout << "Number of steps: " << num_steps << "\n\n";

    // Store initial conserved quantities
    double initial_L{ray.L};
    double initial_E{ray.E};

    std::cout << "Initial state:\n";
    std::cout << "  r = " << (ray.r / rs) << " rs\n";
    std::cout << "  φ = " << ray.phi << " rad\n";
    std::cout << "  v_r = " << ray.v_r << "\n";
    std::cout << "  v_φ = " << ray.v_phi << "\n\n";

    // Run integration steps
    for (int i{0}; i < num_steps; ++i) {
        rk4Step(ray, dlambda);
    }

    std::cout << "After " << num_steps << " RK4 steps:\n";
    std::cout << "  r = " << (ray.r / rs) << " rs\n";
    std::cout << "  φ = " << ray.phi << " rad\n";
    std::cout << "  v_r = " << ray.v_r << "\n";
    std::cout << "  v_φ = " << ray.v_phi << "\n\n";

    // Recompute conserved quantities
    double current_L{ray.r * ray.r * ray.v_phi};
    double f{1.0 - rs / ray.r};
    double dt_dlambda{std::sqrt((ray.v_r * ray.v_r) / (f * f) + (ray.r * ray.r * ray.v_phi * ray.v_phi) / f)};
    double current_E{f * dt_dlambda};

    std::cout << "Conserved quantity check:\n";
    std::cout << "  Initial L = " << initial_L << "\n";
    std::cout << "  Current L = " << current_L << "\n";
    std::cout << "  Initial E = " << initial_E << "\n";
    std::cout << "  Current E = " << current_E << "\n";

    // Calculate drift percentages
    double L_drift{std::abs(current_L - initial_L) / std::abs(initial_L) * 100.0};
    double E_drift{std::abs(current_E - initial_E) / std::abs(initial_E) * 100.0};

    std::cout << "  L drift: " << L_drift << "%\n";
    std::cout << "  E drift: " << E_drift << "%\n";
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

    // Display black hole information
    std::cout << "\n=== Black Hole Properties ===\n";
    std::cout << "Mass: " << BH_MASS << " kg (Sagittarius A*)\n";
    std::cout << "Schwarzschild radius: " << rs << " m\n";
    std::cout << "  = " << rs / 1000.0 << " km\n";
    std::cout << "  = " << rs / 1e9 << " million km\n";
    std::cout << "Photon sphere: " << 1.5 * rs / 1e9 << " million km\n\n";

    // Create test ray - horizontal ray approaching from left
    double initialX{-1e11};
    double initialY{-3.27606302719999999e10};
    double velocityX{c};
    double velocityY{0.0};

    Ray testRay{initialX, initialY, velocityX, velocityY};
    testRay.displayConservedQuantities();

    // Run integration test
    const double dlambda{1.0};
    const int num_steps{10};
    debugRK4Integration(testRay, num_steps, dlambda);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        const float aspect{static_cast<float>(WIDTH) / static_cast<float>(HEIGHT)};

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-10.0 * aspect, 10.0 * aspect, -10.0, 10.0, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Event horizon
        glColor3f(0.0f, 0.0f, 0.0f);
        drawCircle(0.0f, 0.0f, 1.0f, 100);

        // Photon sphere outline
        glColor3f(0.0f, 0.8f, 0.8f);
        glLineWidth(2.0f);
        drawCircleOutline(0.0f, 0.0f, 1.5f, 100);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
