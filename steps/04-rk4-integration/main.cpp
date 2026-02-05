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

// TODO: Implement geodesicRHS() function
// This computes the right-hand side of the geodesic differential equations
// Parameters:
//   ray - the current state
//   rhs - output array for [dr/dλ, dφ/dλ, dv_r/dλ, dv_phi/dλ]

// TODO: Implement addState() helper function
// This adds two state vectors: out = a + b * factor
// Parameters:
//   a[4] - first state vector
//   b[4] - second state vector
//   factor - scaling factor for b
//   out[4] - result

// TODO: Implement rk4Step() function
// This performs one RK4 integration step
// Parameters:
//   ray - the ray to update (modified in place)
//   dlambda - step size along affine parameter

// TODO: Implement debugRK4Integration() function
// This function tests the RK4 integration and verifies conserved quantities
// Parameters:
//   ray - the ray to integrate (passed by value, so original is not modified)
//   num_steps - number of integration steps to perform
//   dlambda - step size for each integration step

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

    // TODO: Verify the rays conserved quantities remain constant (ie low % drift) after integration
    // by using the debugging code provided.

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
