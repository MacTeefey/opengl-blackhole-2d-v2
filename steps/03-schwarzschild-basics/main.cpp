#include <GL/glew.h>
#include <GLFW/glfw3.h>
// TODO: Add #include <glm/glm.hpp> for GLM math library

#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const int WIDTH{800};
const int HEIGHT{600};

// TODO: Define physical constants here
// - Speed of light (c)
// - Gravitational constant (G)
// - Black hole mass (BH_MASS) for Sagittarius A*
// - Schwarzschild radius (rs = 2GM/c²)

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

// TODO: Create Ray struct
// The Ray struct should store:
// - Position in polar coordinates (r, phi)
// - Velocity in polar coordinates (v_r, v_phi)
// - Conserved quantities (E, L)
// - Constructor that takes Cartesian position (x, y) and velocity (vx, vy)

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

    // TODO: Print black hole properties
    // - Mass
    // - Schwarzschild radius in meters and kilometers

    // TODO: Create a test Ray
    // - Initialize with some test position and velocity
    // - Ray should print its conserved quantities when created

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
