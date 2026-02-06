#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <iostream>

#ifndef M_PI // Surprised that C++ doesnt define Pi, but need to so here it is
#define M_PI 3.14159265358979323846 // This level of specificity is fine for our sim
#endif 

const int WIDTH{800};
const int HEIGHT{600};

// TODO: Implement drawCircle function
// Parameters: float x, float y, float radius, int segments
// Use GL_TRIANGLE_FAN, add center vertex first, then loop through angles

// Due to the inability to draw true circles, need to create method to draw in depth triangles
// x & y - center positon of circle
// radius - size of circle
// segments - Number of triangles drawn to simulate a circle (more = smoother) (should be 50-100)
void drawCircle(float x, float y, float radius, int segments) {
    glBegin(GL_TRIANGLE_FAN); // Begining of triangle shapes
    glVertex2f(x, y); // Center of triangle
    // Generate triangles to draw circle
    for (int i{0}; i <= segments; ++i) {
        // Divides angles into equal segments so even number of triangles per circle
        const float angle{(static_cast<float>(i) * 2.0f * static_cast<float>(M_PI)) / static_cast<float>(segments)};
        const float px{x + radius * std::cos(angle)}; // px is the x coordinate on circle's edge
        const float py{y + radius * std::sin(angle)}; // py is the y coordinate on circle's edge
        glVertex2f(px, py); // Forms triangle from coordinates
    }

    glEnd();
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

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-10.0, 10.0, -10.0, 10.0, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        glColor3f(0.2f, 0.4f, 0.8f); // Blue color
        drawCircle(0.0f, 0.0f, 3.0f, 100); // circle at origin w/ radius 3.0

        glColor3f(1.0f, 0.5f, 0.0f); // Orange color
        drawCircle(5.0f, 0.0f, 1.5f, 100); // circle at (5,0) w/ radius 1.5

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}