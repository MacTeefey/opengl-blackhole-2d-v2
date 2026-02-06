#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const int WIDTH{800};
const int HEIGHT{600};

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


// No center vertex needed
// Loop condition: i < segments (not i <= segments)
void drawCircleOutline(float x, float y, float radius, int segments) {
    glBegin(GL_LINE_LOOP); // Uses GL_LINE_LOOP to automatically close loop on itself

    for (int i{0}; i < segments; ++i) {
        const float angle{(static_cast<float>(i) * 2.0f * static_cast<float>(M_PI)) / static_cast<float>(segments)}; // Also splits the circle into small angles based on number of segments
        const float ox{x + radius * std::cos(angle)}; // calculates the parameter x
        const float oy{y + radius * std::sin(angle)}; // and y using math functions
        glVertex2f(ox, oy); // Draws a line instead of a triangle
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
        
        glColor3f(0.0f, 0.0f, 0.0f); // black color to simulate even horizon
        drawCircle(0.0f, 0.0f, 1.0f, 100); // circle around origin w/ radius of 1

        glColor3f(0.0f, 0.8f, 0.8f); // cyan color to simulate photon sphere
        glLineWidth(2.0f);  // 2 pixels wide to see outline
        drawCircleOutline(0.0f, 0.0f, 1.5f, 100); // circle around origin w/ radius of size 1.5

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
