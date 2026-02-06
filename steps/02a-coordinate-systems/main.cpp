#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>

const int WIDTH{800};
const int HEIGHT{600};


int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    // Uses constants to allow adjustability at top of file
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

    // Set up viewport using glViewport()
    int fbWidth, fbHeight;
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    glViewport(0, 0, fbWidth, fbHeight);

    // Set up projection matrix with glOrtho()
    glMatrixMode(GL_PROJECTION); // Switch to projection mode
    glLoadIdentity(); // Reset matrix
    glOrtho(-10.0, 10.0, -10.0, 10.0, -1.0, 1.0); // Set orthographic projection from -10 to 10 in X and Y
    

    glMatrixMode(GL_MODELVIEW); // Switch back to modelview matrix
    glLoadIdentity(); // Reset matrix

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw a test rectangle using glBegin(GL_QUADS)
        glBegin(GL_QUADS);
        // Use glColor3f() to set color
        glColor3f(1.0f, 0.5f, 0.0f);  // Orange
        // Use glVertex2f() to add 4 corners
        glVertex2f(-2.0f, -2.0f);      // Bottom-left
        glVertex2f( 2.0f, -2.0f);      // Bottom-right
        glVertex2f( 2.0f,  2.0f);      // Top-right
        glVertex2f(-2.0f,  2.0f);      // Top-left
        glEnd();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
