#include <GL/glew.h>
#include <GLFW/glfw3.h>

// TODO: Add #include <cmath> for std::cos and std::sin
#include <iostream>

// TODO: Define M_PI (use #ifndef guard since some compilers already define it)

const int WIDTH{800};
const int HEIGHT{600};

// TODO: Implement drawCircle function
// Parameters: float x, float y, float radius, int segments
// Use GL_TRIANGLE_FAN, add center vertex first, then loop through angles

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

        // TODO: Delete this code before drawing the circles.
        glBegin(GL_QUADS);
            glColor3f(1.0f, 0.5f, 0.0f);
            glVertex2f(-2.0f, -2.0f);
            glVertex2f( 2.0f, -2.0f);
            glVertex2f( 2.0f,  2.0f);
            glVertex2f(-2.0f,  2.0f);
        glEnd();

        // TODO: Draw a large blue circle at origin (0, 0) with radius 3.0
        // Use glColor3f(0.2f, 0.4f, 0.8f) for blue color

        // TODO: Draw a smaller orange circle at (5, 0) with radius 1.5
        // Use glColor3f(1.0f, 0.5f, 0.0f) for orange color

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
