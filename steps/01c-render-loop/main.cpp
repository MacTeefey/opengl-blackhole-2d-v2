#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    GLFWwindow* window{ glfwCreateWindow(800, 600, "2D Black Hole Simulator", nullptr, nullptr) };
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window); // Make the OpenGL context current for this window

    glfwSwapInterval(1);  // Enable VSync to limit frame rate to monitor refresh rate

    // TODO: Initialize GLEW
    // Set glewExperimental to GL_TRUE first
    // Then call glewInit and check if it returns GLEW_OK
    // If it fails, destroy window, terminate GLFW, and return -1
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return -1;
    }

    // TODO: Create render loop
    // Create a loop that runs until the window should close
    // Inside the loop:
    //   - Clear the screen (color buffer bit)
    //   - Swap the front and back buffers
    //   - Process pending events
    glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
    }
    glfwSwapBuffers(window);
    glfwPollEvents(); // Checks & Processes events (keyboard presses, clicks, window resize, close button, etc.)

    // TODO: Destroy the window before terminating GLFW
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
