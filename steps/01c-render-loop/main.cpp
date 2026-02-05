// TODO: Add GLEW include here
#include <GLFW/glfw3.h>
#include <iostream>

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    GLFWwindow* window{glfwCreateWindow(800, 600, "2D Black Hole Simulator", nullptr, nullptr)};
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }

    // TODO: Make the OpenGL context current for this window
    // TODO: Enable VSync to limit frame rate to monitor refresh rate
    // Call glfwSwapInterval(1)

    // TODO: Initialize GLEW
    // Set glewExperimental to GL_TRUE first
    // Then call glewInit and check if it returns GLEW_OK
    // If it fails, destroy window, terminate GLFW, and return -1

    // TODO: Create render loop
    // Create a loop that runs until the window should close
    // Inside the loop:
    //   - Clear the screen (color buffer bit)
    //   - Swap the front and back buffers
    //   - Process pending events

    // TODO: Destroy the window before terminating GLFW

    glfwTerminate();
    return 0;
}
