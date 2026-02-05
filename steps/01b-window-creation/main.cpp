#include <GLFW/glfw3.h>
#include <iostream>

int main() {
    // Initialize GLFW library
    // Check if initialization failed and return -1 if so
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW";
        return -1;
    }

    // Create a window (800x600) titled "2D Black Hole Simulator"
    // Store it in a GLFWwindow pointer
    // Check if window creation failed, clean up GLFW, and return -1 if so
    GLFWwindow* window{ glfwCreateWindow(800, 600, "2D Black Hole Simulator", nullptr, nullptr) };
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    // Clean up GLFW resources before exiting
    glfwTerminate();
    return 0;

}
