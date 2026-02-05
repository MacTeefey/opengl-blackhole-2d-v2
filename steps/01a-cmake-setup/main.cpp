#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

int main() {
    // This verifies all libraries are properly linked
    // We test linking but don't actually initialize (Docker has no display)

    // These symbols must resolve during linking:
    // - glfwInit from GLFW
    // - glGetString from OpenGL
    // - glm::vec3 from GLM (header-only)

    glm::vec3 testVector{1.0f, 2.0f, 3.0f};
    (void)testVector; // Suppress unused warning

    // Reference the symbols so linker checks them
    void* glfwPtr = (void*)&glfwInit;
    void* glPtr = (void*)&glGetString;
    (void)glfwPtr;
    (void)glPtr;

    return 0;
}
