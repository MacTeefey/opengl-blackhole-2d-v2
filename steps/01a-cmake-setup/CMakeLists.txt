cmake_minimum_required(VERSION 3.21)

project(blackhole)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenGL package (REQUIRED)
find_package(OpenGL REQUIRED)
# Find glfw3 package (CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)

# Find GLEW package (REQUIRED)
find_package(GLEW REQUIRED)
# Find glm package (CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)

# Create a DEPS variable with all library targets:
#   glfw GLEW::GLEW glm::glm OpenGL::GL
set(DEPS glfw GLEW::GLEW glm::glm OpenGL::GL)

add_executable(blackhole main.cpp)

# Link libraries to the blackhole executable using ${DEPS}
target_link_libraries(blackhole PRIVATE ${DEPS})