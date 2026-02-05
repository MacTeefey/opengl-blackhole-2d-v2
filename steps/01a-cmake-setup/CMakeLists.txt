cmake_minimum_required(VERSION 3.21)

project(blackhole)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# TODO: Find OpenGL package (REQUIRED)
# TODO: Find glfw3 package (CONFIG REQUIRED)

# TODO: Find GLEW package (REQUIRED)
# TODO: Find glm package (CONFIG REQUIRED)

# TODO: Create a DEPS variable with all library targets:
#   glfw GLEW::GLEW glm::glm OpenGL::GL

add_executable(blackhole main.cpp)

# TODO: Link libraries to the blackhole executable using ${DEPS}
