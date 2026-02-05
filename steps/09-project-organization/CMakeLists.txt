cmake_minimum_required(VERSION 3.21)

project(blackhole)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm CONFIG REQUIRED)

set(DEPS glfw GLEW::GLEW glm::glm OpenGL::GL)

# TODO: Add ray.cpp, physics.cpp, rendering.cpp to SOURCES
set(SOURCES
    main.cpp
)

# TODO: Add constants.h, ray.h, physics.h, rendering.h to HEADERS
set(HEADERS
)

# TODO: Replace with: add_executable(blackhole ${SOURCES} ${HEADERS})
add_executable(blackhole main.cpp)

target_link_libraries(blackhole PRIVATE ${DEPS})
