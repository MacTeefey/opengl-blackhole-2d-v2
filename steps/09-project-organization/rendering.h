#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>
#include "ray.h"

#error "TODO: Create Rendering namespace with RenderEngine and function declarations - remove this line when done"

// TODO: Create Rendering namespace

// TODO: Create RenderEngine struct with the following members:
// - GLFWwindow* window
// - std::vector<glm::vec2> stars
// - std::vector<Ray>* rays (pointer to external rays vector)
// - int frame

// TODO: Add RenderEngine constructor
// - Parameters: int width, int height, const char* title, std::vector<Ray>& raysRef
// - Initializes GLFW/GLEW, creates window, generates stars

// TODO: Add RenderEngine destructor
// - Cleans up GLFW resources

// TODO: Add RenderEngine methods:
// - void beginFrame(float viewWidth, float viewHeight) - Setup projection and clear screen
// - void updatePhysics() - Integrate all rays using adaptive budget-based stepping
// - void drawFrame() - Draw entire scene (stars, black hole, rays, point source)
// - void endFrame() - Swap buffers and poll events
// - bool shouldClose() const - Check if window should close

// TODO: Declare generateStars function
// - Returns vector<glm::vec2>
// - Takes int count, float viewWidth, float viewHeight

// TODO: Declare drawStars function
// - Takes const vector<glm::vec2>& parameter

// TODO: Declare drawCircle function
// - Takes x, y, radius, segments

// TODO: Declare drawCircleOutline function
// - Takes x, y, radius, segments

// TODO: Declare drawDashedCircle function
// - Takes x, y, radius, segments

// TODO: Declare drawRays function
// - Takes const vector<Ray>& and currentFrame

// TODO: Declare drawPointSource function
// - Takes x, y coordinates
