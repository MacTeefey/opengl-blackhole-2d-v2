#pragma once

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Geometric units: c = G = 1, Schwarzschild radius = 1
namespace BlackHole {
    constexpr float RS{1.0f};  // Schwarzschild radius (geometric units)
}

// Visual configuration
namespace Visual {
    constexpr int NUM_STARS{200};
    constexpr int CIRCLE_SEGMENTS{100};

    // Window settings
    constexpr int WINDOW_WIDTH{800};
    constexpr int WINDOW_HEIGHT{600};

    // Viewport in geometric units (units of RS)
    constexpr float VIEW_WIDTH{8.0f};     // 8 Schwarzschild radii
    constexpr float VIEW_HEIGHT{6.0f};    // 6 Schwarzschild radii

    // Point source position (geometric units)
    constexpr float POINT_SOURCE_X{-0.95f * VIEW_WIDTH};
    constexpr float POINT_SOURCE_Y{0.85f * VIEW_HEIGHT};

    // Ray timing (frames at 60 FPS with VSync)
    constexpr int ORBITING_START{0};       // Immediately
    constexpr int POINT_SOURCE_START{180}; // ~3 seconds
    constexpr int PARALLEL_START{360};     // ~6 seconds
}

// Simulation parameters
namespace Simulation {
    constexpr float MAX_DISTANCE{30.0f};       // Maximum escape distance (30 RS)
    constexpr float BUDGET_PER_FRAME{0.08f};   // Integration budget per frame
    constexpr float MIN_STEP{0.01f};           // Minimum adaptive step size
    constexpr float MAX_STEP{0.1f};            // Maximum adaptive step size
}