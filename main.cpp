#include <cmath>
#include <iostream>

// Project headers
#include "constants.h"
#include "ray.h"
#include "rendering.h"

void generateOrbitingRay(std::vector<Ray>& rays) {
    const float orbitStartX{-0.9f * Visual::VIEW_WIDTH};
    const float orbitStartY{2.577934f};  // Just below critical 2.598 RS
    rays.emplace_back(orbitStartX, orbitStartY, 1.0f, 0.0f,
                      RayScenario::ORBITING, Visual::ORBITING_START);
}

void generatePointSourceRays(std::vector<Ray>& rays) {
    const int numPointRays{25};
    const float toBlackHoleX{0.0f - Visual::POINT_SOURCE_X};
    const float toBlackHoleY{0.0f - Visual::POINT_SOURCE_Y};
    const float baseAngle{std::atan2(toBlackHoleY, toBlackHoleX)};
    const float spreadAngle{60.0f * static_cast<float>(M_PI) / 180.0f};

    for (int i{0}; i < numPointRays; ++i) {
        const float angleOffset{-spreadAngle / 2.0f + spreadAngle * static_cast<float>(i) / static_cast<float>(numPointRays - 1)};
        const float angle{baseAngle + angleOffset};
        const float vx{std::cos(angle)};
        const float vy{std::sin(angle)};
        rays.emplace_back(Visual::POINT_SOURCE_X, Visual::POINT_SOURCE_Y, vx, vy,
                          RayScenario::POINT_SOURCE, Visual::POINT_SOURCE_START);
    }
}

void generateParallelRays(std::vector<Ray>& rays) {
    const int numParallelRays{70};
    const float parallelStartX{-Visual::VIEW_WIDTH};

    for (int i{0}; i < numParallelRays; ++i) {
        const float t{static_cast<float>(i) / static_cast<float>(numParallelRays - 1)};
        const float startY{-Visual::VIEW_HEIGHT + t * 2.0f * Visual::VIEW_HEIGHT};
        rays.emplace_back(parallelStartX, startY, 1.0f, 0.0f,
                          RayScenario::PARALLEL, Visual::PARALLEL_START);
    }
}

std::vector<Ray> generateRays() {
    std::vector<Ray> rays;

    generateOrbitingRay(rays);
    generatePointSourceRays(rays);
    generateParallelRays(rays);

    std::cout << "Total rays: " << rays.size() << "\n";
    std::cout << "  - Orbiting: 1\n";
    std::cout << "  - Point source: 25\n";
    std::cout << "  - Parallel: 70\n";

    return rays;
}

int main() {
    std::cout << "\n=== Black Hole Simulation (Geometric Units) ===\n";
    std::cout << "Schwarzschild radius: " << BlackHole::RS << " (RS = 1)\n\n";

    std::vector<Ray> rays{generateRays()};

    Rendering::RenderEngine engine{
        Visual::WINDOW_WIDTH,
        Visual::WINDOW_HEIGHT,
        "2D Black Hole Simulator - Organized Project",
        rays
    };

    while (!engine.shouldClose()) {
        engine.beginFrame(Visual::VIEW_WIDTH, Visual::VIEW_HEIGHT);
        engine.updatePhysics();
        engine.drawFrame();
        engine.endFrame();
    }

    return 0;
}