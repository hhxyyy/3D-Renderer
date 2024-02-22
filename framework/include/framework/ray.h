#pragma once
#include "disable_all_warnings.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <limits>
#include <vector>

struct Ray {
    glm::vec3 origin { 0.0f };
    glm::vec3 direction { 0.0f, 0.0f, -1.0f };
    float t { std::numeric_limits<float>::max() };
    // Index of a mesh that the ray hits
    int meshIndex {-1};
    // List of adjacent rays for a ray (max 4)
    std::vector<Ray> adjacentRays;
};
