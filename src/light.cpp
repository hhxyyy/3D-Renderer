#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>


// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color, float sample)
{
    float N = 10.0f;
    glm::vec3 segmentDirection = glm::normalize(segmentLight.endpoint1 - segmentLight.endpoint0);
    float segmentLength = glm::length(segmentLight.endpoint1 - segmentLight.endpoint0) / N;
    
    if (sample == 0.0f) {
        position = segmentLight.endpoint0;
        color = segmentLight.color0;
        return;
    }

    if (sample > N) {
        position = segmentLight.endpoint1;
        color = segmentLight.color1;
        return;
    } else {
        float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        position = segmentLight.endpoint0 + (sample - 1.0f) * segmentLength * segmentDirection + r * segmentDirection * segmentLength;
        float ratio = glm::length(segmentLight.endpoint1 - position) / glm::length(segmentLight.endpoint1 - segmentLight.endpoint0);
        color = (ratio * segmentLight.color0 + (1 - ratio) * segmentLight.color1);
    }
}


// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color, float rowInd, float columnInd)
{
    float N = 10.0f;
    float rh;
    float rv;
    glm::vec3 horizontal(0.0f);
    glm::vec3 vertical(0.0f);
    float horizontalSampleLength = glm::length(parallelogramLight.edge02) / N;
    float verticalSampleLength = glm::length(parallelogramLight.edge01) / N;
    glm::vec3 horizontalColor1(0.0f);
    glm::vec3 horizontalColor2(0.0f);
    if (rowInd == 0.0f) {
        horizontal = glm::vec3 {0.0f, 0.0f, 0.0f};
        horizontalColor1 = parallelogramLight.color0;
        horizontalColor2 = parallelogramLight.color1;
    } else {
        if (rowInd >= N) {
            horizontal = parallelogramLight.edge02;
            horizontalColor1 = parallelogramLight.color2;
            horizontalColor2 = parallelogramLight.color3;
        } else {
            rh = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
            horizontal = (rowInd - 1.0f) * horizontalSampleLength * glm::normalize(parallelogramLight.edge02) + rh * horizontalSampleLength * glm::normalize(parallelogramLight.edge02);
            float ratio = glm::length( horizontal) / glm::length(parallelogramLight.edge02);
            horizontalColor1 = ratio * parallelogramLight.color2 + (1 - ratio) * parallelogramLight.color0;
            horizontalColor2 = ratio * parallelogramLight.color3 + (1 - ratio) * parallelogramLight.color1;
        }
    }
    if (columnInd == 0.0f) {
        vertical = glm::vec3 { 0.0f, 0.0f, 0.0f };
    } else {
        if (columnInd >= N) {
            vertical = parallelogramLight.edge01;
        } else {
            rv = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
            vertical = (columnInd - 1.0f) * verticalSampleLength * glm::normalize(parallelogramLight.edge01) + rv * verticalSampleLength * glm::normalize(parallelogramLight.edge01);
        }
    }
    position = parallelogramLight.v0 + horizontal + vertical;
    float ratio = glm::length( vertical) / glm::length(parallelogramLight.edge01);
    color = ratio * horizontalColor2 + (1 - ratio) * horizontalColor1;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement this function.
    // sld refers to the direction from sample position to light position
    if (!features.enableHardShadow)
        return 1.0f;
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 lightVector = samplePos - intersection;
    intersection = intersection + 0.000001f * glm::normalize(lightVector);
    Ray r1 { intersection, lightVector, 1.0f };
    if (bvh.intersect(r1, hitInfo, features) && 1.0f - r1.t > 1e-9) {
        drawRay(Ray { intersection, lightVector, 1.0f }, 1.0f * debugColor);
        return 0.0f;
    }
    r1.t = 1.0f;
    drawRay(Ray { intersection, lightVector, 1.0f }, glm::vec3 { 1.0f, 1.0f, 1.0f });
    return 1.0f;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 finalShadow { 0.0f, 0.0f, 0.0f };
        glm::vec3 intersection = ray.origin + ray.t * ray.direction;
        
        // TODO: replace this by your own implementation of shading
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                // Perform your calculations for a point light.
                glm::vec3 colorFromThisLight = computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                finalShadow += colorFromThisLight * testVisibilityLightSample(pointLight.position, colorFromThisLight, bvh, features, ray, hitInfo);
            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                // Perform your calculations for a segment light
                if (features.enableSoftShadow) {
                    glm::vec3 position = glm::vec3(0.0);
                    glm::vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);
                    float counter = 0.0f;
                    
                    while (position != segmentLight.endpoint1) {
                        sampleSegmentLight(segmentLight, position, color, counter);
                        glm::vec3 shading = computeShading(position, color, features, ray, hitInfo);
                        finalShadow += shading * testVisibilityLightSample(position, shading, bvh, features, ray, hitInfo);
                        counter += 1.0f;
                        
                        // debug
                        glm::vec3 L = glm::normalize(position - intersection);
                        float rayLength = glm::length(position - intersection);
                        intersection = intersection + 0.000001f * glm::normalize(L);
                        Ray r1 { intersection, L, rayLength };
                        
                        HitInfo h = hitInfo;
                        if (bvh.intersect(r1, hitInfo, features)) {
                            glm::vec3 direction = glm::vec3((r1.origin + r1.direction * r1.t) - position);
                            drawRay(r1, glm::vec3(1.0f, 0.0f, 0.0f));
                        } else {
                            drawRay(Ray { intersection, L, rayLength }, color);
                        }
                        hitInfo = h;
                    }
                    finalShadow = finalShadow / counter;
                } 
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                // Perform your calculations for a parallelogram light.
                if (features.enableSoftShadow) {
                    glm::vec3 position = glm::vec3(0.0);
                    glm::vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);
                    float N = 10.0f;
                    
                    for (float columnIndex = 0.0f; columnIndex <= N; columnIndex += 1.0f) {
                        for (float rowIndex = 0.0f; rowIndex <= N; rowIndex += 1.0f) {
                            sampleParallelogramLight(parallelogramLight, position, color, rowIndex, columnIndex);
                            glm::vec3 shading = computeShading(position, color, features, ray, hitInfo);
                            finalShadow += shading * testVisibilityLightSample(position, shading, bvh, features, ray, hitInfo);

                            // debug
                            glm::vec3 L = glm::normalize(position - intersection);
                            float rayLength = glm::length(position - intersection);
                            intersection = intersection + 0.000001f * glm::normalize(L);
                            Ray r1 { intersection, L, rayLength };
                            HitInfo h = hitInfo;
                            if (bvh.intersect(r1, hitInfo, features)) {
                                drawRay(r1, glm::vec3(1.0f, 0.0f, 0.0f));
                            } else {
                                drawRay(Ray { intersection, L, rayLength }, color);
                            }
                            hitInfo = h;
                        }
                    }
                    finalShadow = finalShadow / (N * N);
                }
            }
            return finalShadow;
        }

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
