#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>

// Forward declarations.
struct Scene;
class Screen;
class Trackball;
class BvhInterface;
struct Features;
struct HitInfo;
struct Image;

// Main rendering function.
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);

// Get the color of a ray.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth = 0);

// Generate adjacents rays for a given ray
void generateAdjacentRays(const Trackball& camera, Ray& ray, const glm::vec2& position, const glm::ivec2& resolution);

// Get level of MipMap
float getLevel(const BvhInterface& bvh, HitInfo& hitInfo, Ray& ray, const glm::vec2 imageDim, const Features& features);


void setDivision(int newDivide);
void setMotionDistance(glm::vec3 distance);
void setBloomFilterSize(int size);
glm::vec3 motionBlur(const Scene& scene, const Ray& cameraRay, const BvhInterface& bvh, const Features& features, const glm::vec3 color, glm::vec3 distance, int divide);
glm::vec3 multipleRays(int x, int y, const Scene& scene, const Trackball& camera, const BvhInterface& bvh, const Features& features, glm::ivec2 windowResolution);
void drawMRDebug(Ray ray, const BvhInterface& bvh, const Features& features);
void debugMotion(glm::vec3 direction);
void debugBloom(Screen& screen);


void debugDOF(Screen& screen, const Trackball& camera);
glm::vec3 dof(const Scene& scene, const BvhInterface& bvh, const Features& features, const Ray& cameraRay, float focal, float apertue);
void setAperture(float a);
void setFocal(float f);