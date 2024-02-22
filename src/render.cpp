#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "bloom.h"
#include "texture.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>
#include <stb/stb_vorbis.c>

int motionSampleNumber = 100;
glm::vec3 motionDistance = { 0.3f, 0.0f, 0.0f };
int filterSize = 0;

float focal = 3.0f;
float aperture = 0.1f;

static const std::filesystem::path dataDirPath { DATA_DIR };
Image posx(dataDirPath / "uwp.jpeg");

std::vector<glm::vec3> filtered;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {
        // We need to implement mipmapping here, and not at bvh.intersect, as getLevel calls bvh.intersect, so we would end up with infinite loop
        if (features.enableTextureMapping && hitInfo.material.kdTexture) {
            // Update hitInfo.material.kd with appropriate function
            if (features.extra.enableMipmapTextureFiltering) {
                Image image = scene.meshes[ray.meshIndex].mipmaps[0];
                float level = getLevel(bvh, hitInfo, ray, glm::vec2 { image.width, image.height }, features);
                hitInfo.material.kd = trilinearInterpolation(scene.meshes[ray.meshIndex], *hitInfo.material.kdTexture, hitInfo.texCoord, level);
            } else
                hitInfo.material.kd = acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
        }

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            // TODO: put your own implementation of recursive ray tracing here.
            // According to chapter 4.8 in the book
            // Ray may bounce forever, therefore we add contraint of a maximim recursion depth
            if (glm::length(hitInfo.material.ks) > 0 && rayDepth <= 3) {
                Lo += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
            }
        }

        // Draw a white debug ray if the ray hits.
        drawRay(ray, Lo);

        // Debugger for multiple rays per pixels
        if (features.extra.enableDebugMultipleRays) {
            drawMRDebug(ray, bvh, features);
        }
        
        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        glm::vec3 color = { 1.0f, 0.0f, 0.0f };

        if (features.extra.enableEnvironmentMapping) {
            glm::vec3 d = ray.direction;
            float v = ((M_PI - acosf(-d.y)) / (M_PI));
            float u = ((M_PI + atan2(d.x, d.z)) / (2.0f * M_PI));
            int xu = (int)(u * posx.width);
            int xv = (int)(v * posx.height);
            int index = xv * posx.width + xu;
            color = posx.pixels[index];
            drawRay(ray, color);
            return color;
        }

        // Draw a red debug ray if the ray missed.
        drawRay(ray, color);
        // Set the color of the pixel to black if the ray misses/background color if enable environment map

        
        return { 0.0f, 0.0f, 0.0f };
        
    }
}


void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            Ray cameraRay = camera.generateRay(normalizedPixelPos);
            // It has to be here before setting the pixel
            if (features.enableTextureMapping && features.extra.enableMipmapTextureFiltering) {
                // We need to generate adjacent rays for a given ray here, as here we have access to the camera
                generateAdjacentRays(camera, cameraRay, glm::vec2 { x, y }, windowResolution);
            }
            glm::vec3 finalColor = getFinalColor(scene, bvh, cameraRay, features);

            if (features.extra.enableMotionBlur) {             
                glm::vec3 color = finalColor;
                glm::vec3 originColor = color;
                color = motionBlur(scene, cameraRay, bvh, features, color, motionDistance, motionSampleNumber);
                color = originColor + color / (float)motionSampleNumber * 2.0f;
                finalColor = color;
            } 
            if (features.extra.enableMultipleRaysPerPixel) {
                finalColor = multipleRays(x, y, scene, camera, bvh, features, windowResolution);
            }
            if (features.extra.enableDepthOfField) {
                finalColor = dof(scene, bvh, features, cameraRay, focal, aperture);
            }

            screen.setPixel(x, y, finalColor);
        }
    }
    if (features.extra.enableBloomEffect) {
        /*
        std::vector<glm::vec3>  filtered = bloomFilter(screen, 5);
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                int index = screen.indexAt(x, y);
                screen.setPixel(x, y, screen.pixels()[index] + filtered[index]);
                //screen.setPixel(x, y, filtered[index]);
            }
        }*/
        filtered = bloomFilter(screen, filterSize, 1);
    }
}

glm::vec3 multipleRays(int x, int y, const Scene& scene, const Trackball& camera, const BvhInterface& bvh, const Features& features, glm::ivec2 windowResolution)
{
    int samples = 20;
    glm::vec2 pixel;
    glm::vec3 newColor(0.0f);

    for (int i = 0; i < samples; i++) {
        for (int j = 0; j < samples; j++) {
            float randI = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
            float randJ = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

            pixel.x = (float)(float(x) + (i + randI) / (samples * 1.0f)) / float(windowResolution.x) * 2.0f - 1.0f;
            pixel.y = (float)(float(y) + (j + randJ) / (samples * 1.0f)) / float(windowResolution.y) * 2.0f - 1.0f;

            Ray newRay = camera.generateRay(pixel);
            newColor += getFinalColor(scene, bvh, newRay, features);
        }
    }
    glm::vec3 finalColor = newColor / (samples * samples * 1.0f);
    return finalColor;
}

void drawMRDebug(Ray ray, const BvhInterface& bvh, const Features& features)
{
    HitInfo hitInfo;
    int samples = 5;
    for (int i = 0; i < samples; ++i) {
        for (int j = 0; j < samples; ++j) {
            float randX = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f - 1.0f) / 100.0f;
            float randY = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f - 1.0f) / 100.0f;
            float randZ = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 2.0f - 1.0f) / 100.0f;

            glm::vec3 rayDirection = glm::normalize(glm::vec3(ray.direction.x + randX, ray.direction.y + randY, ray.direction.z + randZ));
            Ray newRay { ray.origin, rayDirection, std::numeric_limits<float>::max() - 1.0f };

            HitInfo h = hitInfo;
            bvh.intersect(newRay, hitInfo, features);
            drawRay(newRay, glm::vec3 { 1.0f, 0.0f, 1.0f });
            hitInfo = h;
        }
    }
}

glm::vec3 dof(const Scene& scene, const BvhInterface& bvh, const Features& features, const Ray& cameraRay, float focal, float apertue) {
    float ffcoal = (float)focal;
    float fapertue = (float)apertue;
    
    glm::vec3 centrolPoint = cameraRay.origin + ffcoal * cameraRay.direction;
    glm::vec3 color = { 0.0f, 0.0f, 0.0f };
    for (int i = 0; i < 100; i++) {
        float offsetx = std::rand() / double(RAND_MAX) * 2.0f - 1.0f;
        float shiftx = offsetx * fapertue;
        float offsety = std::rand() / double(RAND_MAX) * 2.0f - 1.0f;
        float shifty = offsety * fapertue;
        float offsetz = std::rand() / double(RAND_MAX) * 2.0f - 1.0f;
        float shiftz = offsetz * fapertue;
        glm::vec3 secondOrigin = cameraRay.origin + glm::vec3 { shiftx, shifty, shiftz };
        glm::vec3 secondDirection = centrolPoint - secondOrigin;
        secondDirection = glm::normalize(secondDirection);
        Ray secondRay { secondOrigin, secondDirection };
        color = color + getFinalColor(scene, bvh, secondRay, features);
    }
    color = color / 100.0f;
    return color;
}

glm::vec3 motionBlur(const Scene& scene, const Ray& cameraRay, const BvhInterface& bvh, const Features& features, const glm::vec3 color, glm::vec3 distance, int divide)
{
    glm::vec3 cc = color;

    for (int i = 0; i < divide; i++) {
        float offset = std::rand() / double(RAND_MAX);
        float shiftx = offset * distance.x;
        float shifty = offset * distance.y;
        float shiftz = offset * distance.z;
        glm::vec3 o;
        /*
        if (i < (int)(0.1f * divide)) {
             o = { cameraRay.origin.x - shift, cameraRay.origin.y, cameraRay.origin.z };
        } else {
             o = { cameraRay.origin.x + shift, cameraRay.origin.y, cameraRay.origin.z };
        }*/
        o = { cameraRay.origin.x + shiftx, cameraRay.origin.y + shifty, cameraRay.origin.z + shiftz };
        Ray r1 { o, cameraRay.direction, cameraRay.t };
        glm::vec3 scolor = getFinalColor(scene, bvh, r1, features);
        cc = cc + (1.0f - offset) * scolor;
    }
    return cc;

}

void setDivision(int newDivide) {
    motionSampleNumber = newDivide;
    //std::cout << motionSampleNumber << std::endl;
}

void setMotionDistance(glm::vec3 distance) {
    motionDistance = distance;
}

void setBloomFilterSize(int size) {
    filterSize = size;
}

void debugMotion(glm::vec3 direction) {
    Ray r1 { { 0.0f, 0.0f, 0.0f }, direction, 2.0f };
    drawRay(r1, { 1.0f, 0.0f, 0.0f });
}

void debugBloom(Screen& screen) {
    glm::ivec2 windowResolution = screen.resolution();

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            int colorIndex = screen.indexAt(x, y);
            screen.setPixel(x, y, filtered[colorIndex]);
        }
    }
}

void debugDOF(Screen& screen, const Trackball& camera) {
    glm::ivec2 windowResolution = screen.resolution();
    const glm::vec2 normalizedPixelPos {
        float(0) / float(windowResolution.x) * 2.0f - 1.0f,
        float(0) / float(windowResolution.y) * 2.0f - 1.0f
    };
    Ray cameraRay = camera.generateRay(normalizedPixelPos);

    float ffcoal = 10.0f;
    float fapertue = (float)aperture;
    Ray r1 { cameraRay.origin, cameraRay.direction, ffcoal };
    drawRay(r1, { 1.0f, 0.0f, 0.0f });


    glm::vec3 centrolPoint = cameraRay.origin + ffcoal * cameraRay.direction;
    glm::vec3 color = { 1.0f, 0.0f, 0.0f };
    
        float offsetx = 0.5f;
        float shiftx = offsetx * fapertue;
        float offsety = 0.5f;
        float shifty = offsety * fapertue;
        float offsetz = 0.5f;
        float shiftz = offsetz * fapertue;
        glm::vec3 secondOrigin = cameraRay.origin + glm::vec3 { shiftx, shifty, shiftz };
        glm::vec3 secondDirection = centrolPoint - secondOrigin;
        secondDirection = glm::normalize(secondDirection);
        Ray secondRay { secondOrigin, secondDirection };
        drawRay(secondRay, { 1.0f, 1.0f, 1.0f });

}

void setFocal(float f) {
    focal = f;
}

void setAperture(float a) {
    aperture = a;
}


    /*
void motionBlur(const Scene& scene, Screen& screen, const Trackball& camera, const BvhInterface& bvh, const Features& features, float distance, int divide)
{
    glm::ivec2 windowResolution = screen.resolution();
    std::vector<glm::vec3> originColor = screen.pixels();
    std::vector<glm::vec3> finalColor = screen.pixels();

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            float offset = distance / (1.0f * divide);
            float lightDecrease = 0.9f;
            int index = screen.indexAt(x, y);
            for (int i = 1; i < divide; i++) {
                const glm::vec2 normalizedPixelPos {
                    (float(x) + offset) / float(windowResolution.x) * 2.0f - 1.0f,
                    (float(y)) / float(windowResolution.y) * 2.0f - 1.0f
                };
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                glm::vec3 color = getFinalColor(scene, bvh, cameraRay, features);

                finalColor[index] = finalColor[index] + lightDecrease * color;
                offset += offset;
                lightDecrease = lightDecrease - 0.9f / (1.0f * divide);
            }
            finalColor[index] = finalColor[index];
            screen.setPixel(x, y, finalColor[index]);
        }
    }
}*/



glm::vec2 getNormalizedPixel(glm::vec2& pos, const glm::ivec2& windowResolution) {
    return glm::vec2 {
        float(pos.x) / float(windowResolution.x) * 2.0f - 1.0f,
        float(pos.y) / float(windowResolution.y) * 2.0f - 1.0f
    };
}


// Checks if a position is on the screen
bool inRange(const glm::vec2& pos, const glm::ivec2& res)
{
    return pos.x > -1 && pos.y > -1 && pos.x < res.x && pos.y < res.y;
}

// Generates adjacent rays in all four direction
void generateAdjacentRays(const Trackball& camera, Ray& ray, const glm::vec2& position, const glm::ivec2& resolution) {
    std::array direction { -1, 0, 1, 0 };
    // For all directions it generates a ray if it fits on the screen
    for (int i = 0; i < direction.size(); i++) {
        int dirX = direction[i];
        int dirY = direction[(i + 1) % direction.size()];
        int x = position.x + dirX;
        int y = position.y + dirY;
        glm::vec2 newPixel { x, y };
        if (inRange(newPixel, resolution)) {
            ray.adjacentRays.push_back(camera.generateRay(getNormalizedPixel(newPixel, resolution)));
        }
    }
}

// Gets level of a MipMap based on the distance between texels
float getLevel(const BvhInterface& bvh, HitInfo& hitInfo, Ray& ray, const glm::vec2 imageDim, const Features& features)
{
    float d = 0.0f;
    for (Ray adjacentRay : ray.adjacentRays) {
        HitInfo hitAdjacent;
        // We first need to make sure that an adjacent ray actually hits something and it hits the same mesh
        if (bvh.intersect(adjacentRay, hitAdjacent, features) && adjacentRay.meshIndex == ray.meshIndex) {
            // Substract the texCoord and later rescale to the actual pixel dimensions
            float x = std::abs(hitInfo.texCoord.x - hitAdjacent.texCoord.x) * imageDim.x;
            float y = std::abs(hitInfo.texCoord.y - hitAdjacent.texCoord.y) * imageDim.y;
            float dis = x;
            if (y > dis)
                dis = y;
            if (dis > d)
                d = dis;
        }
    }
    if (d < 1.0f)
        return 0.0f;
    return std::log2(d);
}