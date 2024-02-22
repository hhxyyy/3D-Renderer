#include "texture.h"
#include <framework/image.h>
#include <iostream>
#include "render.h"

// Makes sure that t is within [0, 1]
// Applies repeat mode for a texture
float repeatMode(float t) {
    if (t < 0) {
        t += std::floor(std::abs(t));
        t++;
    } else
        t -= std::floor(std::abs(t));
    return t;
}

// Linear interpolation of two colors based on the distance
glm::vec3 linearInterpolation(glm::vec3 v0, glm::vec3 v1, float distance)
{
    return ((1 - distance) * v0) + (distance * v1);
}

// Bilinear interpolation of a texCoord
glm::vec3 bilinearInterpolation(const Image& image, const glm::vec2& texCoord) {
    float i = repeatMode(texCoord.x) * image.width - 0.5f;
    float j = (1 - repeatMode(texCoord.y)) * image.height - 0.5f;
    // If i or j is smaller than 0, apply repeat mode
    if (i < 0)
        i += image.width;
    if (j < 0)
        j += image.height;
    int xLeft = std::floor(i);
    int yTop = std::floor(j);
    int xRight = (xLeft + 1) % image.width;
    int yBottom = (yTop + 1) % image.height;
    float alpha = i - xLeft;
    float beta = j - yTop;

    // Get all four vertices 
    glm::vec3 v0 = image.pixels[yTop * image.width + xLeft];
    glm::vec3 v1 = image.pixels[yTop * image.width + xRight];
    glm::vec3 v2 = image.pixels[yBottom * image.width + xLeft];
    glm::vec3 v3 = image.pixels[yBottom * image.width + xRight];

    glm::vec3 colorTop = linearInterpolation(v0, v1, alpha);
    glm::vec3 colorBottom = linearInterpolation(v2, v3, alpha);
    
    return linearInterpolation(colorTop, colorBottom, beta);
}

// Trilinear interpolation for a given level
glm::vec3 trilinearInterpolation(const Mesh& mesh, const Image& image, const glm::vec2& texCoord, float level) {
    int lower = std::floor(level);
    int higher = lower + 1;
    float alpha = level - lower;
    // Color at top for lower mipmap
    glm::vec3 v0 = bilinearInterpolation(mesh.mipmaps[lower], texCoord);
    // Color at bottom for upper mipmap
    glm::vec3 v1 = bilinearInterpolation(mesh.mipmaps[higher], texCoord);
    // Interpolate both colors
    return linearInterpolation(v0, v1, alpha);
}

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    if (features.extra.enableBilinearTextureFiltering) {
        return bilinearInterpolation(image, texCoord);
    } else {
        int i = repeatMode(texCoord.x) * image.width;
        int j = repeatMode(texCoord.y) * image.height;
        int ind = (image.height - 1 - j) * image.width + i;
        return image.pixels[ind];
    }
}