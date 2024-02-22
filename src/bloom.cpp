#include "bloom.h"
#include "intersect.h"
#include "light.h"
#include "render.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

std::vector<glm::vec3> bloomFilter(Screen& screen, int filterSize)
{
    Filter filter { filterSize };
    glm::ivec2 windowResolution = screen.resolution();
    float weight[5] = { 0.227027f, 0.1945946f, 0.1216216f, 0.054054f, 0.016216f };

    std::vector<glm::vec3> originColors = screen.pixels();
    std::vector<glm::vec3> threshold;
    for (auto& color : originColors) {
        float gamma = color.x * 0.2126f + color.y * 0.7152f + color.z * 0.0722f;
        if (gamma >= 0.6f) {
            threshold.push_back(color);
        } else {
            threshold.push_back({ 0.0f, 0.0f, 0.0f });
        }
    }

    std::vector<glm::vec3> filteredColors(threshold.size());

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            
            glm::vec3 sum = { 0.0f, 0.0f, 0.0f };
            int originIndex = screen.indexAt(x, y);
            sum = sum + threshold[originIndex] * weight[0];

            for (int fx = 1; fx < filter.size; fx++) {
                int xx1 = x + fx;
                int xx2 = x - fx;
                glm::vec3 c1;
                glm::vec3 c2;
                if (xx1 < 0 || xx1 >= windowResolution.x) {
                    c1 = { 0.0f, 0.0f, 0.0f };
                    sum = sum + c1;
                } else {
                    int newIndex = screen.indexAt(xx1, y);
                    sum = sum + threshold[newIndex] * weight[fx];
                }
                if (xx2 < 0 || xx2 >= windowResolution.x) {
                    c2 = { 0.0f, 0.0f, 0.0f };
                    sum = sum + c2;
                } else {
                    int newIndex = screen.indexAt(xx2, y);
                    sum = sum + threshold[newIndex] * weight[fx];
                }
            }
            threshold[originIndex] = sum;

            

           
        }
    }

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            glm::vec3 sum = { 0.0f, 0.0f, 0.0f };
            int originIndex = screen.indexAt(x, y);
            sum = sum + threshold[originIndex] * weight[0];

            for (int fy = 1; fy < filter.size; fy++) {
                int yy1 = y + fy;
                int yy2 = y - fy;
                glm::vec3 c1;
                glm::vec3 c2;
                if (yy1 < 0 || yy1 >= windowResolution.y) {
                    c1 = { 0.0f, 0.0f, 0.0f };
                    sum = sum + c1;
                } else {
                    int newIndex = screen.indexAt(x, yy1);
                    sum = sum + threshold[newIndex] * weight[fy];
                }
                if (yy2 < 0 || yy2 >= windowResolution.y) {
                    c2 = { 0.0f, 0.0f, 0.0f };
                    sum = sum + c2;
                } else {
                    int newIndex = screen.indexAt(x, yy2);
                    sum = sum + threshold[newIndex] * weight[fy];
                }
            }

            threshold[originIndex] = sum;
        }
    }

    return threshold;

}



std::vector<glm::vec3> bloomFilter(Screen& screen, int filterSize, int type)
{
    Filter filter { filterSize };
    glm::ivec2 windowResolution = screen.resolution();
    std::vector<glm::vec3> filtered(screen.pixels().size());
    std::vector<glm::vec3> originColors = screen.pixels();
    std::vector<glm::vec3> threshold;
    for (auto& color : originColors) {
        float gamma = color.x * 0.2126f + color.y * 0.7152f + color.z * 0.0722f;
        if (gamma >= 0.6f) {
            threshold.push_back(color);
        } else {
            threshold.push_back({ 0.0f, 0.0f, 0.0f });
        }
    }

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };

            glm::vec3 sum = { 0.0f, 0.0f, 0.0f };

            for (int fx = -filter.size; fx <= filter.size; ++fx) {
                for (int fy = -filter.size; fy <= filter.size; ++fy) {
                    int xx = x + fx;
                    int yy = y + fy;
                    float fill = filter.fill;
                    if (xx < 0 || yy < 0 || xx >= windowResolution.x || yy >= windowResolution.y) {
                        fill = 0.0f;
                        sum = sum + glm::vec3(fill, fill, fill);
                    } else {
                        int thresholdIndex = screen.indexAt(xx, yy);
                        glm::vec3 filterValue = fill * threshold[thresholdIndex];
                        sum = sum + filterValue;
                    }
                }
            }
            // sum = sum * filter.scale;
            int colorIndex = screen.indexAt(x, y);
            sum = sum  + originColors[colorIndex];
            screen.setPixel(x, y, sum);
            filtered[colorIndex] = sum;
        }
    }
    return filtered;
}

void bloomImage(Image& image, Screen& screen)
{
    Filter filter {};

    std::vector<glm::vec3> originColors = image.pixels;
    std::vector<glm::vec3> threshold;
    std::vector<glm::vec3> final(image.pixels.size());
    for (auto& color : originColors) {
        float gamma = color.x * 0.2126f + color.y * 0.7152f + color.z * 0.0722f;
        if (gamma >= 0.6f) {
            threshold.push_back(color);
        } else {
            threshold.push_back({ 0.0f, 0.0f, 0.0f });
        }
    }

    for (int y = 0; y < image.height; y++) {
        for (int x = 0; x != image.width; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.

            glm::vec3 sum = { 0.0f, 0.0f, 0.0f };

            for (int fx = -filter.size; fx <= filter.size; ++fx) {
                for (int fy = -filter.size; fy <= filter.size; ++fy) {
                    int xx = x + fx;
                    int yy = y + fy;
                    float fill = filter.fill;
                    if (xx < 0 || yy < 0 || xx >= image.width || yy >= image.height) {
                        fill = 0.0f;
                        sum = sum + glm::vec3(fill, fill, fill);
                    } else {
                        int thresholdIndex = (yy * image.width + xx);
                        glm::vec3 filterValue = fill * threshold[thresholdIndex];
                        sum = sum + filterValue;
                    }
                }
            }
            // sum = sum * filter.scale;
            int colorIndex = (y * image.width + x);
            sum = sum + originColors[colorIndex];
            // screen.setPixel(x, y, sum);
            final[colorIndex] = sum;
        }
    }

    for (int y = 0; y < image.height; y++) {
        for (int x = 0; x != image.width; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            // int index = (y * image.width + x);
            int index = screen.indexAt(x, y);
            screen.setPixel(x, y, final[index]);
        }
    }
}