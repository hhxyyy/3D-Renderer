#include "intersect.h"
#include "light.h"
#include "render.h"
#include "screen.h"

struct Filter {

    int size = 10;
    float fill = 1.0f / ((2 * size + 1) * (2 * size + 1));
    float scale = 100.0f;
};

std::vector<glm::vec3> bloomFilter(Screen& screen, int filterSize);
void bloomImage(Image& image, Screen& screen);
std::vector<glm::vec3> bloomFilter(Screen& screen, int filterSize, int type);