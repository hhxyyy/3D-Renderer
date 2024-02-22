#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    // TODO: implement this function.
    glm::vec3 v01 = -v0 + v1;
    glm::vec3 v12 = -v1 + v2;
    glm::vec3 v20 = -v2 + v0;
    glm::vec3 v0p = p - v0;
    glm::vec3 v1p = p - v1;
    glm::vec3 v2p = p - v2;
    glm::vec3 c1 = glm::cross(v01, v0p);
    glm::vec3 c2 = glm::cross(v12, v1p);
    glm::vec3 c3 = glm::cross(v20, v2p);
    float a = glm::dot(c1, c2);
    float b = glm::dot(c2, c3);
    float c = glm::dot(c3, c1);
    if (a >= 0 && b >= 0 && c >= 0) {
        return true;
    } else if (a <= 0 && b <= 0 && c <= 0) {
        return true;
    }

    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    // TODO: implement this function.
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal);
    if (t <= 0.0f) {
        return false;
    }
    if (t > ray.t) {
        return false;
    }
    ray.t = t;
    return true;
    // return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    // TODO: implement this function.
    Plane plane;
    glm::vec3 v01 = v1 - v0;
    glm::vec3 v02 = v2 - v0;
    glm::vec3 n = glm::cross(v01, v02);
    n = glm::normalize(n);
    float d = glm::dot(n, v0);
    plane = { d, n };
    return plane;
    // return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    float oldT = ray.t;
    Plane plane = trianglePlane(v0, v1, v2);
    if (intersectRayWithPlane(plane, ray)) {
        glm::vec3 p = ray.origin + ray.t * ray.direction;
        if (pointInTriangle(v0, v1, v2, plane.normal, p)) {
            if (oldT < ray.t) {
                ray.t = oldT;
                return false;
            }
            return true;
        }
        ray.t = oldT;
        return false;
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    glm::vec3 morigin = ray.origin - sphere.center;
    float a = ray.direction.x * ray.direction.x + ray.direction.y * ray.direction.y + ray.direction.z * ray.direction.z;
    float b = 2 * (ray.direction.x * morigin.x + ray.direction.y * morigin.y + ray.direction.z * morigin.z);
    float c = morigin.x * morigin.x + morigin.y * morigin.y + morigin.z * morigin.z - sphere.radius * sphere.radius;

    float delta = b * b - 4 * a * c;

    if (delta < 0.0f) {
        return false;
    }

    float t;

    if (delta > 0.0f) {
        float t0 = (-b + sqrt(delta)) / (2 * a);
        float t1 = (-b - sqrt(delta)) / (2 * a);
        if (t0 > 0.0f && t1 > 0.0f)
            t = (t0 < t1) ? t0 : t1;
        else if (t0 < 0.0f && t1 < 0.0f) {
            return false;
        } else {
            if (t0 > 0.0f)
                t = t0;
            else if (t1 > 0.0f)
                t = t1;
            else
                return false;
        }
    } else if (delta == 0.0f) {
        t = -b / (2 * a);
        if (t <= 0.0f)
            return false;
    }

    if (t < ray.t) {
        ray.t = t;
        return true;
    }

    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    // TODO: implement this function.
    float txmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;

    float txmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float tymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;

    float tinx = (txmin < txmax) ? txmin : txmax;
    float tiny = (tymin < tymax) ? tymin : tymax;
    float tinz = (tzmin < tzmax) ? tzmin : tzmax;

    float toutx = (txmin > txmax) ? txmin : txmax;
    float touty = (tymin > tymax) ? tymin : tymax;
    float toutz = (tzmin > tzmax) ? tzmin : tzmax;

    float tin;
    if (tinx > tiny) {
        if (tinx > tinz) {
            tin = tinx;
        } else {
            tin = tinz;
        }
    } else {
        if (tiny > tinz) {
            tin = tiny;
        } else {
            tin = tinz;
        }
    }

    float tout;
    if (toutx < touty) {
        if (toutx < toutz) {
            tout = toutx;
        } else {
            tout = toutz;
        }
    } else {
        if (touty < toutz) {
            tout = touty;
        } else {
            tout = toutz;
        }
    }

    if (tout < 0.0f || tin > tout) {
        return false;
    }

    if (tin < ray.t && tin >= 0.0f) {
        ray.t = tin;
        return true;
    }
    return false;
}
