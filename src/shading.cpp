#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.
    // Computes diffuse term according to Phong shading model
    glm::vec3 intersection = ray.origin + ray.t * ray.direction;
    glm::vec3 L = glm::normalize(lightPosition - intersection);
    glm::vec3 n = glm::normalize(hitInfo.normal);
    float dot = glm::dot(n, L);
    glm::vec3 diffuse;
    if (dot <= 0)
        diffuse = { 0, 0, 0 };
    else
        diffuse = lightColor * hitInfo.material.kd * dot;
    
    glm::vec3 V = glm::normalize(ray.origin - intersection);
    glm::vec3 R = glm::normalize(2 * glm::dot(L, n) * n - L);

    // Compute the glossy reflections if enabled
    if (features.extra.enableGlossyReflection) {
        glm::vec3 gloss(0.0f);
        glm::vec3 reflectedViewer = glm::normalize(2 * glm::dot(V, n) * n - V);
        float width = 1.0f / hitInfo.material.shininess;
        glm::vec3 t { 0.0f, 0.0f, 1.0f };
        if (glm::dot(reflectedViewer, t) == glm::length(reflectedViewer) * glm::length(t)) {
            t = glm::vec3(0.0f, 1.0f, 0.0f);
        }
        glm::vec3 u = glm::normalize(glm::cross(reflectedViewer, t));
        glm::vec3 v = glm::normalize(glm::cross(reflectedViewer, u));

        int samples = 20;
        for (int i = 0; i < samples; ++i) {

            float randU = -width / 2.0f + (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * width;
            float randV = -width / 2.0f + (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * width;
            glm::vec3 perturbedD = glm::normalize(reflectedViewer + randU * u + randV * v);
            if (glm::dot(perturbedD, hitInfo.normal) < 0.0f)
                continue;
            Ray pRay { intersection, perturbedD, 0.4f };
            drawRay(pRay, glm::vec3 { 1.0f, 0.0f, 0.4f });
            Ray rRay { intersection, reflectedViewer, 1.0f };
            drawRay(rRay, glm::vec3 { 1.0f, 0.0f, 0.0f });
            float dot = glm::dot(reflectedViewer, perturbedD);
            gloss += lightColor * hitInfo.material.ks * dot;
        }
        return diffuse + gloss / (samples * 1.0f);
    }

    // Computes specular term according to Phong shading model
    glm::vec3 specular;
    if (glm::dot(n, L) <= 0 || glm::dot(R, V) <= 0)
        specular = { 0, 0, 0 };
    else {
        dot = glm::dot(V, R);
        specular = lightColor * hitInfo.material.ks * glm::pow(dot, hitInfo.material.shininess);
    }

    // Returns sum of both
    return diffuse + specular;
}

const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    // TODO: implement the reflection ray computation.
    glm::vec3 origin = ray.origin + ray.t * ray.direction;
    glm::vec3 direction = glm::normalize(ray.direction - 2 * glm::dot(ray.direction, hitInfo.normal) * hitInfo.normal);
    origin += 0.0001f * direction;
    Ray reflectionRay { origin, direction };
    return reflectionRay;
}