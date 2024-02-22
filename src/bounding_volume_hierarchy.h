#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

struct Node {
    // 0 is interior, 1 is leaves
    int type;

    // AABB contains all triangles in the node
    AxisAlignedBox aabb;

    // the child contains the indices (from the nodeList) of children node if node type == 0, the x in uvec2 is left node , the y in uvec2 is right node
    // the child contains all triangles of the node if node type == 1, the x in uvec2 is the mesh index, the y in uvec2 is triangle index ( for example: a triangle(x, y) is the y-th contained in the x-th mesh£©
    std::vector<glm::uvec2> child;

    // the level of the node
    int level;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level, const Features& features);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

    glm::vec2 getSplitPointSAH(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes, const AxisAlignedBox& box, int bin);
    Node buildNode(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes, std::vector<Node>& nodes, int maxLevel, int splitAxis, bool sahBinning);
    float getSplitPoint(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes, int splitAxis);
    AxisAlignedBox getAABB(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes);
    bool traverseBVH(Ray& ray, HitInfo& hitInfo, const Features& features) const;

private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    int level = 10;
    Node headNode;
    // If SAH is enabled then creates normal BVH for debug
    Node headNodeBVH;

    std::vector<Node> nodeList;

    // If SAH is enabled then creates normal BVH for debug
    std::vector<Node> nodeListBVH;
};

// Set number of bins to the number chosen by a user
void setBinsNumber(int bins);