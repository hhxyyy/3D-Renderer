#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <queue>
#include <iostream>

int number_bins = 0;

void setBinsNumber(int bins)
{
    number_bins = bins;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{
    // TODO: implement BVH construction.
    int maxLevel = 0;
    int splitAxis = 0;
    Scene scene = *m_pScene;
    std::vector<glm::uvec2> triangleIndex;
    for (int i = 0; i < scene.meshes.size(); i++) {
        Mesh mesh = scene.meshes[i];
        for (int j = 0; j < mesh.triangles.size(); j++) {
            triangleIndex.push_back({ i, j });
        }
    }
    // Creates BVH for normal BVH without SAH
    std::vector<Node> nodesBVH;
    Node headBVH = buildNode(pScene, triangleIndex, nodesBVH, maxLevel, splitAxis, false);
    nodesBVH.push_back(headBVH);
    // If SAH is enabled, create BVH with SAH and set headNode to nodeSAH
    // For a visual debug keep a list of nodes without SAH enabled
    if (features.extra.enableBvhSahBinning) {
        std::vector<Node> nodesSAH;
        Node headSAH = buildNode(pScene, triangleIndex, nodesSAH, maxLevel, splitAxis, true);
        nodesSAH.push_back(headSAH);
        headNode = headSAH;
        nodeList = nodesSAH;
        
        headNodeBVH = headBVH;
        nodeListBVH = nodesBVH;
    } else {
        headNode = headBVH;
        nodeList = nodesBVH;
    }
}

// build nodes recursively
// splitAxis 0 is x, 1 is y, 2 is z
Node BoundingVolumeHierarchy::buildNode(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes, std::vector<Node>& nodes, int maxLevel, int splitAxis, bool sahBinning)
{
    if (maxLevel < level && triangleIndexes.size() > 1) {
        int type = 0; // interior
        AxisAlignedBox aabb = getAABB(pScene, triangleIndexes);
        float centroid = -1;
        // SAH changes splitAxis so remember it, if SAH is not enabled
        int oldSplitAxis = splitAxis;
        // In both cases getSplitPoint
        if (sahBinning) {
            glm::vec2 res = getSplitPointSAH(pScene, triangleIndexes, aabb, number_bins);
            centroid = res.x;
            splitAxis = res.y;
        } else
            centroid = getSplitPoint(pScene, triangleIndexes, splitAxis);
        std::vector<glm::uvec2> left;
        std::vector<glm::uvec2> right;
        Scene scene = *pScene;
        for (auto& triangle : triangleIndexes) {
            int meshIndex = triangle.x;
            int triangleIndex = triangle.y;
            glm::uvec3& triangles = scene.meshes[meshIndex].triangles[triangleIndex];
            glm::vec3 v1 = scene.meshes[meshIndex].vertices[triangles.x].position;
            glm::vec3 v2 = scene.meshes[meshIndex].vertices[triangles.y].position;
            glm::vec3 v3 = scene.meshes[meshIndex].vertices[triangles.z].position;
            float baryPoint = (v1[splitAxis] + v2[splitAxis] + v3[splitAxis]) / 3.0f;
            if (baryPoint <= centroid)
                left.push_back(triangle);
            else
                right.push_back(triangle);
        }
        splitAxis = oldSplitAxis;

        splitAxis = (splitAxis + 1) % 3;
        maxLevel++;
        Node lnode = buildNode(pScene, left, nodes, maxLevel, splitAxis, sahBinning);
        Node rnode = buildNode(pScene, right, nodes, maxLevel, splitAxis, sahBinning);

        std::vector<glm::uvec2> child;
        glm::uvec2 childIndex;
        nodes.push_back(lnode);
        childIndex.x = nodes.size() - 1;
        nodes.push_back(rnode);
        childIndex.y = nodes.size() - 1;
        child.push_back(childIndex);

        Node node { type, aabb, child, maxLevel - 1 };
        // nodes.push_back(node);
        return node;
    }

    int type = 1; // leaves
    AxisAlignedBox aabb = getAABB(pScene, triangleIndexes);
    std::vector<glm::uvec2> child = triangleIndexes;
    Node node { type, aabb, child, maxLevel };
    // nodes.push_back(node);
    return node;
}

// Calculates of the surface of the given box
float getSurface(const AxisAlignedBox& box) {
    float front = std::abs(box.upper[1] - box.lower[1]) * std::abs(box.upper[0] - box.lower[0]);
    float up = std::abs(box.upper[2] - box.lower[2]) * std::abs(box.upper[0] - box.lower[0]);
    float right = std::abs(box.upper[2] - box.lower[2]) * std::abs(box.upper[1] - box.lower[1]);
    return front * 2 + up * 2 + right * 2;
}

// Calculates cost of SAH splitting
float calculateCost(const AxisAlignedBox& leftBox, int leftTriangles, const AxisAlignedBox& rightBox, int rightTriangles)
{
    float leftSurface = getSurface(leftBox);
    float rightSurface = getSurface(rightBox);
    return leftSurface * leftTriangles + rightSurface * rightTriangles;
}

// Get split point according to SAH + Binning
glm::vec2 BoundingVolumeHierarchy::getSplitPointSAH(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes, const AxisAlignedBox& box, int bin)
{
    Scene scene = *pScene;
    float minimalCost { std::numeric_limits<float>::max() };
    float centroidRes = -1.0f;
    int splitAxisRes = -1;
    // Check all axis
    for (int splitAxis = 0; splitAxis < 2; splitAxis++) {
        float distance = std::abs(box.upper[splitAxis] - box.lower[splitAxis]);
        // Computes the width of singular bin
        float singularBin = distance / bin;
        for (int i = 1; i < bin; i++) {
            std::vector<glm::uvec2> left;
            std::vector<glm::uvec2> right;
            // Get centroid point
            float centroid = box.lower[splitAxis] + i * singularBin;
            for (auto& triangle : triangleIndexes) {
                int meshIndex = triangle.x;
                int triangleIndex = triangle.y;
                glm::uvec3& triangles = scene.meshes[meshIndex].triangles[triangleIndex];
                glm::vec3 v1 = scene.meshes[meshIndex].vertices[triangles.x].position;
                glm::vec3 v2 = scene.meshes[meshIndex].vertices[triangles.y].position;
                glm::vec3 v3 = scene.meshes[meshIndex].vertices[triangles.z].position;
                float baryPoint = (v1[splitAxis] + v2[splitAxis] + v3[splitAxis]) / 3.0f;
                if (baryPoint <= centroid)
                    left.push_back(triangle);
                else
                    right.push_back(triangle);
            }
            AxisAlignedBox leftBox = getAABB(pScene, left);
            AxisAlignedBox rightBox = getAABB(pScene, right);
            // Calculate cost of this splitting
            float cost = calculateCost(leftBox, left.size(), rightBox, right.size());
            // Check if is smaller, if so, update centroidRes and splitAxisRes
            if (cost < minimalCost) {
                minimalCost = cost;
                centroidRes = centroid;
                splitAxisRes = splitAxis;
            }
        }
    }
    // SAH changes splitAxis and centroid
    return glm::vec2 { centroidRes, splitAxisRes };
}

// find the median of all triangles
float BoundingVolumeHierarchy::getSplitPoint(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes, int splitAxis)
{

    std::vector<float> splits;
    Scene scene = *pScene;

    for (auto& triangle : triangleIndexes) {
        int meshIndex = triangle.x;
        int triangleIndex = triangle.y;
        glm::uvec3 triangles = scene.meshes[meshIndex].triangles[triangleIndex];
        glm::vec3 v1 = scene.meshes[meshIndex].vertices[triangles.x].position;
        glm::vec3 v2 = scene.meshes[meshIndex].vertices[triangles.y].position;
        glm::vec3 v3 = scene.meshes[meshIndex].vertices[triangles.z].position;
        splits.push_back((v1[splitAxis] + v2[splitAxis] + v3[splitAxis]) / 3.0f);
    }

    // get median
    int n = splits.size() / 2;
    std::nth_element(splits.begin(), splits.begin() + n, splits.end());
    float xn = splits[n];

    if (splits.size() % 2 != 0)
        return xn;
    else {
        std::nth_element(splits.begin(), splits.begin() + n - 1, splits.end());
        return (xn + splits[n - 1]) / 2.0f;
    }
    // return xn;
}

// get AABB contains all triangles stored in the node
AxisAlignedBox BoundingVolumeHierarchy::getAABB(Scene* pScene, std::vector<glm::uvec2>& triangleIndexes)
{
    float xmin = 1e9;
    float ymin = 1e9;
    float zmin = 1e9;
    float xmax = -1e9;
    float ymax = -1e9;
    float zmax = -1e9;
    Scene scene = *pScene;

    for (auto& triangle : triangleIndexes) {
        int meshIndex = triangle.x;
        int triangleIndex = triangle.y;
        glm::uvec3 triangles = scene.meshes[meshIndex].triangles[triangleIndex];
        glm::vec3 v1 = scene.meshes[meshIndex].vertices[triangles.x].position;
        glm::vec3 v2 = scene.meshes[meshIndex].vertices[triangles.y].position;
        glm::vec3 v3 = scene.meshes[meshIndex].vertices[triangles.z].position;
        std::vector<glm::vec3> positions;
        positions.push_back(v1);
        positions.push_back(v2);
        positions.push_back(v3);
        for (auto& position : positions) {
            if (position.x < xmin)
                xmin = position.x;
            if (position.y < ymin)
                ymin = position.y;
            if (position.z < zmin)
                zmin = position.z;
            if (position.x > xmax)
                xmax = position.x;
            if (position.y > ymax)
                ymax = position.y;
            if (position.z > zmax)
                zmax = position.z;
        }
    }
    glm::vec3 lower { xmin, ymin, zmin };
    glm::vec3 upper { xmax, ymax, zmax };
    AxisAlignedBox aabb { lower, upper };
    return aabb;
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    // return 1;
    int level = std::numeric_limits<int>::min();
    for (Node node : nodeList) {
        if (node.level > level)
            level = node.level;
    }
    return level;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    // return 1;
    int count = 0;
    for (Node node : nodeList) {
        if (node.type == 1)
            count++;
    }
    return count;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level, const Features& features)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    // drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    glm::vec3 color {1.0f, 1.0f, 1.0f};
    if (features.extra.enableBvhSahBinning) {
        color = { 1.0f, 0.0f, 0.0f };
    }
    for (Node node : nodeList) {
        if (node.level == level) {
            drawAABB(node.aabb, DrawMode::Wireframe, color);
        }
    }
    // If BVH enabled draw AABB of normal BVH with color white
    if (features.extra.enableBvhSahBinning) {
        color = { 1.0f, 1.0f, 1.0f };
        for (Node node : nodeListBVH) {
            if (node.level == level) {
                drawAABB(node.aabb, DrawMode::Wireframe, color);
            }
        }
    }
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    // AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    // drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    // drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    /*
    std::vector<Node> leaves;
    for (Node& node : nodeList) {
        if (node.type == 1)
            leaves.push_back(node);
    }
    if (leafIdx < leaves.size())
        drawAABB(leaves[leafIdx].aabb, DrawMode::Wireframe);
    return;*/

    Scene scene = *m_pScene;
    std::vector<Node> leaves;
    for (Node& node : nodeList) {
        if (node.type == 1)
            leaves.push_back(node);
    }
    if (leafIdx < leaves.size()) {
        drawAABB(leaves[leafIdx].aabb, DrawMode::Wireframe);
        Node leaf = leaves[leafIdx];
        for (auto& triangle : leaf.child) {
            int meshIndex = triangle.x;
            int triangleIndex = triangle.y;
            glm::uvec3& triangles = scene.meshes[meshIndex].triangles[triangleIndex];
            Vertex v1 = scene.meshes[meshIndex].vertices[triangles.x];
            Vertex v2 = scene.meshes[meshIndex].vertices[triangles.y];
            Vertex v3 = scene.meshes[meshIndex].vertices[triangles.z];
            drawTriangle(v1, v2, v3);
        }
    }
    return;

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

// Check if origin is inside the origin
bool insideAABB(const AxisAlignedBox& box, const glm::vec3& point)
{
    if (point.x >= box.lower.x && point.x <= box.upper.x && point.y >= box.lower.y && point.y <= box.upper.y && point.z >= box.lower.z && point.z <= box.upper.z) {
        return true;
    }
    return false;
}

// Updates hitInfo based on features
void setHitInfo(Ray& ray, HitInfo& hitInfo, const Features& features, Vertex& a, Vertex& b, Vertex& c)
{
    glm::vec3 origin = ray.origin + ray.t * ray.direction;
    hitInfo.barycentricCoord = computeBarycentricCoord(a.position, b.position, c.position, origin);
    hitInfo.normal = glm::normalize(glm::cross(b.position - a.position, c.position - a.position));

    if (features.enableTextureMapping && hitInfo.material.kdTexture) {
        hitInfo.texCoord = interpolateTexCoord(a.texCoord, b.texCoord, c.texCoord, hitInfo.barycentricCoord);
    }

    if (features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(a.normal, b.normal, c.normal, hitInfo.barycentricCoord);
        Ray aNormal { a.position, a.normal };
        Ray bNormal { b.position, b.normal };
        Ray cNormal { c.position, c.normal };
        Ray normal { origin, hitInfo.normal };
        drawRay(aNormal, glm::vec3(1.0f, 0.0f, 0.0f));
        drawRay(bNormal, glm::vec3(1.0f, 0.0f, 0.0f));
        drawRay(cNormal, glm::vec3(1.0f, 0.0f, 0.0f));
        // Draw interpolated normal with white color
        drawRay(normal, glm::vec3(1.0f, 1.0f, 1.0f));
    }
}


bool BoundingVolumeHierarchy::traverseBVH(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    float t = ray.t;
    // If there is no hit and ray.origin is not inside aabb then stop traversing that aabb
    if (!intersectRayWithShape(headNode.aabb, ray) && !insideAABB(headNode.aabb, ray.origin))
        return false;
    ray.t = t;
    bool hit = false;

    std::priority_queue<std::pair<float, int>> q;
    q.push(std::make_pair(t, nodeList.size() - 1));

    // For debug
    Vertex a;
    Vertex b;
    Vertex c;
    std::vector<Mesh>& meshes = m_pScene->meshes; 

    while (!q.empty()) {
        auto tmp = q.top();
        q.pop();
        Node node = nodeList[tmp.second];
        if (features.enableBVHDebug) {
            drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3 { 0.0f, 1.0f, 0.0f });
        }
        // Node is a leaf
        if (node.type == 1) {
            // We are traversing along node's children
            for (glm::vec2 pairChild : node.child) {
                int meshIndex = pairChild.x;
                int triangleIndex = pairChild.y;

                Mesh& mesh = meshes[meshIndex];
                glm::vec3 triangle = mesh.triangles[triangleIndex];

                Vertex v0 = mesh.vertices[triangle.x];
                Vertex v1 = mesh.vertices[triangle.y];
                Vertex v2 = mesh.vertices[triangle.z];

                t = ray.t;
                // Check if there is intersection with this triangle
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo) && ray.t < t) {
                    hit = true;
                    hitInfo.material = mesh.material;
                    // Important for mipmapping
                    // For a mesh it throws an exeption
                    ray.meshIndex = meshIndex;
                    setHitInfo(ray, hitInfo, features, v0, v1, v2);
                    
                    // For debug
                    a = v0;
                    b = v1;
                    c = v2;
                } else {
                    // If no hit, update t value to the old one
                    ray.t = t;
                }
            }
            // If there is a hit and the next node's distance in the queue is bigger than the hit distance
            // Then we can just return
            if (hit && (q.size() == 0 || ray.t < std::abs(q.top().first))) {
                drawTriangle(a, b, c);
                return true;
            } 
        } else {
            // Node is an interiar
            // For each child check if it intersects and push to queue distance of this node
            t = ray.t;
            if (insideAABB(nodeList[node.child[0].x].aabb, ray.origin) || intersectRayWithShape(nodeList[node.child[0].x].aabb, ray)) {
                float distanceX = ray.t;
                q.push(std::make_pair(-distanceX, node.child[0].x));
                if (features.enableBVHDebug) {
                    drawAABB(nodeList[node.child[0].x].aabb, DrawMode::Wireframe, glm::vec3 { 1.0f, 0.0f, 0.0f });
                } else {
                    drawAABB(nodeList[node.child[0].x].aabb, DrawMode::Wireframe);
                }
            }
            ray.t = t;
            if (insideAABB(nodeList[node.child[0].y].aabb, ray.origin) || intersectRayWithShape(nodeList[node.child[0].y].aabb, ray)) {
                float distanceY = ray.t;
                q.push(std::make_pair(-distanceY, node.child[0].y));
                if (features.enableBVHDebug) {
                    drawAABB(nodeList[node.child[0].y].aabb, DrawMode::Wireframe, glm::vec3 { 1.0f, 0.0f, 0.0f });
                } else {
                    drawAABB(nodeList[node.child[0].y].aabb, DrawMode::Wireframe);
                }
            }
            ray.t = t;
        }
    }
    return hit;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        Vertex a;
        Vertex b;
        Vertex c;
        bool hit = false;
        int meshIndex = 0;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    // for mesh it throws an exeption
                    ray.meshIndex = meshIndex;
                    hit = true;
                    a = v0;
                    b = v1;
                    c = v2;
                }
            }
            meshIndex++;
        }

        setHitInfo(ray, hitInfo, features, a, b, c);

        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return traverseBVH(ray, hitInfo, features);
    }


}