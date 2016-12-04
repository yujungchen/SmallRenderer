
#include <stdlib.h>
#include <vector>
#include "utility.h"
#include "core.h"
#include "glm.h"
#include <glm/glm.hpp>


#pragma once

extern GLMmodel* model;


struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() { }
    BVHPrimitiveInfo(int pn, BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
};

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};

// 2 4Bytes
struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};





// BVHAccel Declarations
class BVHAccel {
public:
    // BVHAccel Public Methods
    BVHAccel(uint32_t maxPrims = 1);
    BBox WorldBound();
    ~BVHAccel();

	BVHBuildNode *recursiveBuild(std::vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end, uint32_t *totalNodes, std::vector<Primitive > &orderedPrims);

	uint32_t flattenBVHTree(BVHBuildNode *node, uint32_t *offset);
	bool IntersectP(Ray &ray);
	bool Intersect(Ray &ray, Intersection *isect);
    //bool IntersectGeo(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, float &Ns, std::vector<Primitive> m_PrimList);
    void InterpolateGeo(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, float &Ns, float &Eta, std::vector<Primitive> &m_PrimList);
    void InterpolateGeoV2(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, glm::vec3 &Emission, float &Ns, float &Eta, std::vector<Primitive> &m_PrimList);

    void IsectGeometry(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, glm::vec3 &Kb, bool &hasBump, glm::vec3 &Emission, 
                        MicroFacetType &MicroFacet, DistributionType &Distribution, float &Roughness, float &Ns, float &Eta, std::vector<Primitive> &m_PrimList, 
                        glm::vec3 &Sigma_a, glm::vec3 &Sigma_s);

    char *GetMatName();

public:
    
    // BVHAccel Private Data
    std::vector<Primitive> primitives;
    uint32_t maxPrimsInNode;
    LinearBVHNode *nodes;

    char *m_MatName;
};

