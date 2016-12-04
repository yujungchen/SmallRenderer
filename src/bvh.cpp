#include "bvh.h"




struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) {
        return a.centroid[dim] < b.centroid[dim];
    }
};



struct CompareToBucket {
    CompareToBucket(int split, int num, int d, BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(BVHPrimitiveInfo &p);

    int splitBucket, nBuckets, dim;
    BBox &centroidBounds;
};


bool CompareToBucket::operator()(BVHPrimitiveInfo &p) {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    return b <= splitBucket;
}


static inline bool IntersectBox(BBox &bounds, Ray &ray, Vector &invDir, uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}



BVHBuildNode *BVHAccel::recursiveBuild(std::vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end, uint32_t *totalNodes, std::vector<Primitive > &orderedPrims) {
	(*totalNodes)++;
	BVHBuildNode *node = new BVHBuildNode();
	  // Compute bounds of all primitives in BVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
	uint32_t nPrimitives = end - start;
	if (nPrimitives == 1) {
		// Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    } 
	else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
		if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // If nPrimitives is no greater than maxPrimsInNode,
            // then all the nodes can be stored in a compact bvh node.
            if (nPrimitives <= maxPrimsInNode) {
                // Create leaf _BVHBuildNode_
                uint32_t firstPrimOffset = orderedPrims.size();
                for (uint32_t i = start; i < end; ++i) {
                    uint32_t primNum = buildData[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                return node;
            }
            else {
                // else if nPrimitives is greater than maxPrimsInNode, we
                // need to split it further to guarantee each node contains
                // no more than maxPrimsInNode primitives.
                node->InitInterior(dim,
                                   recursiveBuild(buildData, start, mid,
                                                  totalNodes, orderedPrims),
                                   recursiveBuild(buildData, mid, end,
                                                  totalNodes, orderedPrims));
                return node;
            }
        }

		// SAH split: Partition primitives using approximate SAH
		if (nPrimitives <= 4) {
			// Partition primitives into equally-sized subsets
			mid = (start + end) / 2;
			std::nth_element(&buildData[start], &buildData[mid],
							 &buildData[end-1]+1, ComparePoints(dim));
		}
		else {
			// Allocate _BucketInfo_ for SAH partition buckets
			const int nBuckets = 12;
			struct BucketInfo {
				BucketInfo() { count = 0; }
				int count;
				BBox bounds;
			};
			BucketInfo buckets[nBuckets];

			// Initialize _BucketInfo_ for SAH partition buckets
			for (uint32_t i = start; i < end; ++i) {
				int b = nBuckets * ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) / (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
				if (b == nBuckets) b = nBuckets-1;
					buckets[b].count++;
					buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
			}

			// Compute costs for splitting after each bucket
			float cost[nBuckets-1];
			for (int i = 0; i < nBuckets-1; ++i) {
				BBox b0, b1;
				int count0 = 0, count1 = 0;
				for (int j = 0; j <= i; ++j) {
					b0 = Union(b0, buckets[j].bounds);
					count0 += buckets[j].count;
				}
				for (int j = i+1; j < nBuckets; ++j) {
					b1 = Union(b1, buckets[j].bounds);
					count1 += buckets[j].count;
				}
				cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) / bbox.SurfaceArea();
			}

			// Find bucket to split at that minimizes SAH metric
			float minCost = cost[0];
			uint32_t minCostSplit = 0;
			for (int i = 1; i < nBuckets-1; ++i) {
				if (cost[i] < minCost) {
					minCost = cost[i];
					minCostSplit = i;
				}
			}

			// Either create leaf or split primitives at selected SAH bucket
			if (nPrimitives > maxPrimsInNode || minCost < nPrimitives) {
				BVHPrimitiveInfo *pmid = std::partition(&buildData[start], &buildData[end-1]+1, CompareToBucket(minCostSplit, (uint32_t)nBuckets, dim, centroidBounds));
				mid = pmid - &buildData[0];
			}
           
			else {
				// Create leaf _BVHBuildNode_
				uint32_t firstPrimOffset = orderedPrims.size();
				for (uint32_t i = start; i < end; ++i) {
					uint32_t primNum = buildData[i].primitiveNumber;
					orderedPrims.push_back(primitives[primNum]);
				}
				node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
				return node;
			}
		}
		 node->InitInterior(dim,
                           recursiveBuild(buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
	return node;
}

// BVHAccel Method Definitions
BVHAccel::BVHAccel(uint32_t mp) {
	uint32_t numtriangles = model->numtriangles;
	if (numtriangles == 0) {
        nodes = NULL;
        return;
    }


	// Initialize _buildData_ array for primitives
	std::vector<BVHPrimitiveInfo> buildData;
    buildData.reserve(numtriangles);
	for (uint32_t i = 0; i < numtriangles; i++) {
		Primitive *pri = new Primitive(i);
		primitives.push_back(*pri);
        BBox bbox = pri->GetWorldBound();
        buildData.push_back(BVHPrimitiveInfo(i, bbox));
    }

	//for (uint32_t i = 0; i < buildData.size(); i++) {
	//	printf("%f %f %f\n", buildData[i].centroid.x, buildData[i].centroid.y, buildData[i].centroid.z);
   // }
	//printf("Number of primitives: %d\n", primitives.size());

    // Recursively build BVH tree for primitives
    uint32_t totalNodes = 0;
    std::vector<Primitive> orderedPrims;
    orderedPrims.reserve(primitives.size());

	
    BVHBuildNode *root = recursiveBuild(buildData, 0, primitives.size(), &totalNodes, orderedPrims);
	primitives.swap(orderedPrims);
	nodes = (LinearBVHNode *)malloc(sizeof(LinearBVHNode) * totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBVHNode;
	uint32_t offset = 0;
    flattenBVHTree(root, &offset);

}


uint32_t BVHAccel::flattenBVHTree(BVHBuildNode *node, uint32_t *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1],
                                                       offset);
    }
    return myOffset;
}



bool BVHAccel::IntersectP(Ray &ray) {
	if (!nodes) 
		return false;
	Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	uint32_t todo[64];
	uint32_t todoOffset = 0, nodeNum = 0;
	while (true) {
		LinearBVHNode *node = &nodes[nodeNum];
		if (::IntersectBox(node->bounds, ray, invDir, dirIsNeg)) {
			// Process BVH node _node_ for traversal
			if (node->nPrimitives > 0) {
				for (uint32_t i = 0; i < node->nPrimitives; ++i) {
					if (primitives[node->primitivesOffset+i].IntersectP(ray)) {
						return true;
					}
				}
				if (todoOffset == 0) 
					break;
				nodeNum = todo[--todoOffset];
			} else {
				if (dirIsNeg[node->axis]) {
					// second child first
					todo[todoOffset++] = nodeNum + 1;
					nodeNum = node->secondChildOffset;
				} else {
					todo[todoOffset++] = node->secondChildOffset;
					nodeNum = nodeNum + 1;
				}
			}
		} else {
			if (todoOffset == 0) 
				break;
			nodeNum = todo[--todoOffset];
		}
	}
    
	return false;
}

bool BVHAccel::Intersect(Ray &ray, Intersection *isect) {
	if (!nodes) 
		return false;
    
	bool hit = false;
	Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
	uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	// Follow ray through BVH nodes to find primitive intersections
	uint32_t todoOffset = 0, nodeNum = 0;
	uint32_t todo[64];

	while (true) {
		LinearBVHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectBox(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {// leaf node
                // Intersect ray with primitives in leaf BVH node
                for (uint32_t i = 0; i < node->nPrimitives; ++i) {
                    if (primitives[node->primitivesOffset+i].Intersect(ray, isect)) {
						hit = true;
					} 
				}
                if (todoOffset == 0) 
					break;
				nodeNum = todo[--todoOffset];
			} else {// internal node
				// Put far BVH node on _todo_ stack, advance to near node
				if (dirIsNeg[node->axis]) {
					todo[todoOffset++] = nodeNum + 1;
					nodeNum = node->secondChildOffset;
				} else {
					todo[todoOffset++] = node->secondChildOffset;
					nodeNum = nodeNum + 1;
				}
			}
		} else {
			if (todoOffset == 0) 
				break;
		nodeNum = todo[--todoOffset];
		}
	}
    return hit;
}


void BVHAccel::InterpolateGeo(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, float &Ns, float &Eta, std::vector<Primitive> &m_PrimList){
	// Store Geo
	// Position
	Pos = glm::vec3(0.0f, 0.0f, 0.0f);
	float t = isect->uvt[2];
	Point hitP = ray.o + t * ray.d;
	Pos = glm::vec3(hitP.x, hitP.y, hitP.z); 

	// Normal
	N = glm::vec3(0.0f, 0.0f, 0.0f);
	Normal _N;
	Normal n0(model->normals[model->triangles[isect->triId].nindices[0] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[0] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[0] * 3 + 2]);
	Normal n1(model->normals[model->triangles[isect->triId].nindices[1] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[1] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[1] * 3 + 2]);
	Normal n2(model->normals[model->triangles[isect->triId].nindices[2] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[2] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[2] * 3 + 2]);
	_N = BaryInterpolationN(n0, n1, n2, isect->uvt[0], isect->uvt[1]);
	N = glm::vec3(_N.x, _N.y, _N.z);


	// Material and Albedo
	Kd = glm::vec3(0.0f, 0.0f, 0.0f);
	Ks = glm::vec3(0.0f, 0.0f, 0.0f);
	Ns = 0.0f;
	
	Vector _Kd, _Ks;
	float _Ns;	
	float _Eta;
	m_PrimList[isect->triId].GetKdKsNsEta(isect->uvt[0], isect->uvt[1], _Kd, _Ks, _Ns, _Eta);

	Kd = glm::vec3(_Kd.x, _Kd.y, _Kd.z);
	Ks = glm::vec3(_Ks.x, _Ks.y, _Ks.z);
	Ns = _Ns;
	Eta = _Eta;
}

void BVHAccel::InterpolateGeoV2(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, glm::vec3 &Emission, float &Ns, float &Eta, std::vector<Primitive> &m_PrimList){
	// Store Geo
	Vector _Kd, _Ks, _Emission;
	float _Ns;	
	float _Eta;
	m_PrimList[isect->triId].GetKdKsNsEtaEmission(isect->uvt[0], isect->uvt[1], _Kd, _Ks, _Emission, _Ns, _Eta);

	// Position
	Pos = glm::vec3(0.0f, 0.0f, 0.0f);
	float t = isect->uvt[2];
	Point hitP = ray.o + t * ray.d;
	Pos = glm::vec3(hitP.x, hitP.y, hitP.z); 

	// Normal
	N = glm::vec3(0.0f, 0.0f, 0.0f);
	Normal _N;
	Normal n0(model->normals[model->triangles[isect->triId].nindices[0] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[0] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[0] * 3 + 2]);
	Normal n1(model->normals[model->triangles[isect->triId].nindices[1] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[1] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[1] * 3 + 2]);
	Normal n2(model->normals[model->triangles[isect->triId].nindices[2] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[2] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[2] * 3 + 2]);
	_N = BaryInterpolationN(n0, n1, n2, isect->uvt[0], isect->uvt[1]);
	N = glm::vec3(_N.x, _N.y, _N.z);


	// Material and Albedo
	Kd = glm::vec3(0.0f, 0.0f, 0.0f);
	Ks = glm::vec3(0.0f, 0.0f, 0.0f);
	Ns = 0.0f;
	
	Kd = glm::vec3(_Kd.x, _Kd.y, _Kd.z);
	Ks = glm::vec3(_Ks.x, _Ks.y, _Ks.z);
	Emission = glm::vec3(_Emission.x, _Emission.y, _Emission.z);
	Ns = _Ns;

	Eta = _Eta;
}


char *BVHAccel::GetMatName(){
	return m_MatName;
}

void BVHAccel::IsectGeometry(Ray &ray, Intersection *isect, glm::vec3 &Pos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, glm::vec3 &Kb, bool &hasBump, glm::vec3 &Emission, 
	MicroFacetType &MicroFacet, DistributionType &Distribution, float &Roughness, float &Ns, float &Eta, std::vector<Primitive> &m_PrimList, 
	glm::vec3 &Sigma_a, glm::vec3 &Sigma_s) {
	
	// Store Geo
	Vector _Kd, _Ks, _Kb, _Emission;
	float _Ns;	
	float _Eta;

	Vector _Sigma_a;
	Vector _Sigma_s;

	m_PrimList[isect->triId].GetMaterial(isect->uvt[0], isect->uvt[1], _Kd, _Ks, _Kb, hasBump, _Emission, _Ns, _Eta, MicroFacet, Distribution, Roughness, 
		_Sigma_a, _Sigma_s);

	Sigma_a = glm::vec3(_Sigma_a.x, _Sigma_a.y, _Sigma_a.z);
	Sigma_s = glm::vec3(_Sigma_s.x, _Sigma_s.y, _Sigma_s.z);

	//printf("%s\n", m_PrimList[isect->triId].GetMaterialName());
	m_MatName = m_PrimList[isect->triId].GetMaterialName();

	// Position
	Pos = glm::vec3(0.0f, 0.0f, 0.0f);
	float t = isect->uvt[2];
	Point hitP = ray.o + t * ray.d;
	Pos = glm::vec3(hitP.x, hitP.y, hitP.z); 

	// Normal
	N = glm::vec3(0.0f, 0.0f, 0.0f);
	Normal _N;
	Normal n0(model->normals[model->triangles[isect->triId].nindices[0] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[0] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[0] * 3 + 2]);
	Normal n1(model->normals[model->triangles[isect->triId].nindices[1] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[1] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[1] * 3 + 2]);
	Normal n2(model->normals[model->triangles[isect->triId].nindices[2] * 3 + 0], model->normals[model->triangles[isect->triId].nindices[2] * 3 + 1], model->normals[model->triangles[isect->triId].nindices[2] * 3 + 2]);
	_N = BaryInterpolationN(n0, n1, n2, isect->uvt[0], isect->uvt[1]);
	N = glm::vec3(_N.x, _N.y, _N.z);


	// Material and Albedo
	Kd = glm::vec3(0.0f, 0.0f, 0.0f);
	Ks = glm::vec3(0.0f, 0.0f, 0.0f);
	Kb = glm::vec3(0.0f, 0.0f, 0.0f);
	Ns = 0.0f;
	
	Kd = glm::vec3(_Kd.x, _Kd.y, _Kd.z);
	Ks = glm::vec3(_Ks.x, _Ks.y, _Ks.z);
	Kb = glm::vec3(_Kb.x, _Kb.y, _Kb.z);

	if(hasBump){
		glm::vec3 T = glm::vec3(0.0);
		glm::vec3 B = glm::vec3(0.0);

		LocalBasis(N, &T, &B);
		N = T * Kb.x + B * Kb.y + N * Kb.z;
		
	}


	Emission = glm::vec3(_Emission.x, _Emission.y, _Emission.z);
	Ns = _Ns;

	Eta = _Eta;
}