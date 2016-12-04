#include "utility.h"
#include "define.h"

#pragma once

// Intersection Declarations
struct Intersection {
    // Intersection Public Methods
    Intersection() {
        triId = -1;
        uvt[0] = 0;
		uvt[1] = 0;
		uvt[2] = 0;
    }

    // Intersection Public Data
    int triId;
    float uvt[3];	
};



// Primitive Declarations
class Primitive {
public:
    // Primitive Interface
    Primitive(uint32_t triId) : primitiveId(triId) { }
    ~Primitive();
    BBox GetWorldBound();

	bool IntersectP(Ray &ray);
	bool Intersect(Ray &ray, Intersection *isect);
    
	void GetUVs(float uv[3][2]); 

	Vector GetKd(float u, float v);
	
	void GetKdKsNs(float u, float v, Vector &Kd, Vector &Ks, float &Ns);
    void GetKdKsNsEta(float u, float v, Vector &Kd, Vector &Ks, float &Ns, float &Eta);
	void GetKdKsNsEtaEmission(float u, float v, Vector &Kd, Vector &Ks, Vector &Emission, float &Ns, float &Eta);
    void GetMaterial(float u, float v, Vector &Kd, Vector &Ks, Vector &Kb, bool &hasBump, Vector &Emission, float &Ns, float &Eta, 
        MicroFacetType &MicroFacet, DistributionType &Distribution, float &Roughness, 
        Vector &Sigmna_a, Vector &Sigmna_s);

    char *GetMaterialName();
	
    // Primitive Public Data
    uint32_t primitiveId;
};


void UnitizeScene(float Scale);
bool BrouteForceTracing(Ray &ray, Intersection *insect);