#include "microfacet.h"
#include "utility.h"


glm::vec3 SampleMicroNormal(glm::vec3 &N, float Ns, DistributionType &Distribution) {
	glm::vec3 MicroNormal = glm::vec3(0.0f);

	// Compute Normal using BlinnPhong Distribution Function
	glm::vec3 U = glm::vec3(0.0f);
	glm::vec3 V = glm::vec3(0.0f);
	LocalBasis(N, &U, &V);
	glm::vec3 BlinnPhongSample = BlinnPhongSampler(Ns);
	MicroNormal = U * BlinnPhongSample.x + V * BlinnPhongSample.y + N * BlinnPhongSample.z;
	MicroNormal = glm::normalize(MicroNormal);

	return MicroNormal;
}


glm::vec3 ComputeTorranceMicroFacetBRDF(glm::vec3 &Pos0, glm::vec3 &Pos1, glm::vec3 &Pos2, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, 
	float Ns, DistributionType &Distribution, float &Roughness,
	glm::vec3 &MicroNormal){
	
	glm::vec3 BRDF = glm::vec3(0.0f);

	// Compute Torrance Sparrow Microfacet Model
	// Use Blinn-Phong Distribution Function
	// BRDF = D * F * G / (4 * Dot(N, L) * Dot(N, V));
	MicroNormal = SampleMicroNormal(N, Ns, Distribution);
	
	// Compute Outgoing Vector
	glm::vec3 negwi = glm::normalize(Pos1 - Pos0);
	glm::vec3 wo = glm::normalize(Pos2 - Pos1);
	glm::vec3 wi = negwi * (-1.0f);

	// Compute HalfVector
	glm::vec3 w_half = glm::normalize(wo + wi);

	float NdotL = glm::dot(N, wo); // N dot wo

	if(NdotL == 0.0f)
		return BRDF;
	
	float NdotH = fmax(glm::dot(N, w_half), 0.0f);
	float NdotV = fmax(glm::dot(N, wi), 0.0f);
	float VdotH = fmax(glm::dot(wi, w_half), 0.0f);

	if(NdotV == 0.0f)
		return BRDF;
	
	// Compute D
	float D = (Ns + 2.0f) * 0.5f * INV_PI * pow(glm::dot(MicroNormal, w_half), Ns);
	
	// Compute G
	float G1 = 2.0f * NdotH * NdotV / VdotH;
	float G2 = 2.0f * NdotH * NdotL / VdotH;
	float G = fmin(1.0, fmin(G1, G2));

	// Compute F
	float OneMinusVdotH = (1.0f - VdotH);
	float F = OneMinusVdotH * OneMinusVdotH * OneMinusVdotH * OneMinusVdotH * OneMinusVdotH;
	float F0 = 0.8;
	F = F * (1.0 - F0);
	F = F + F0;


	BRDF = Kd * INV_PI + Ks * F * D * G / (4.0f * NdotV * NdotL);

	return BRDF;
}
