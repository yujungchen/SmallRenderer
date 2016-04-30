#include "light.h"
#include "sampler.h"
#include "utility.h"

PointLight::PointLight(glm::vec3 _lpos, glm::vec3 _lemission){
	m_lpos = _lpos;
	m_lemission = _lemission;
}

PointLight::~PointLight(){}

glm::vec3 PointLight::sampleL(){
	return m_lemission;
}

glm::vec3 PointLight::getlpos(){
	return m_lpos;
}


float AreaLight::ComputeArea(glm::vec3 A, glm::vec3 B, glm::vec3 C) {
	float area = 0.0f;
	
	float d_a = glm::length(A - B);	
	float d_b = glm::length(B - C);
	float d_c = glm::length(C - A);

	float s = (d_a + d_b + d_c) * 0.5f;
	area = sqrt(s * (s - d_a) * (s - d_b) * (s - d_c));

	return area;
}

float AreaLight::getPdf(int TriIdx) {
	return m_Tri_l[TriIdx].pdf;
}

/* Define Area Light Primitive*/
AreaLight::AreaLight(GLMmodel *model){
	m_Num_lTri = model->numLightTri;
	m_Tri_l = (LightTri*)malloc(sizeof(LightTri) * m_Num_lTri);

	for(int Idx = 0 ; Idx < model->numLightTri ; Idx++){
		m_Tri_l[Idx].Kd.x = model->materials[model->m_areaLight[Idx].matId].diffuse[0];
		m_Tri_l[Idx].Kd.y = model->materials[model->m_areaLight[Idx].matId].diffuse[1];
		m_Tri_l[Idx].Kd.z = model->materials[model->m_areaLight[Idx].matId].diffuse[2];
		m_Tri_l[Idx].Emission.x = model->materials[model->m_areaLight[Idx].matId].emissive[0];
		m_Tri_l[Idx].Emission.y = model->materials[model->m_areaLight[Idx].matId].emissive[1];
		m_Tri_l[Idx].Emission.z = model->materials[model->m_areaLight[Idx].matId].emissive[2];
		//printf("%f ", m_Tri_l[Idx].Kd.x);	printf("%f ", m_Tri_l[Idx].Kd.y);	printf("%f\n", m_Tri_l[Idx].Kd.z);	
		//printf("%f ", m_Tri_l[Idx].Emission.x);	printf("%f ", m_Tri_l[Idx].Emission.y);	printf("%f\n", m_Tri_l[Idx].Emission.z);
		
		//printf("Area light %d\n", Idx);
		for(int vIdx = 0 ; vIdx < 3 ; vIdx++){
			m_Tri_l[Idx].V[vIdx].x = model->vertices[model->m_areaLight[Idx].vindices[vIdx] * 3 + 0];
			m_Tri_l[Idx].V[vIdx].y = model->vertices[model->m_areaLight[Idx].vindices[vIdx] * 3 + 1];
			m_Tri_l[Idx].V[vIdx].z = model->vertices[model->m_areaLight[Idx].vindices[vIdx] * 3 + 2];

			m_Tri_l[Idx].N[vIdx].x = model->normals[model->m_areaLight[Idx].nindices[vIdx] * 3 + 0];
			m_Tri_l[Idx].N[vIdx].y = model->normals[model->m_areaLight[Idx].nindices[vIdx] * 3 + 1];
			m_Tri_l[Idx].N[vIdx].z = model->normals[model->m_areaLight[Idx].nindices[vIdx] * 3 + 2];
			//printf("[V%d] %f ", vIdx, m_Tri_l[Idx].V[vIdx].x);	printf("%f ", m_Tri_l[Idx].V[vIdx].y);	printf("%f\n", m_Tri_l[Idx].V[vIdx].z);
			//printf("[N%d] %f ", vIdx, m_Tri_l[Idx].N[vIdx].x);	printf("%f ", m_Tri_l[Idx].N[vIdx].y);	printf("%f\n", m_Tri_l[Idx].N[vIdx].z);
		}
		//printf("\n");

		m_Tri_l[Idx].Area = ComputeArea(m_Tri_l[Idx].V[0], m_Tri_l[Idx].V[1], m_Tri_l[Idx].V[2]);
		m_Tri_l[Idx].pdf = 1.0f / m_Tri_l[Idx].Area;
		m_Tri_l[Idx].pdf = m_Tri_l[Idx].pdf/ (float)model->numLightTri;
		printf("Light Tri %d Area : %f Pdf : %f\n", Idx, m_Tri_l[Idx].Area, m_Tri_l[Idx].pdf);
	}


}

AreaLight::~AreaLight(){}

glm::vec3 AreaLight::sampleL(glm::vec3 &Pos, glm::vec3 &N){
	glm::vec3 lemission = glm::vec3(0.0, 0.0, 0.0);
	int TriNum = rand() % m_Num_lTri;

	// Sample the triangle
	if(TriNum >= m_Num_lTri){
		printf("Sample the wrong index.\n");
		exit(0);
	}

	LightTri SampleTri = m_Tri_l[TriNum];

	// Barycentric sample coefficient
	glm::vec3 BaryCoef = TriBaryUVSampler();
	Pos = BarycentricInterpolation(SampleTri.V[0], SampleTri.V[1], SampleTri.V[2], BaryCoef);
	// Trick : Assuming the Normals are the same of a primitive.
	//N = BarycentricInterpolation(SampleTri.N[0], SampleTri.N[1], SampleTri.N[2], BaryCoef);
	N = SampleTri.N[0];

/*
	glm::vec3 LocalP = glm::vec3(0.0f);
	glm::vec3 U = glm::vec3(0.0f);
	glm::vec3 V = glm::vec3(0.0f);
	LocalP = CosHemiSampler();
	LocalBasis(TempN, &U, &V);
	N = V * LocalP.x + U * LocalP.y + TempN * LocalP.z;
*/
	lemission = SampleTri.Kd * SampleTri.Emission;
	return lemission;
}