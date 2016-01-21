#include "light.h"

TestLight::TestLight(glm::vec3 _lpos, glm::vec3 _lemission){
	m_lpos = _lpos;
	m_lemission = _lemission;
}

TestLight::~TestLight(){}

glm::vec3 TestLight::sampleL(){
	return m_lemission;
}

glm::vec3 TestLight::getlpos(){
	return m_lpos;
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
		
		printf("Area light %d\n", Idx);
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
	}

}

AreaLight::~AreaLight(){}
