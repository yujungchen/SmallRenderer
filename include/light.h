#include <stdio.h>
#include <stdlib.h>
#include <glm/glm.hpp>
#include "glm.h"

#pragma once

class TestLight{

public:
	TestLight(glm::vec3 _lpos, glm::vec3 _lemission);
	~TestLight();
	glm::vec3 sampleL();
	glm::vec3 getlpos();
private:
	glm::vec3 m_lpos;
	glm::vec3 m_lemission;
};

typedef struct{
	glm::vec3 V[3];
	glm::vec3 N[3];
	glm::vec3 Kd;
	glm::vec3 Emission;
	float Area;
} LightTri;

class AreaLight{
public:
	AreaLight(GLMmodel *_model);
	~AreaLight();
	
private:
	int m_Num_lTri;
	LightTri *m_Tri_l;
};