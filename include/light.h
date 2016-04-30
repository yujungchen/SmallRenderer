#include <stdio.h>
#include <stdlib.h>
#include <glm/glm.hpp>
#include "glm.h"

#pragma once

class PointLight{

public:
	PointLight(glm::vec3 _lpos, glm::vec3 _lemission);
	~PointLight();
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
	float pdf;
} LightTri;

class AreaLight{
public:
	AreaLight(GLMmodel *_model);
	~AreaLight();
	float ComputeArea(glm::vec3 A, glm::vec3 B, glm::vec3 C);
	float getPdf(int TriIdx);
	glm::vec3 sampleL(glm::vec3 &Pos, glm::vec3 &N);

private:
	int m_Num_lTri;
	LightTri *m_Tri_l;
	

};