#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "bvh.h"

#pragma once


class Volume{

public:
	Volume();
	~Volume();

	void VolumeInit(char *_Name, glm::vec3 _Sigma_a, glm::vec3 _Sigma_s, float _N, GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList);
	glm::vec3 VolumeTracing(glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos);
	float MonoChrome(glm::vec3 refract, glm::vec3 pos, glm::vec3 N, float sigma_a, float sigma_s, float sigma_t);
	float SampleDist(float sigma_t);

private:
	GLMmodel *m_model;
	BVHAccel *m_bvh;
	std::vector<Primitive> m_PrimList;

	char m_MatName[20];
	glm::vec3 m_Sigma_a;
	glm::vec3 m_Sigma_s;
	glm::vec3 m_Sigma_t;

	float m_eta;
};