#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "bvh.h"
#include "light.h"
#include "camera.h"
#include "radiometry.h"

#pragma once


class PathIntegrator{
  
public:
	PathIntegrator(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
				   TestLight *_l, Camera *_camera, bool _NEE_Enable);
	~PathIntegrator();

	glm::vec3 ComputeRadiance(int sample_x, int sample_y, int PathDepth);
	glm::vec3 NEE(glm::vec3 lPos, glm::vec3 Pos, glm::vec3 PrevPos, glm::vec3 N, glm::vec3 Kd, glm::vec3 Ks, float Ns, glm::vec3 lEmission);

private:

	GLMmodel *m_model;
	BVHAccel *m_bvh;
	std::vector<Primitive> m_PrimList;

	TestLight *m_l;
	Camera *m_camera;

	bool m_NEE_Enable;

};