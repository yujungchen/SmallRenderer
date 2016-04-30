#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "bvh.h"
#include "light.h"
#include "dipole.h"
#include "camera.h"
#include "radiometry.h"

#pragma once


class PathIntegrator{
  
public:
	PathIntegrator(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
				   PointLight *_l, AreaLight *_al, Camera *_camera, bool _NEE_Enable, bool _UseAreaLight, 
				   Dipole *_Dipole);
	~PathIntegrator();

	glm::vec3 ComputeRadiance(int sample_x, int sample_y, int PathDepth);
	glm::vec3 NEE(glm::vec3 &Pos, glm::vec3 &PrevPos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, float Ns, float Eta,
				MicroFacetType &MicroFacetModel, DistributionType &Distribution, float &Roughness);

private:

	GLMmodel *m_model;
	BVHAccel *m_bvh;
	std::vector<Primitive> m_PrimList;

	PointLight *m_l;
	AreaLight *m_al;
	Camera *m_camera;
	Dipole *m_Dipole;

	bool m_NEE_Enable;
	bool m_UseAreaLight;

};