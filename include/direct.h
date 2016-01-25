#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "bvh.h"
#include "light.h"
#include "camera.h"
#include "radiometry.h"

#pragma once

class DirectIllumination{

public:
	DirectIllumination(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
					   PointLight *_l, AreaLight *_al, Camera *_camera, 
					   int _Width, int _Height, int _PathSample);
	~DirectIllumination();
	void Render(glm::vec3 *m_Img, int SampleNumber);


private:

	GLMmodel *m_model;
	BVHAccel *m_bvh;
	std::vector<Primitive> m_PrimList;

	PointLight *m_l;
	AreaLight *m_al;
	Camera *m_camera;

	int m_Width;
	int m_Height;
	int m_PathSample;

	bool m_useArealight;

};