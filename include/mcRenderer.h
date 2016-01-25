#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "bvh.h"
#include "light.h"
#include "camera.h"
#include "path.h"
#include "direct.h"

#pragma once

class MCRenderer{
  
public:
	MCRenderer(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
			   PointLight *_l, AreaLight *_al, Camera *_camera, 
			   int _Width, int _Height, float _AspectRatio, int _PathSample, float _FocusDist,
			   int _PathDepth, bool _NEE_Enable);
	~MCRenderer();

	void Render();

private:

	GLMmodel *m_model;
	BVHAccel *m_bvh;
	std::vector<Primitive> m_PrimList;

	PointLight *m_l;
	AreaLight *m_al;
	Camera *m_camera;
	DirectIllumination *m_Direct;

	int m_Width;
	int m_Height;
	int m_DirectSampleNum;
	float m_AspectRatio;
	unsigned char *m_ColorImg;
	glm::vec3 *m_Img;
	glm::vec3 *m_DirectImg;
	int m_PathSample;
	float m_FocusDist;
	int m_PathDepth;
	bool m_NEE_Enable;

};