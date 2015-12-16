#include <stdio.h>
#include <stdlib.h>
#include <glm/glm.hpp>
#include "utility.h"

#pragma once

class Camera{

public:
	Camera(glm::vec3 _CameraPos, glm::vec3 _CameraLookat, glm::vec3 _CameraUp,
		   float _Fstop, float _FocalLength, float _Aptertrue, float _SensorWidth, float _SensorHeight,
		   float _AspectRatio, int _ResWidth, int _ResHeight, float _FovV, float _FovH, 
		   bool _AllinFocus);
	~Camera();

	Ray CameraRay(int Pxl_x, int Pxl_y);
	
	glm::vec3 m_CameraPos;
	glm::vec3 m_CameraLookat;
	glm::vec3 m_CameraUp;

private:
	
	float m_Fstop; 
	float m_FocalLength; 
	float m_Aperture;
	float m_SensorWidth; 
	float m_SensorHeight;
	float m_HalfSensorWidth; 
	float m_HalfSensorHeight;
	float m_AspectRatio;
	int m_ResWidth;
	int m_ResHeight; 
	float m_FovV;
	float m_FovH;

	float m_PxlSizeW;
	float m_PxlSizeH;
	float m_PlaneScaleW;
	float m_PlaneScaleH;

	bool m_isPionhole;

	// Image Axis
	glm::vec3 m_U;
	glm::vec3 m_V;
	glm::vec3 m_W;

	float m_FocusDist;

	glm::vec3 m_ImgCenter;

	float m_tScale;

	// Camera Axis
	glm::vec3 m_CameraU;
	glm::vec3 m_CameraV;
	glm::vec3 m_CameraW;
	
};
