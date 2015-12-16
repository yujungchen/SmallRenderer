#include <iostream>
#include <string>
#include <glm/glm.hpp>
#pragma once

using namespace std;

void ReadConfigure(char *SceneFile,
				   //Scene Configuration
				   string &m_Model,
				   float &m_SceneScale, 
				   //Camera Configuration
				   glm::vec3 &m_CameraPos, glm::vec3 &m_CameraUp, glm::vec3 &m_CameraLookat,
				   float &m_Fstop, bool &m_AllinFocus, float &m_FocalLength, float &m_Aperture, 
				   float &m_SensorWidth, float &m_SensorHeight, float &m_AspectRatio, 
				   int &m_ResWidth, int &m_ResHeight, float &m_FovV, float &m_FovH,
				   int &m_SPP, 
				   //Light Configuration
				   bool &m_TestLight, glm::vec3 &m_TestLightPos, glm::vec3 &m_TestLightEmission,
				   int &m_PathDepth, bool &m_NEE_Enable
				   );