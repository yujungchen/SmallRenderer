

#pragma once
char *SceneFile,
				   //Scene Configuration
				   float &m_SceneScale, 
				   //Camera Configuration
				   glm::vec3 &m_CameraPos, glm::vec3 &m_CameraUp, glm::vec3 &m_CameraLookat,
				   float &m_Fstop, float &m_FocalLength, float &m_Aperture, 
				   float &m_SensorWidth, float &m_SensorHeight, float &m_AspectRatio, 
				   int &m_ResWidth, int &m_ResHeight, float &m_FovV, float &m_FovH,
				   //Light Configuration
				   bool &m_TestLight, glm::vec3 &m_TestLightPos, glm::vec3 &m_TestLightEmission