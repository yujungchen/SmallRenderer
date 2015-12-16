#include "camera.h"
#include "sampler.h"
#define EPSILON 0.0001f

Camera::Camera(glm::vec3 _CameraPos, glm::vec3 _CameraLookat, glm::vec3 _CameraUp,
		   float _Fstop, float _FocalLength, float _Aptertrue, float _SensorWidth, float _SensorHeight,
		   float _AspectRatio, int _ResWidth, int _ResHeight, float _FovV, float _FovH, 
		   bool _AllinFocus){

	m_CameraPos = _CameraPos;
	m_CameraLookat = _CameraLookat;
	m_CameraUp = _CameraUp;

	m_Fstop = _Fstop;
	m_FocalLength = _FocalLength;
	m_Aperture = _Aptertrue;
	m_SensorWidth = _SensorWidth;
	m_SensorHeight = _SensorHeight;
	m_AspectRatio = _AspectRatio;
	m_ResWidth = _ResWidth;
	m_ResHeight = _ResHeight;
	m_FovV = _FovV;
	m_FovH = _FovH;

	m_PxlSizeW = m_SensorWidth / (float)m_ResWidth;
	m_PxlSizeH = m_SensorHeight / (float)m_ResHeight;

	m_HalfSensorWidth = m_SensorWidth * 0.5f;
	m_HalfSensorHeight = m_SensorHeight * 0.5f;

	m_isPionhole = false;
	// 0.0001mm
	if(m_Aperture < 0.0000001)
		m_isPionhole = true;

	if(_AllinFocus)
		m_isPionhole = true;

	m_U = glm::vec3(0.0, 0.0, 0.0);
	m_V = glm::vec3(0.0, 0.0, 0.0);
	m_W = glm::vec3(0.0, 0.0, 0.0);

	glm::vec3 m_W = m_CameraLookat - m_CameraPos;
	m_FocusDist = glm::length(m_W);
	LocalBasis(glm::normalize(m_W), &m_U, &m_V);

	m_W = glm::normalize(m_W) * m_FocalLength;
	m_ImgCenter = m_CameraPos + m_W;

	m_tScale = m_FocusDist / m_FocalLength;

	m_CameraW = glm::normalize(m_CameraLookat - m_CameraPos);
	LocalBasis(m_CameraW, &m_CameraU, &m_CameraV);
	m_CameraU = glm::normalize(m_CameraU);
	m_CameraV = glm::normalize(m_CameraV);

}

Camera::~Camera(){

}

Ray Camera::CameraRay(int Pxl_x, int Pxl_y){
	Vector Direction(0.0, 0.0f, 0.0f);
	Point Eye = Point(m_CameraPos.x, m_CameraPos.y, m_CameraPos.z);

	glm::vec3 Dir(0.0, 0.0, 0.0);
	Dir.x = -m_HalfSensorWidth + (RandomNumber() + Pxl_x) * m_PxlSizeW;
	Dir.y = -m_HalfSensorHeight + (RandomNumber() + Pxl_y) * m_PxlSizeH;

	glm::vec3 ws_PlaneSamplePos = m_ImgCenter + Dir.x * m_U + Dir.y * m_V;
	glm::vec3 RayDir = (ws_PlaneSamplePos - m_CameraPos);
	RayDir = glm::normalize(RayDir);

	if(!m_isPionhole){
		glm::vec3 FocusP = m_CameraPos + m_tScale * RayDir;
		
		// Lens Sample
		glm::vec2 UnitSample = UniformDiskSampling();
		float lens_u = UnitSample.x * m_Aperture * 0.5f;
		float lens_v = UnitSample.y * m_Aperture * 0.5f;
		glm::vec3 lensOffset(0.0f, 0.0f, 0.0f);

		lensOffset = m_CameraU * lens_u + m_CameraV * lens_v;
		glm::vec3 lensSample = m_CameraPos + lensOffset;
		RayDir = (FocusP - lensSample);
		RayDir = glm::normalize(RayDir);

		Eye = Point(lensSample.x, lensSample.y, lensSample.z);
	}

	Direction = Vector(RayDir.x, RayDir.y, RayDir.z);

	return Ray(Eye, Direction, EPSILON);
}