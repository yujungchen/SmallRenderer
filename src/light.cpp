#include "light.h"

TestLight::TestLight(glm::vec3 _lpos, glm::vec3 _lemission){
	m_lpos = _lpos;
	m_lemission = _lemission;
}

TestLight::~TestLight(){}

glm::vec3 TestLight::sampleL(){
	return m_lemission;
}

glm::vec3 TestLight::getlpos(){
	return m_lpos;
}

