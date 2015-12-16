#include <stdio.h>
#include <stdlib.h>
#include <glm/glm.hpp>

#pragma once

class TestLight{

public:
	TestLight(glm::vec3 _lpos, glm::vec3 _lemission);
	~TestLight();
	glm::vec3 sampleL();
	glm::vec3 getlpos();


private:

	glm::vec3 m_lpos;
	glm::vec3 m_lemission;

};