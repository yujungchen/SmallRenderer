#include "sampler.h"

double RandomNumber(){
	return (double) rand() / (RAND_MAX + 1.0 );
}

glm::vec2 UniformDiskSampling(){
	float r = RandomNumber();
	float theta = RandomNumber();
	glm::vec2 UnitSample = glm::vec2(0.0f, 0.0f);

	float rnd_r = sqrt(r);
	float rnd_angle = 2.0f * M_PI * theta;

	UnitSample.x = rnd_r * cos(rnd_angle);
	UnitSample.y = rnd_r * sin(rnd_angle);
	
	return UnitSample;
}

glm::vec3 CosHemiSampler(){
	float u1 = RandomNumber();
	float u2 = RandomNumber();
	float r = sqrt(u1);
	float theta = 2 * M_PI * u2;
 
	float x = r * cos(theta);
	float y = r * sin(theta);
 
	return glm::vec3(x, y, sqrt(glm::max(0.0f, 1 - u1)));
}