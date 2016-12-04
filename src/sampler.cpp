#include "stdio.h"
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

glm::vec3 SphereSampler(){

	float theta = 2 * M_PI * RandomNumber();
	float phi = acos(1 - 2 * RandomNumber());
	float x = sin(phi) * cos(theta);
	float y = sin(phi) * sin(theta);
	float z = cos(phi);

	return glm::vec3(x, y, z);	
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

glm::vec3 TriBaryUVSampler(){
	float u1 = RandomNumber();
	float u2 = RandomNumber();
	float u1sqrt = sqrt(u1);

	float U = 1.0f - u1sqrt;
	float V = (1.0f - u2) * u1sqrt;
	float T = u2 * u1sqrt;

	if(U + V + T - 1.0f > 0.01f){
		printf("Bad bary sample %f %f %f = %.10f\n", U, V, T, U + V + T);
	}

	return glm::vec3(U, V, T);
}

glm::vec3 BlinnPhongSampler(float Ns){
	float u1 = RandomNumber();
	float u2 = RandomNumber();
	// PBRT and Mitsuba use (Ns + 1.0f)
	float CosTheta = powf(u1, 1.0f / (Ns + 1.0f) );
	float Phi = u2 * 2.0f * (float)M_PI;
	float SinTheta = sqrt(1.0f - CosTheta * CosTheta);
  
	float x = cos(Phi) * SinTheta;
	float y = sin(Phi) * SinTheta;
	float z = CosTheta;
	return glm::vec3(x, y, z);
}