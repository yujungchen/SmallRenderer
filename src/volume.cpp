#include <math.h>
#include <string.h>
#include "utility.h"
#include "core.h"
#include "sampler.h"
#include "volume.h"

#define EPSILON 0.00001f

Volume::Volume(){

}

Volume::~Volume(){
	
}

void Volume::VolumeInit(char *_Name, glm::vec3 _Sigma_a, glm::vec3 _Sigma_s, float _N, GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList){
	m_Sigma_a = _Sigma_a;
	m_Sigma_s = _Sigma_s;
	m_Sigma_t = m_Sigma_s + m_Sigma_a;
	m_eta = 1.0f / _N;
	//printf("%f\n", m_eta);
	m_model = _model;
	m_bvh = _bvh;
	m_PrimList = _PrimList;
}

float Volume::SampleDist(float sigma_t){
	float dist = 0.0f;

	dist = (-1.0f) * log((float)RandomNumber()) / sigma_t;
	return dist;

}

float Volume::MonoChrome(glm::vec3 refract, glm::vec3 pos, glm::vec3 N, float sigma_a, float sigma_s, float sigma_t) {
	float Channel = 0.0f;

	float t = 0.0f;
	float Diffusion = 0.0f;
	float Source = 0.0f;
	glm::vec3 PhaseDir = glm::vec3(0.0f);
	glm::vec3 ScatterDir = glm::vec3(refract.x, refract.y, refract.z);
	Diffusion = sigma_s / sigma_t;


	t = SampleDist(sigma_t);	//mm-1
	t = t * 0.0001f;			//m
	Ray v_ray = Ray(Point(pos.x, pos.y, pos.z), Vector(ScatterDir.x, ScatterDir.y, ScatterDir.z), EPSILON);
	Point NextDest = v_ray.o + t * v_ray.d;
	PhaseDir = SphereSampler();
	
	printf("%f %f %f %f\n", t, PhaseDir.x, PhaseDir.y, PhaseDir.z);

	return Channel;
}

glm::vec3 Volume::VolumeTracing(glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos) {
	glm::vec3 irad = glm::vec3(0.0f);

	glm::vec3 incident = glm::normalize(prev_pos - pos);
	glm::vec3 refract = glm::refract(incident, N, m_eta);
	float r = MonoChrome(refract, pos, N, m_Sigma_a.x, m_Sigma_s.x, m_Sigma_t.x);

	return irad;
}