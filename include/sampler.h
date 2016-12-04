#include <glm/glm.hpp>
#include <math.h>

#pragma once

double RandomNumber();
glm::vec2 UniformDiskSampling();
glm::vec3 CosHemiSampler();
glm::vec3 SphereSampler();
glm::vec3 TriBaryUVSampler();

glm::vec3 BlinnPhongSampler(float Ns);