#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "sampler.h"
#include "define.h"

#pragma once

glm::vec3 ComputeTorranceMicroFacetBRDF(glm::vec3 &Pos0, glm::vec3 &Pos1, glm::vec3 &Pos2, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, 
	float Ns, DistributionType &Distribution, float &Roughness,
	glm::vec3 &MicroNormal);


glm::vec3 SampleMicroNormal(glm::vec3 &N, float Ns, DistributionType &Distribution);