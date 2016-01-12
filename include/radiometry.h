#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "sampler.h"


#pragma once

typedef enum{
	Diffuse = 0, 
	Glossy,
	Phong,
	Misc
}MaterialType;

glm::vec3 EvalPhongBRDF(glm::vec3 Pos0, glm::vec3 Pos1, glm::vec3 Pos2, glm::vec3 N, glm::vec3 Kd, glm::vec3 Ks, float Ns);
float ComputeG(glm::vec3 Pos0, glm::vec3 Pos1, glm::vec3 N0, glm::vec3 N1);
float ComputeG2PLight(glm::vec3 Pos0, glm::vec3 Pos1, glm::vec3 N0);
MaterialType DetermineMat(glm::vec3 Kd, glm::vec3 Ks);

class PathVtx{

public:
	PathVtx(){
		Pos = glm::vec3(0.0f);
		N = glm::vec3(0.0f);
		Kd = glm::vec3(0.0f);
		Ks = glm::vec3(0.0f);
		Ns = 0.0f;
		Mat = Misc;
	}
	PathVtx(glm::vec3 _Pos, glm::vec3 _N, glm::vec3 _Kd, glm::vec3 _Ks, float _Ns){
		Pos = _Pos;
		N = _N;
		Kd = _Kd;
		Ks = _Ks;
		Ns = _Ns;
		Mat = DetermineMat(_Kd, _Ks);
	}
	~PathVtx();

	glm::vec3 Pos;
	glm::vec3 N;
	glm::vec3 Kd;
	glm::vec3 Ks;
	float Ns;
	MaterialType Mat;
};

glm::vec3 LocalDirSampling(glm::vec3 PrevPos, glm::vec3 Pos, glm::vec3 N, glm::vec3 Kd, glm::vec3 Ks, float Ns, double &Pdf_W_proj, glm::vec3 &Throughput);