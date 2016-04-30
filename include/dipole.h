#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <glm/glm.hpp>
#include "bvh.h"


#pragma once


class Dipole{

public:
	Dipole();
	~Dipole();
	void DipoleConfiguration(char *_Name, glm::vec3 _Sigma_a, glm::vec3 _Sigma_s, float _N, GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList);
	void PrintDipoleInfo();
	void ComputeFdr();

	glm::vec3 ComputeRadiance(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos, glm::vec3 &Kd);
	glm::vec3 ComputeSingleScatter(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos);
	float ComputeSSFactor(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos, float &Sigma_s, float &Sigma_t, float &Sigma_tr, 
		glm::vec3 &To, bool &isValid, glm::vec3 &BoundaryPt);
	glm::vec3 ComputeMultipleScatter(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos);
	float ComputeRd(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos, float &Sigma_s, float &Sigma_t, float &Sigma_tr, 
		bool &isValid, glm::vec3 &U, glm::vec3 &V, float zr, float zv);
	float FrDiel(float cosi, float cost, float etai, float etat);
	float EvalDieletric(float cosi, float eta0, float eta1);
	glm::vec3 Schlick(glm::vec3 V, glm::vec3 N);

	glm::vec3 GetShiftPos();


private:

	GLMmodel *m_model;
	BVHAccel *m_bvh;
	std::vector<Primitive> m_PrimList;

	char m_MatName[20];
	glm::vec3 m_Sigma_a;
	glm::vec3 m_Sigma_s;
	glm::vec3 m_Sigma_t;
	glm::vec3 m_Sigma_tr;
	glm::vec3 m_Aplha_Prime;
	glm::vec3 m_zr;
	glm::vec3 m_zv;
	glm::vec3 m_ShiftPos;
	float m_N;
	float m_OneoverEta;
	float m_Fdr;
	float m_A;
	float m_Phase;
	float m_PhaseSqr;

	int m_SSNum;
	int m_MSNum;
};

