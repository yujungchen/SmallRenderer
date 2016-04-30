#include <math.h>
#include <string.h>
#include "utility.h"
#include "core.h"
#include "sampler.h"
#include "dipole.h"

#define EPSILON 0.00001f

Dipole::Dipole(){
	m_Sigma_a = glm::vec3(0.0f, 0.0f, 0.0f);
	m_Sigma_s = glm::vec3(0.0f, 0.0f, 0.0f);
	m_Sigma_t = glm::vec3(0.0f, 0.0f, 0.0f);
	m_Sigma_tr = glm::vec3(0.0f, 0.0f, 0.0f);
	m_Aplha_Prime = glm::vec3(0.0f, 0.0f, 0.0f);
	m_zr = glm::vec3(0.0f, 0.0f, 0.0f);
	m_zv = glm::vec3(0.0f, 0.0f, 0.0f);
	m_N = 1.0f;
	m_Fdr = 0.0f;
	m_A = 0.0f;
}

void Dipole::ComputeFdr(){
	m_Fdr = -1.4399f / (m_N * m_N) + 0.7099f / m_N + 0.6681f + 0.0636f * m_N;
}

void Dipole::DipoleConfiguration(char *_Name, glm::vec3 _Sigma_a, glm::vec3 _Sigma_s, float _N, GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList) {

	m_Sigma_a = _Sigma_a;
	m_Sigma_s = _Sigma_s;
	m_Sigma_t = m_Sigma_s + m_Sigma_a;
	glm::vec3 Temp = m_Sigma_a * m_Sigma_t * 3.0f;
	m_Sigma_tr = glm::vec3(sqrt(Temp.x), sqrt(Temp.y), sqrt(Temp.z));
	m_Aplha_Prime = m_Sigma_s / m_Sigma_t;
	m_zr = glm::vec3(1.0f) / m_Sigma_t;
	strcpy(m_MatName, _Name);
	m_N = _N;
	m_OneoverEta = 1.0f / m_N;
	ComputeFdr();
	m_A = (1.0f + m_Fdr) / (1.0f - m_Fdr);
	m_zv = m_zr * (1.0f + 4.0f / 3.0f * m_A);

	m_model = _model;
	m_bvh = _bvh;
	m_PrimList = _PrimList;

	m_SSNum = 8;
	m_MSNum = 8;

	m_Phase = 0.8f;
	m_PhaseSqr = m_Phase * m_Phase;

	m_ShiftPos = glm::vec3(0.0f);
}

void Dipole::PrintDipoleInfo(){
	printf("\nMaterial %s\n", m_MatName);
	printf("Eta %f\n", m_N);
	printf("Sigma_a (%f %f %f)\n", m_Sigma_a.x, m_Sigma_a.y, m_Sigma_a.z);
	printf("Sigma_s (%f %f %f)\n", m_Sigma_s.x, m_Sigma_s.y, m_Sigma_s.z);
	printf("Sigma_t (%f %f %f)\n", m_Sigma_t.x, m_Sigma_t.y, m_Sigma_t.z);
	printf("Sigma_tr (%f %f %f)\n", m_Sigma_tr.x, m_Sigma_tr.y, m_Sigma_tr.z);
	printf("Sigma_tr (%f %f %f)\n", 1.0f / m_Sigma_tr.x, 1.0f / m_Sigma_tr.y, 1.0f / m_Sigma_tr.z);
	printf("Alpha_Prime (%f %f %f)\n", m_Aplha_Prime.x, m_Aplha_Prime.y, m_Aplha_Prime.z);
	printf("zr (%f %f %f)\n", m_zr.x, m_zr.y, m_zr.z);
	printf("zv (%f %f %f)\n", m_zv.x, m_zv.y, m_zv.z);
	printf("Fdr %f\n", m_Fdr);
	printf("A %f\n", m_A);
	printf("\n");
}

Dipole::~Dipole(){

}

float Dipole::FrDiel(float cosi, float cost, float etai, float etat) {
	float Rparl = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
	float Rperp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
	return (Rparl * Rparl + Rperp * Rperp) / 2.0f;
}

/* Reference code for computing the frensnel term
inline double fresnel(const double cos_theta, const double eta) {
	const double sin_theta_t_sqr = 1.0 / (eta * eta) * (1.0 - cos_theta * cos_theta);
	if (sin_theta_t_sqr >= 1.0) return 1.0;
	const double cos_theta_t = sqrt(1.0 - sin_theta_t_sqr);
	const double r_s = (cos_theta - eta * cos_theta_t) / (cos_theta + eta * cos_theta_t);
	const double r_p = (eta * cos_theta - cos_theta_t) / (eta * cos_theta + cos_theta_t);
	return (r_s * r_s + r_p * r_p) * 0.5;
}
*/

float Dipole::EvalDieletric(float cosi, float eta0, float eta1) {	
	float Dielectric = 0.0f;
	float sint = eta0 / eta1 * sqrtf(fmax(0.0f, 1.0f - cosi * cosi));

	if(sint > 1.0f) {
		return 1.0f;
	}
	else {
		float cost = sqrtf(fmax(0.0f, 1.0f - sint * sint)); 
		Dielectric = FrDiel(fabs(cosi), cost, eta0, eta1); 
	
		return Dielectric;
	}
}

float Dipole::ComputeSSFactor(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos, 
	float &Sigma_s, float &Sigma_t, float &Sigma_tr, glm::vec3 &To, 
	bool &isValid, glm::vec3 &BoundaryPt) {
	
	isValid = true;
	float SS = 0.0f;

	float Seed = RandomNumber();
	float SampleT = -log(1.0f - Seed) / Sigma_t;
	glm::vec3 PSamle = pos + SampleT * To;
	float so = glm::length(PSamle - pos);

	glm::vec3 toL = glm::normalize(l_pos - PSamle);
	
	// Sample Pi
	Vector toL_D = Vector(toL.x, toL.y, toL.z);
	Point O = Point(PSamle.x, PSamle.y, PSamle.z);
	Ray MarchRay = Ray(O, toL_D, EPSILON); 
	Intersection *March_insect = new Intersection();
	
	Point Pi_v = Point(0.0f, 0.0f, 0.0f);
	Normal ni_v = Normal(0.0f, 0.0f, 0.0f);
		
	if(m_bvh->Intersect(MarchRay, March_insect)){
		float t = March_insect->uvt[2];
		Pi_v = MarchRay.o + t * MarchRay.d;
	}
	else {
		delete March_insect;
		return SS;
	}
		
	if(l_pos.x == Pi_v.x && l_pos.y == Pi_v.y && l_pos.z == Pi_v.z) {
		delete March_insect;
		isValid = true;
		return SS;
	}

	Normal n0(m_model->normals[m_model->triangles[March_insect->triId].nindices[0] * 3 + 0], m_model->normals[m_model->triangles[March_insect->triId].nindices[0] * 3 + 1], m_model->normals[m_model->triangles[March_insect->triId].nindices[0] * 3 + 2]);
	Normal n1(m_model->normals[m_model->triangles[March_insect->triId].nindices[1] * 3 + 0], m_model->normals[m_model->triangles[March_insect->triId].nindices[1] * 3 + 1], m_model->normals[m_model->triangles[March_insect->triId].nindices[1] * 3 + 2]);
	Normal n2(m_model->normals[m_model->triangles[March_insect->triId].nindices[2] * 3 + 0], m_model->normals[m_model->triangles[March_insect->triId].nindices[2] * 3 + 1], m_model->normals[m_model->triangles[March_insect->triId].nindices[2] * 3 + 2]);
	ni_v = BaryInterpolationN(n0, n1, n2, March_insect->uvt[0], March_insect->uvt[1]);

	glm::vec3 Pi = glm::vec3(Pi_v.x, Pi_v.y, Pi_v.z);
	glm::vec3 ni = glm::vec3(ni_v.x, ni_v.y, ni_v.z);
	// Sample Pi


	float si = glm::length(Pi - PSamle);
	float LdotNi = fabs(glm::dot(toL, ni));
	float si_prime = si * LdotNi / sqrt(1.0f - m_OneoverEta * m_OneoverEta * (1.0f - LdotNi * LdotNi));

	glm::vec3 LtoPi = glm::normalize(Pi - l_pos);	
	glm::vec3 Ti = glm::refract(LtoPi, ni, 1.0f / m_N);
	Ti = glm::normalize(Ti);

	glm::vec3 w_in = glm::normalize(pos - prev_pos);
	float Ft = 1.0f - EvalDieletric(fabs(glm::dot(w_in, N)), 1.0f, m_N);
	float Fti = 1.0f - EvalDieletric(LdotNi, 1.0f, m_N);

	float Phase = (1.0f - m_PhaseSqr) / (pow(1.0f - 2.0f * m_Phase * glm::dot(Ti, To) + m_PhaseSqr, 1.5) * 4.0f * M_PI);
	float G = fabs(glm::dot(ni, To)) / LdotNi;
	float Sigma_tc = Sigma_t + G * Sigma_t;

	SS = Ft * Fti * Phase * Sigma_s * exp(-1.0f * si_prime * Sigma_t) / (Sigma_tc * Sigma_t);

	delete March_insect;
	return SS;
}

glm::vec3 Dipole::ComputeSingleScatter(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos) {
	
	glm::vec3 SingleScatter = glm::vec3(0.0f);
	glm::vec3 RefractRay = glm::refract(glm::normalize(pos - prev_pos), N, 1.0f / m_N);
	RefractRay = glm::normalize(RefractRay);

	// Refract Ray
	bool isReIntersect = true;
	Point r_O = Point(pos.x, pos.y, pos.z);
	Vector r_dir = Vector(RefractRay.x, RefractRay.y, RefractRay.z);
	Ray TestRay = Ray(r_O, r_dir, EPSILON); 
	Intersection *Test_insect = new Intersection();
	Point Test_P = Point(0.0f, 0.0f, 0.0f);
	
	if(m_bvh->Intersect(TestRay, Test_insect)) {
		float t = Test_insect->uvt[2];
		Test_P = TestRay.o + t * TestRay.d;
	}
	else {
		isReIntersect = false;
		//printf("No Intersect?\n");
	}
	delete Test_insect;
	glm::vec3 BoundaryPt = glm::vec3(Test_P.x, Test_P.y, Test_P.z);
	// Refract Ray

	int ValidSampleNum = 0;
	for(int Idx = 0 ; Idx < m_SSNum ; Idx++) {
		float R = 0.0f;
		float G = 0.0f;
		float B = 0.0f;
		bool isValid = true;
		R = ComputeSSFactor(l_pos, l_N, pos, N, prev_pos, m_Sigma_s.x, m_Sigma_t.x, m_Sigma_tr.x, RefractRay, isValid, BoundaryPt);
		G = ComputeSSFactor(l_pos, l_N, pos, N, prev_pos, m_Sigma_s.y, m_Sigma_t.y, m_Sigma_tr.y, RefractRay, isValid, BoundaryPt);
		B = ComputeSSFactor(l_pos, l_N, pos, N, prev_pos, m_Sigma_s.z, m_Sigma_t.z, m_Sigma_tr.z, RefractRay, isValid, BoundaryPt);

		SingleScatter = SingleScatter + glm::vec3(R, G, B);
	}
	
	SingleScatter = SingleScatter / (float)m_SSNum;

	return SingleScatter;
}

float Dipole::ComputeRd(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos, float &Sigma_s, float &Sigma_t, float &Sigma_tr, 
		bool &isValid, glm::vec3 &U, glm::vec3 &V, 
		float zr, float zv) {
	float Rd = 0.0f;
	
	float offset = 0.005f;
	float Seed = RandomNumber();
	float DiscR = -log(1.0f - Seed) / Sigma_tr;

	float xangle = 2.0f * M_PI * RandomNumber();
	float yangle = 2.0f * M_PI * RandomNumber();
	float xt = DiscR * cos(xangle);
	float yt = DiscR * sin(yangle);
	

	//glm::vec3 SampleP = pos + N * offset + (1.0f / Sigma_tr) * (U * xt + V * yt);

	glm::vec3 SampleP = pos + N * offset + (U * xt + V * yt);
	
	Point Sample = Point(SampleP.x, SampleP.y, SampleP.z);
	Vector ToMed = -1.0f * Vector(N.x, N.y, N.z);

	Ray SampleRay = Ray(Sample, ToMed, EPSILON);
	Intersection *Sample_isect = new Intersection();
	Point Pi_v = Point(0.0f, 0.0f, 0.0f);
	Normal ni_v = Normal(0.0f, 0.0f, 0.0f);
		
	if(m_bvh->Intersect(SampleRay, Sample_isect)){
		float t = Sample_isect->uvt[2];
		Pi_v = SampleRay.o + t * SampleRay.d;
	}
	else {
		delete Sample_isect;
		return Rd;
	}
	Normal n0(m_model->normals[m_model->triangles[Sample_isect->triId].nindices[0] * 3 + 0], m_model->normals[m_model->triangles[Sample_isect->triId].nindices[0] * 3 + 1], m_model->normals[m_model->triangles[Sample_isect->triId].nindices[0] * 3 + 2]);
	Normal n1(m_model->normals[m_model->triangles[Sample_isect->triId].nindices[1] * 3 + 0], m_model->normals[m_model->triangles[Sample_isect->triId].nindices[1] * 3 + 1], m_model->normals[m_model->triangles[Sample_isect->triId].nindices[1] * 3 + 2]);
	Normal n2(m_model->normals[m_model->triangles[Sample_isect->triId].nindices[2] * 3 + 0], m_model->normals[m_model->triangles[Sample_isect->triId].nindices[2] * 3 + 1], m_model->normals[m_model->triangles[Sample_isect->triId].nindices[2] * 3 + 2]);
	ni_v = BaryInterpolationN(n0, n1, n2, Sample_isect->uvt[0], Sample_isect->uvt[1]);

	glm::vec3 Pi = glm::vec3(Pi_v.x, Pi_v.y, Pi_v.z);
	glm::vec3 ni = glm::vec3(ni_v.x, ni_v.y, ni_v.z);
	
	float r = glm::length(Pi - pos);

	float dr = sqrt(r * r + zr * zr);
	float dv = sqrt(r * r + zv * zv);

	float Sigma_tr_dr = Sigma_tr * dr;
	float Sigma_tr_dv = Sigma_tr * dv;
	
	//float pdf = INV_PI * Sigma_tr * exp(-1.0f * Sigma_tr * (xt * xt + yt * yt));

	Rd = (Sigma_tr_dr + 1.0f) * exp(-1.0f * Sigma_tr_dr) * zr / (dr * dr * dr) +
	     (Sigma_tr_dv + 1.0f) * exp(-1.0f * Sigma_tr_dv) * zv / (dv * dv * dv);
	Rd = Rd / Sigma_t;

	Rd = INV_PI * 0.25f * Rd / (Sigma_tr * Sigma_tr * Sigma_tr * exp(-1.0f * Sigma_tr * r));
	
	glm::vec3 toL = glm::normalize(l_pos - Pi);
	float LdotNi = fabs(glm::dot(toL, ni));
	glm::vec3 w_in = glm::normalize(pos - prev_pos);
	
	float Ft = 1.0f - EvalDieletric(fabs(glm::dot(w_in, N)), 1.0f, m_N);
	float Fti = 1.0f - EvalDieletric(LdotNi, 1.0f, m_N);

	Rd = INV_PI * Rd * Ft * Fti;

	delete Sample_isect;
	return Rd;
}
	
glm::vec3 Dipole::ComputeMultipleScatter(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos) {
	glm::vec3 MultipleScatter = glm::vec3(0.0f);

	glm::vec3 U = glm::vec3(0.0f);
	glm::vec3 V = glm::vec3(0.0f);
	LocalBasis(N, &U, &V);

	for(int Idx = 0 ; Idx < m_MSNum ; Idx++) {
		float R = 0.0f;
		float G = 0.0f;
		float B = 0.0f;
		bool isValid = true;
		R = ComputeRd(l_pos, l_N, pos, N, prev_pos, m_Sigma_s.x, m_Sigma_t.x, m_Sigma_tr.x, isValid, U, V, m_zr.x, m_zv.x);
		G = ComputeRd(l_pos, l_N, pos, N, prev_pos, m_Sigma_s.y, m_Sigma_t.y, m_Sigma_tr.y, isValid, U, V, m_zr.y, m_zv.y);
		B = ComputeRd(l_pos, l_N, pos, N, prev_pos, m_Sigma_s.z, m_Sigma_t.z, m_Sigma_tr.z, isValid, U, V, m_zr.z, m_zv.z);
		
		MultipleScatter = MultipleScatter + glm::vec3(R, G, B) * m_Aplha_Prime;
	}
	MultipleScatter = MultipleScatter / (float)m_MSNum;

	return MultipleScatter;
}

glm::vec3 Dipole::Schlick(glm::vec3 V, glm::vec3 N){
	glm::vec3 Schlick = glm::vec3(1.0f);

	float Val = glm::dot(V, N);
	float f = (1.0f - 1.0f / m_N) * (1.0f - 1.0f / m_N) / ((1.0f + 1.0f / m_N) * (1.0f + 1.0f / m_N));
	float Factor = f + (1.0f - f) * Val * Val * Val * Val * Val;
	Factor = (Factor >= 1.0f) ? 1.0f : Factor;


	return Schlick * Factor;
}


glm::vec3 Dipole::ComputeRadiance(glm::vec3 &l_pos, glm::vec3 &l_N, glm::vec3 &pos, glm::vec3 &N, glm::vec3 &prev_pos, glm::vec3 &Kd) {
	
	glm::vec3 Radiance = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 SS = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 MS = glm::vec3(0.0f, 0.0f, 0.0f);
	
	SS = ComputeSingleScatter(l_pos, l_N, pos, N, prev_pos);
	MS = ComputeMultipleScatter(l_pos, l_N, pos, N, prev_pos);
	//glm::vec3 V = glm::normalize(prev_pos - pos);
	//glm::vec3 L = glm::normalize(l_pos - pos);
	//glm::vec3 H = glm::normalize(V + L);

	Radiance = MS * 300000.0f + SS * 3000.0f;// + Schlick(V, H);
	Radiance = Radiance * Kd;
	return Radiance;
}