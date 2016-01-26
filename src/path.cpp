#include "path.h"
#define EPSILON 0.0001f

PathIntegrator::PathIntegrator(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
							   PointLight *_l, AreaLight *_al, Camera *_camera, bool _NEE_Enable, bool _UseAreaLight){
	m_model = _model;
	m_bvh = _bvh;
	m_l = _l;
	m_al = _al;
	m_camera = _camera;
	m_PrimList = _PrimList;
	m_NEE_Enable = _NEE_Enable;
	m_UseAreaLight = _UseAreaLight;
}

PathIntegrator::~PathIntegrator(){

}

glm::vec3 PathIntegrator::NEE(glm::vec3 &Pos, glm::vec3 &PrevPos, glm::vec3 &N, glm::vec3 &Kd, glm::vec3 &Ks, float Ns, float Eta){

	glm::vec3 NEERad = glm::vec3(0.0f);
	
	if(Eta > 0.0)
		return NEERad;

	// Sample the light source
	glm::vec3 l_Pos = glm::vec3(0.0f);
	glm::vec3 l_N = glm::vec3(0.0f);
	glm::vec3 l_emission = glm::vec3(0.0f);

	if(m_UseAreaLight){
		l_emission = m_al->sampleL(l_Pos, l_N);
	}
	else{
		l_Pos = m_l->getlpos();
		l_N = glm::normalize(Pos - l_Pos);
		l_emission = m_l->sampleL();
	}

	glm::vec3 dir2Light = l_Pos - Pos;
	dir2Light = glm::normalize(dir2Light);

	// Backface of the light source
	if(glm::dot(l_N, dir2Light) >= 0.0f)
		return NEERad;

	Vector shaod_dir = Vector(dir2Light.x, dir2Light.y, dir2Light.z);
	Point hitP = Point(Pos.x, Pos.y, Pos.z);
	Ray ShadowRay(hitP, shaod_dir, EPSILON, glm::length(l_Pos - Pos) - 0.0001f);
	
	if(m_bvh->IntersectP(ShadowRay))
		return NEERad;

	NEERad = EvalPhongBRDF(PrevPos, Pos, l_Pos, N, Kd, Ks, Ns) * ComputeG(Pos, l_Pos, N, l_N) * l_emission;
	
	return NEERad;
}

glm::vec3 PathIntegrator::ComputeRadiance(int sample_x, int sample_y, int PathDepth){
	
	glm::vec3 Rad = glm::vec3(0.0f);
	
	Point hitP(0.0, 0.0, 0.0);
	Intersection *insect = new Intersection();
	glm::vec3 Pos = glm::vec3(0.0);
	glm::vec3 N = glm::vec3(0.0);
	glm::vec3 Kd = glm::vec3(0.0);
	glm::vec3 Ks = glm::vec3(0.0);
	float Ns = 0.0f;
	float Eta = 0.0f;
	
	glm::vec3 PrevN = glm::vec3(0.0);

	Ray RaySeg = m_camera->CameraRay(sample_x, sample_y);
	glm::vec3 SampleDir = glm::vec3(0.0f);
	glm::vec3 Throughput = glm::vec3(1.0f, 1.0f, 1.0f);
	glm::vec3 PrevPos = m_camera->m_CameraPos;

	// Pdf Parameter
	double Pdf_A = 1.0;
	double Prev_Pdf_W_proj = 1.0;

	glm::vec3 VtxThroughput = glm::vec3(1.0f, 1.0f, 1.0f);

	for(int Depth = 0 ; Depth < PathDepth ; Depth++){

		Pos = glm::vec3(0.0);
		N = glm::vec3(0.0);
		Kd = glm::vec3(0.0);
		Ks = glm::vec3(0.0);
		Ns = 0.0f;
		Eta = 0.0f;

		if(m_bvh->Intersect(RaySeg, insect)){
			float t = insect->uvt[2];
			hitP = RaySeg.o + t * RaySeg.d;
			Pos = glm::vec3(hitP.x, hitP.y, hitP.z);
		}
		else
			break;

		VtxThroughput = glm::vec3(1.0f, 1.0f, 1.0f);
		double Current_Pdf_W_proj = 1.0;	

		glm::vec3 Emission = glm::vec3(0.0f, 0.0f, 0.0f);
		m_bvh->InterpolateGeoV2(RaySeg, insect, Pos, N, Kd, Ks, Emission, Ns, Eta, m_PrimList);

		// Hit the light source
		if(isLight(Emission)){
			glm::vec3 dir2Light = Pos - PrevPos;
			dir2Light = glm::normalize(dir2Light);

			// Backface of the light source
			if(glm::dot(N, dir2Light) < 0.0f){
				float GTerm = ComputeG(PrevPos, Pos, PrevN, N);
				Rad = Rad + Throughput * Kd * Emission * GTerm;
			}
			break;
		}


		if(m_NEE_Enable){
			if(Depth > 0){
				// Next event estimation
				glm::vec3 Contribution = NEE(Pos, PrevPos, N, Kd, Ks, Ns, Eta);
				Rad = Rad + Throughput * Contribution;
			}
		}
		else{
			if(Depth == (PathDepth - 1)){
				// Next event estimation
				glm::vec3 Contribution = NEE(Pos, PrevPos, N, Kd, Ks, Ns, Eta);		
				Rad = Rad + Throughput * Contribution;
			}
		}

		// Keep sampling the ray
		SampleDir = LocalDirSampling(PrevPos, Pos, N, Kd, Ks, Ns, Eta, Current_Pdf_W_proj, VtxThroughput);
		Point P = Point(Pos.x, Pos.y, Pos.z);
		Vector Dir = Vector(SampleDir.x, SampleDir.y, SampleDir.z);
		RaySeg = Ray(P, Dir, EPSILON);

		Throughput = Throughput * VtxThroughput;


		if(Depth > 0){
			float GTerm = ComputeG(PrevPos, Pos, PrevN, N);
			Pdf_A = Pdf_A * (Prev_Pdf_W_proj * GTerm);
		}

		PrevPos = Pos;
		PrevN = N;
		Prev_Pdf_W_proj = Current_Pdf_W_proj;
	}

	delete insect;

	return Rad;
}