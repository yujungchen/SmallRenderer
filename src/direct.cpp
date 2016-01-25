#include "direct.h"
#include "time.h"

//#define OMP
#define EPSILON 0.0001f

DirectIllumination::DirectIllumination(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
					   PointLight *_l, AreaLight *_al, Camera *_camera, 
					   int _Width, int _Height, int _PathSample){
	m_model = _model;
	m_bvh = _bvh;
	m_l = _l;
	m_al = _al;
	m_camera = _camera;
	m_PrimList = _PrimList;
	m_Width = _Width;
	m_Height = _Height;
	m_PathSample = _PathSample;
	m_useArealight = false;
	if(m_model->hasLight)
		m_useArealight = true;
}

DirectIllumination::~DirectIllumination(){

}


void DirectIllumination::Render(glm::vec3 *m_Img, int SampleNumber){

	printf("Direct Illumination\n");
	// Progress Illustration
	unsigned int TotalTask = m_PathSample * m_Width * m_Height;
	unsigned int CurrentTask = 0;
	float PercentageInterval = (float)TotalTask / 10.0f;
	float Progress = 0.0f;
	printf("Start Rendering\n");
	// Progress Illustration

	Intersection *insect = new Intersection();

	clock_t begin = clock();
	for(int SampleNum = 0 ; SampleNum < m_PathSample ; SampleNum++){	
		for(int h = m_Height - 1; h >= 0 ; h--){
#ifdef OMP
			#pragma omp parallel for
#endif
			for(int w = 0; w < m_Width ; w++){

				// Progress Illustration
				CurrentTask = CurrentTask + 1;
				if((float)CurrentTask / PercentageInterval >= Progress){
					printf("%.0f%%\n", Progress * 10.0f);
					Progress = Progress + 1.0f;
				}
				// Progress Illustration
				int CurrentPxlIdx = h * m_Width + w;

				Ray EyeRay = m_camera->CameraRay(w, h);

				Point hitP(0.0, 0.0, 0.0);
				
				glm::vec3 Contribution = glm::vec3(0.0);
				glm::vec3 Pos = glm::vec3(0.0);
				glm::vec3 N = glm::vec3(0.0);
				glm::vec3 Kd = glm::vec3(0.0);
				glm::vec3 Ks = glm::vec3(0.0);
				glm::vec3 Emission = glm::vec3(0.0);
				float Ns = 0.0f;
				float Eta = 0.0f;

				if(m_bvh->Intersect(EyeRay, insect)){
					float t = insect->uvt[2];
					hitP = EyeRay.o + t * EyeRay.d;
					Pos = glm::vec3(hitP.x, hitP.y, hitP.z);
				}
				else
					continue;

				m_bvh->InterpolateGeoV2(EyeRay, insect, Pos, N, Kd, Ks, Emission, Ns, Eta, m_PrimList);

				// If hit light source, directly return color
				if(isLight(Emission)){
					Contribution = Kd;
					m_Img[CurrentPxlIdx] = m_Img[CurrentPxlIdx] + Contribution;
					continue;
				}

				// Direct Illumination
				// Sample light source
				glm::vec3 l_Pos = glm::vec3(0.0f);
				glm::vec3 l_N = glm::vec3(0.0f);
				glm::vec3 l_emission = glm::vec3(0.0f);
				if(m_useArealight){
					l_emission = m_al->sampleL(l_Pos, l_N);
				}
				else{
					l_Pos = m_l->getlpos();
					l_N = glm::normalize(Pos - l_Pos);
					l_emission = m_l->sampleL();
				}

				// Cast shadow ray
				glm::vec3 dir2Light = l_Pos - Pos;
				dir2Light = glm::normalize(dir2Light);

				// Backface of the light source
				if(glm::dot(l_N, dir2Light) >= 0.0f)
					continue;

				Vector shaod_dir = Vector(dir2Light.x, dir2Light.y, dir2Light.z);
				Ray ShadowRay(hitP, shaod_dir, EPSILON, glm::length(l_Pos - Pos) - 0.0001f);
				if(m_bvh->IntersectP(ShadowRay))
					continue;
				
				if(Eta == 0.0f){
					Contribution = EvalPhongBRDF(m_camera->m_CameraPos, Pos, l_Pos, N, Kd, Ks, Ns) * ComputeG(Pos, l_Pos, N, l_N) * l_emission;
				}
				
				m_Img[CurrentPxlIdx] = m_Img[CurrentPxlIdx] + Contribution;
			}
		}
	}
	printf("Rendering Done.\n");\
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Rendering time:%.3fsec\n\n", elapsed_secs);
}