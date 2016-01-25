#include "mcRenderer.h"
#include "sampler.h"
#include "radiometry.h"
#include <time.h>

#define OMP
#define EPSILON 0.00001f

MCRenderer::MCRenderer(GLMmodel *_model, BVHAccel *_bvh, std::vector<Primitive> &_PrimList, 
					   PointLight *_l, AreaLight *_al, Camera *_camera,
					   int _Width, int _Height, float _AspectRatio, int _PathSample, float _FocusDist,
					   int _PathDepth, bool _NEE_Enable){
	m_model = _model;
	m_bvh = _bvh;
	m_l = _l;
	m_al = _al;
	m_camera = _camera;
	m_Width = _Width;
	m_Height = _Height;
	m_PrimList = _PrimList;
	m_AspectRatio = _AspectRatio;
	m_PathSample = _PathSample;
	m_FocusDist = _FocusDist;
	m_PathDepth = _PathDepth;
	m_NEE_Enable = _NEE_Enable;


	// Allocate Frame
	m_ColorImg = (unsigned char*)malloc(m_Width * m_Height * 3 * sizeof(unsigned char));
	m_Img = (glm::vec3 *)malloc(m_Width * m_Height * sizeof(glm::vec3));
	m_DirectImg = (glm::vec3 *)malloc(m_Width * m_Height * sizeof(glm::vec3));
	for(int idx = 0 ; idx < m_Width * m_Height ; idx++){
		m_Img[idx] = glm::vec3(0.0, 0.0, 0.0);
		m_DirectImg[idx] = glm::vec3(0.0, 0.0, 0.0);
		m_ColorImg[idx * 3] = 0;
		m_ColorImg[idx * 3 + 1] = 0;
		m_ColorImg[idx * 3 + 2] = 0;
	}
	// Allocate Frame
	if(m_model->hasLight)
		m_DirectSampleNum = 128;
	else
		m_DirectSampleNum = 64;
	m_Direct = new DirectIllumination(_model, _bvh, _PrimList, _l, _al, _camera, _Width, _Height, m_DirectSampleNum);
}

MCRenderer::~MCRenderer(){

}



void MCRenderer::Render(){
	printf("Rendering...\n");

	// Direct Illumination
	m_Direct->Render(m_DirectImg, m_DirectSampleNum);

	// Progress Illustration
	unsigned int TotalTask = m_PathSample * m_Width * m_Height;
	unsigned int CurrentTask = 0;
	float PercentageInterval = (float)TotalTask / 10.0f;
	float Progress = 0.0f;
	// Progress Illustration
	printf("Indirect Illumination\n");
	printf("Start Rendering\n");


	PathIntegrator *Path = new PathIntegrator(m_model, m_bvh, m_PrimList, m_l, m_al, m_camera, m_NEE_Enable, m_model->hasLight);
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
				
				// Indirect Illumination
				glm::vec3 Radiance = Path->ComputeRadiance(w, h, m_PathDepth);
				int CurrentPxlIdx = h * m_Width + w;
				m_Img[CurrentPxlIdx] = m_Img[CurrentPxlIdx] + Radiance;
			}
		}
	}

	printf("Rendering Done.\n");
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Rendering time:%.3fsec\n\n", elapsed_secs);



	FILE *OutImg_f;
	OutImg_f = fopen("Temp.raw", "wb");

	for(int h = m_Height - 1; h >= 0 ; h--){
		for(int w = 0; w < m_Width ; w++){
			unsigned char r, g, b;
			int CurrentPxlIdx = h * m_Width + w;

			glm::vec3 ClampDirectColor = m_DirectImg[CurrentPxlIdx] / (float)m_DirectSampleNum;
			glm::vec3 ClampColor = m_Img[CurrentPxlIdx] / (float)m_PathSample;

			glm::vec3 TotalColor = ClampDirectColor + ClampColor;

			if(TotalColor.x > 1.0f)
				TotalColor.x = 1.0;
			if(TotalColor.y > 1.0f)
				TotalColor.y = 1.0;
			if(TotalColor.z > 1.0f)
				TotalColor.z = 1.0;

			r = (unsigned char)(int)(TotalColor.x * 255.0f);
			g = (unsigned char)(int)(TotalColor.y * 255.0f);
			b = (unsigned char)(int)(TotalColor.z * 255.0f);
			m_ColorImg[CurrentPxlIdx * 3] = r;
			m_ColorImg[CurrentPxlIdx * 3 + 1] = g;
			m_ColorImg[CurrentPxlIdx * 3 + 2] = b;
			fputc(r, OutImg_f);
			fputc(g, OutImg_f);
			fputc(b, OutImg_f);
		}
	}

	fclose(OutImg_f);


	char ImgName[100];
	sprintf(ImgName, "Focus_%.2f_Sample%d_Trace.ppm", m_FocusDist, m_PathSample);
	frameToPPM(ImgName, m_ColorImg, m_Width, m_Height);
	printf("Output Reference Done...\n\n");
}