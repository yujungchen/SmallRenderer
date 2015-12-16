#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "configure.h"
#include "bvh.h"
#include "light.h"
#include "camera.h"
#include "mcRenderer.h"


GLMmodel *model;
BVHAccel *bvh;
std::vector<Primitive> triList;

int main(int argc, char **argv){
	
	bool isConfiguration = false;

	if(argc == 2){
		if(strcmp(argv[1], "configure.txt") != 0){
			printf("\nPlease use the correct configuration or directly input parameter.\n");
			printf("[Usage] ./SmallRenderer [Model] [ImageWidth] [SampleNumber] [FocusDistance*10]\n\n");
			return 0;
		}
		else
			isConfiguration = true;
	}
	else{
		if(argc < 5){
	    	printf("\n[Usage]: ./SmallRenderer [Model] [ImageWidth] [SampleNumber] [FocusDistance*10]\n\n");
	    	return 0;
		}
		else
			isConfiguration = false;
	}
	
	srand (time(NULL));
	
	//Scene Configuration
	string m_Model;
	float m_SceneScale = 0.0f;
	
	//Camera Configuration
	glm::vec3 m_CameraPos(0.0f, 0.0f, 0.0f);
	glm::vec3 m_CameraUp(0.0f, 0.0f, 0.0f);
	glm::vec3 m_CameraLookat(0.0f, 0.0f, 0.0f);
	float m_Fstop = 0.0f;
	bool m_AllinFocus = false;
	float m_FocalLength = 0.0f; 
	float m_Aperture = 0.0f;
	float m_SensorWidth = 0.0f; 
	float m_SensorHeight = 0.0f; 
	float m_AspectRatio = 0.0f;
	int m_ResWidth = 0;
	int m_ResHeight = 0; 
	float m_FovV = 0.0f;
	float m_FovH = 0.0f;
	int m_SPP = 0;
	int m_PathDepth = 1;
	bool m_NEE_Enable = false;


	//Light Configuration
	bool m_TestLight = true;
	glm::vec3 m_TestLightPos(0.0f, 0.0f, 0.0f);
	glm::vec3 m_TestLightEmission(0.0f, 0.0f, 0.0f);

	ReadConfigure(argv[1],
				  //Scene Configuration
				  m_Model,
				  m_SceneScale, 
				  //Camera Configuration
				  m_CameraPos, m_CameraUp, m_CameraLookat,
				  m_Fstop, m_AllinFocus, m_FocalLength, m_Aperture, 
				  m_SensorWidth, m_SensorHeight, m_AspectRatio, 
				  m_ResWidth, m_ResHeight, m_FovV, m_FovH,
				  m_SPP,
				  //Light Configuration
				  m_TestLight, m_TestLightPos, m_TestLightEmission, 
				  m_PathDepth, m_NEE_Enable);

	float m_FocusDist = m_CameraLookat.z;	

	if(isConfiguration){
		cout<<"Loading\t"<<m_Model<<endl;
		const char *ModelPath = m_Model.c_str();
		printf("Reading obj file...\n");
		model = glmReadOBJ(ModelPath);
		printf("Finish reading.\n");
	}
	else{
		// Command as Input
		char ModelPath[80];
		char ObjName[80];
		strcpy(ModelPath, "../Model/");
		sprintf(ObjName, "%s/%s.obj", argv[1], argv[1]);
		strcat(ModelPath, ObjName);
		printf("\nLoading %s\n", ModelPath);
		printf("Reading obj file...\n");
		model = glmReadOBJ(ModelPath);
		printf("Finish reading.\n");

		// Setting Frame
		m_ResWidth = atoi(argv[2]);
		m_SPP = atoi(argv[3]);
		m_FocusDist = atof(argv[4]);
		m_AspectRatio = 1.f;
		m_ResHeight = m_ResWidth * m_AspectRatio;
		m_SceneScale = 10.0f;

		printf("Image Size = %d x %d\n",m_ResWidth, m_ResHeight);
	}


	// Read Scene
	UnitizeScene(m_SceneScale);
	for (uint32_t i = 0; i < model->numtriangles; i++) {
		Primitive *pri = new Primitive(i);
		triList.push_back(*pri);
	}
	// Read Scene

	// Build BVH
	printf("\nBuilding AS...\n");
	bvh = new BVHAccel();
	printf("Finish building AS.\n\n");
	// Build BVH

	TestLight *l = new TestLight(m_TestLightPos, m_TestLightEmission);
	Camera *cam = new Camera(m_CameraPos, m_CameraLookat, m_CameraUp,
							 m_Fstop, m_FocalLength, m_Aperture, m_SensorWidth, m_SensorHeight,
							 m_AspectRatio, m_ResWidth, m_ResHeight, m_FovV, m_FovH, m_AllinFocus);


	MCRenderer *Renderer = new MCRenderer(model, bvh, triList, l, cam, m_ResWidth, m_ResHeight, m_AspectRatio, m_SPP, m_FocusDist, m_PathDepth, m_NEE_Enable);
	
	Renderer->Render();

	return 0;
}