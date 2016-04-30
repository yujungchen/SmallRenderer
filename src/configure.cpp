#include "configure.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

void ReadConfigure(char *SceneFile,
				   //Scene Configuration
				   string &m_Model,
				   float &m_SceneScale, 
				   //Camera Configuration
				   glm::vec3 &m_CameraPos, glm::vec3 &m_CameraUp, glm::vec3 &m_CameraLookat,
				   float &m_Fstop, bool &m_AllinFocus, float &m_FocalLength, float &m_Aperture, 
				   float &m_SensorWidth, float &m_SensorHeight, float &m_AspectRatio, 
				   int &m_ResWidth, int &m_ResHeight, float &m_FovV, float &m_FovH, 
				   int &m_SPP, 
				   //Light Configuration
				   bool &m_TestLight, glm::vec3 &m_TestLightPos, glm::vec3 &m_TestLightEmission,
				   int &m_PathDepth, bool &m_NEE_Enable
				   ){


	ifstream f_in(SceneFile);

	if(!f_in.is_open())
		printf("Can't Open File\n");
	
	printf("\nRendering Configuration\n");
	while(!f_in.eof()){
		
		char Item[100];
		
		f_in.getline(Item, sizeof(Item), ' ');

		if(strcmp(Item, "Model") == 0){
			printf("\nScene Configuration\n");
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			//m_Model = Content;
			m_Model = Content;
			printf("%s %s\n", Item, Content);
		}
		else if(strcmp(Item, "SceneScale") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_SceneScale = atof(Content);
			printf("%s %.2fm\n", Item, m_SceneScale);
		}	
		else if(strcmp(Item, "CameraPos") == 0){
			printf("\nCamera Configuration\n");
			char X[10];
			char Y[10];
			char Z[10];
			f_in.getline(X, sizeof(X), ' ');
			m_CameraPos.x = atof(X);
			f_in.getline(Y, sizeof(Y), ' ');
			m_CameraPos.y = atof(Y);
			f_in.getline(Z, sizeof(Z), '\n');
			m_CameraPos.z = atof(Z);
			printf("Camera Pos (%.2f %.2f %.2f)\n", m_CameraPos.x, m_CameraPos.y, m_CameraPos.z);
		}
		else if(strcmp(Item, "CameraUpVector") == 0){
			char X[10];
			char Y[10];
			char Z[10];
			f_in.getline(X, sizeof(X), ' ');
			m_CameraUp.x = atof(X);
			f_in.getline(Y, sizeof(Y), ' ');
			m_CameraUp.y = atof(Y);
			f_in.getline(Z, sizeof(Z), '\n');
			m_CameraUp.z = atof(Z);
			printf("Camera Up Vector (%.2f %.2f %.2f)\n", m_CameraUp.x, m_CameraUp.y, m_CameraUp.z);
		}
		else if(strcmp(Item, "CameraLookat") == 0){
			char X[10];
			char Y[10];
			char Z[10];
			f_in.getline(X, sizeof(X), ' ');
			m_CameraLookat.x = atof(X);
			f_in.getline(Y, sizeof(Y), ' ');
			m_CameraLookat.y = atof(Y);
			f_in.getline(Z, sizeof(Z), '\n');
			m_CameraLookat.z = atof(Z);
			printf("Camera Lookat (%.2f %.2f %.2f)\n", m_CameraLookat.x, m_CameraLookat.y, m_CameraLookat.z);
		}
		else if(strcmp(Item, "Fstop") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_Fstop = atof(Content);
			printf("%s %.1ff\n", Item, m_Fstop);
		}
		else if(strcmp(Item, "AllinFocus") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			if(atoi(Content) == 1){
				m_AllinFocus = true;
				printf("Use Pinhole Camera\n");
			}
			else if(atoi(Content) == 0){
				m_AllinFocus = false;
				printf("Use Thin Lens Camera\n");
			}
			else{
				printf("What the hell!\n");
				break;
			}
			
		}
		else if(strcmp(Item, "FocalLength") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_FocalLength = atof(Content) / 1000.0f;
			printf("%s %.0fmm\n", Item, m_FocalLength * 1000.0f);
			m_Aperture = m_FocalLength / m_Fstop;
			printf("Aperture %.0fmm\n", m_Aperture * 1000.0f);
		}
		else if(strcmp(Item, "SensorWidth") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_SensorWidth = atof(Content) / 1000.0f;
			printf("%s %.0fmm\n", Item, m_SensorWidth * 1000.0f);
		}
		else if(strcmp(Item, "SensorHeight") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_SensorHeight = atof(Content) / 1000.0f;
			printf("%s %.0fmm\n", Item, m_SensorHeight * 1000.0f);

			// Compute FoV
			m_FovV = atan(m_SensorHeight * 0.5f / m_FocalLength) * 2.0f * 180.0f / M_PI;
			printf("FoV (Vertical) %.2f Degree\n", m_FovV);

			m_FovH = atan(m_SensorWidth * 0.5f / m_FocalLength) * 2.0f * 180.0f / M_PI;
			printf("FoV (Horizontal) %.2f Degree\n", m_FovH);

			// Compute Aspect Ratio
			m_AspectRatio = m_SensorWidth / m_SensorHeight;
			printf("Aspect Ratio %.2f\n", m_AspectRatio);

		}
		else if(strcmp(Item, "ResWidth") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_ResWidth = atoi(Content);
			m_ResHeight = (int)((float)m_ResWidth / m_AspectRatio);
			printf("Resolution %dx%d\n", m_ResWidth, m_ResHeight);
		}
		else if(strcmp(Item, "SPP") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_SPP = atoi(Content);
			printf("SPP %d\n", m_SPP);
		}
		// Light Configuration
		else if(strcmp(Item, "TestLight") == 0){
			printf("\nLight Configuration\n");
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			int Val = atoi(Content);
			if(Val > 0){
				m_TestLight = true;
				printf("Use Test Light\n");
			}
			else{
				m_TestLight = false;
				printf("Use Real Light\n");
			}
		}
		else if(strcmp(Item, "TestLightPos") == 0){
			
			char X[10];
			char Y[10];
			char Z[10];
			f_in.getline(X, sizeof(X), ' ');
			m_TestLightPos.x = atof(X);
			f_in.getline(Y, sizeof(Y), ' ');
			m_TestLightPos.y = atof(Y);
			f_in.getline(Z, sizeof(Z), '\n');
			m_TestLightPos.z = atof(Z);
			printf("Test Light Pos (%.2f %.2f %.2f)\n", m_TestLightPos.x, m_TestLightPos.y, m_TestLightPos.z);
			
		}
		else if(strcmp(Item, "TestLightEmission") == 0){
			
			char X[10];
			char Y[10];
			char Z[10];
			f_in.getline(X, sizeof(X), ' ');
			m_TestLightEmission.x = atof(X);
			f_in.getline(Y, sizeof(Y), ' ');
			m_TestLightEmission.y = atof(Y);
			f_in.getline(Z, sizeof(Z), '\n');
			m_TestLightEmission.z = atof(Z);
			printf("Test Light Emission (%.2f %.2f %.2f)\n", m_TestLightEmission.x, m_TestLightEmission.y, m_TestLightEmission.z);
			
		}
		else if(strcmp(Item, "PathDepth") == 0){
			printf("\nPath Integrator Configuration\n");
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			m_PathDepth = atoi(Content);
			printf("Path Depth %d\n", m_PathDepth);
		}	
		else if(strcmp(Item, "NEE_Enable") == 0){
			char Content[100];
			f_in.getline(Content, sizeof(Content), '\n');
			if(atoi(Content) == 1){
				m_NEE_Enable = true;
				printf("Enable Next Event Estimation\n");
			}
			else if(atoi(Content) == 0){
				m_NEE_Enable = false;
				printf("Disable Next Event Estimation\n");
			}
			else{
				printf("What the hell!\n");
				break;
			}
		}	
		else{
			//printf("Seriously!?\n");
			printf("\nPlease Check Your Scene File. (%s)\n", Item);
			break;
		}
	}


	printf("\n");

}