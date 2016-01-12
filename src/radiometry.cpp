#include "radiometry.h"
#include "utility.h"


glm::vec3 EvalPhongBRDF(glm::vec3 Pos0, glm::vec3 Pos1, glm::vec3 Pos2, glm::vec3 N, glm::vec3 Kd, glm::vec3 Ks, float Ns){
	glm::vec3 BRDF = glm::vec3(0.0f);
	
	float GlossyFactor = 0.0f;

	glm::vec3 v2Tov1 = glm::normalize(Pos1 - Pos2);
	glm::vec3 v1Tov0 = glm::normalize(Pos0 - Pos1);
	glm::vec3 v2Reflect = glm::normalize(glm::reflect(v2Tov1, N));
	float CosAlpha = glm::max(glm::dot(v1Tov0, v2Reflect), 0.0f);
	GlossyFactor = (Ns + 2.0f) * 0.5f * INV_PI * pow(CosAlpha, Ns);

	BRDF = Kd * INV_PI + GlossyFactor * Ks;
	return BRDF;
}

float ComputeG(glm::vec3 Pos0, glm::vec3 Pos1, glm::vec3 N0, glm::vec3 N1){
	float GTerm = 0.0f;
	glm::vec3 vConnect = glm::normalize(Pos0 - Pos1);
	float Cos0 = glm::dot(N0, vConnect);
	float Cos1 = glm::dot(N1, vConnect);
	float Dist = glm::length(Pos0 - Pos1);

	GTerm = fabs(Cos0 * Cos1) / (Dist * Dist);

	return GTerm;
}

float ComputeG2PLight(glm::vec3 Pos0, glm::vec3 Pos1, glm::vec3 N0){
	float GTerm = 0.0f;
	glm::vec3 vConnect = glm::normalize(Pos1 - Pos0);
	float Cos0 = std::max(glm::dot(N0, vConnect), 0.0f);
	float Dist = glm::length(Pos0 - Pos1);

	GTerm = Cos0 / (Dist * Dist);

	return GTerm;
}

MaterialType DetermineMat(glm::vec3 Kd, glm::vec3 Ks, float Eta){
	MaterialType Mat = Misc;

	if(Eta != 0.0f){
		Mat = Glass;
	}
	else{
		if((Kd.x > 0.0 || Kd.y > 0.0 || Kd.z > 0.0) && 
		   (Ks.x > 0.0 || Ks.y > 0.0 || Ks.z > 0.0))
		{
			Mat = Phong;
		}
		else if((Kd.x > 0.0 || Kd.y > 0.0 || Kd.z > 0.0) && 
		   		(Ks.x == 0.0 && Ks.y == 0.0 && Ks.z == 0.0)){
			Mat = Diffuse;
		}
		else if((Kd.x == 0.0 && Kd.y == 0.0 && Kd.z == 0.0) && 
		   		(Ks.x > 0.0 || Ks.y > 0.0 || Ks.z > 0.0)){
			Mat = Glossy;
		}
		else{
			Mat = Misc;
			printf("Unknown Material!\n");
		}
	}

	return Mat;
}
Vector Reflect(Vector i, Vector n){
	return i - 2.0f * n * Dot(n,i);
}
	
glm::vec3 LocalDirSampling(glm::vec3 PrevPos, glm::vec3 Pos, glm::vec3 N, glm::vec3 Kd, glm::vec3 Ks, float Ns, float Eta, double &Pdf_W_proj, glm::vec3 &Throughput){
	glm::vec3 LocalDir = glm::vec3(0.0f);

	MaterialType Mat = DetermineMat(Kd, Ks, Eta);

	glm::vec3 U = glm::vec3(0.0f);
	glm::vec3 V = glm::vec3(0.0f);
	
	switch(Mat){
		case Diffuse:{
			glm::vec3 LocalP = glm::vec3(0.0f);
			LocalP = CosHemiSampler();
			LocalBasis(N, &U, &V);
			LocalDir = V * LocalP.x + U * LocalP.y + N * LocalP.z;

			// Compute Projected Angular Pdf
			Pdf_W_proj = 1.0 / M_PI;	// (Cos / PI) * (1 / Cos)
			
			Throughput = Kd;
			break;
		}
		case Glossy:{
			glm::vec3 InVec = glm::normalize(Pos - PrevPos);
			glm::vec3 ReflectVec = glm::reflect(InVec, N);
			LocalBasis(ReflectVec, &U, &V);
			
			float u1 = RandomNumber();
			float u2 = RandomNumber();
			float CosTheta = powf(u1, 1.0f / (Ns + 1.0f) );
			float Phi = u2 * 2.0f * (float)M_PI;
			float SinTheta = sqrt(1.0f - CosTheta * CosTheta);
  
			float x = cos(Phi) * SinTheta;
			float y = sin(Phi) * SinTheta;
			float z = CosTheta;
			LocalDir = U * x + V * y + ReflectVec * z;

			// Compute Projected Angular Pdf
			Pdf_W_proj = (Ns + 1.0f) / (2.0f * M_PI) * pow(fabs(glm::dot(LocalDir, ReflectVec)), Ns - 1.0f);	//Note: The proj. angular pdf needs to be checked.
			// (n + 1) / (2 * PI) * (pow(Cos, Ns)) / Cos = (n + 1) / (2 * PI) * (pow(Cos, Ns - 1))
			float GlossyConst = ((Ns + 2.0f) / (Ns + 1.0f)) * fabs(glm::dot(LocalDir, ReflectVec));				//Note: The cosine term needs to be checked.
			
			Throughput = Ks * GlossyConst;
			break;
		}
		case Glass:{
			glm::vec3 InVec = glm::normalize(Pos - PrevPos);
			glm::vec3 ReflectVec = glm::reflect(InVec, N);

			glm::vec3 N_l = glm::dot(InVec, N) < 0.0f ? N : N * (-1.0f);
			bool into = glm::dot(N_l, N) > 0.0f;
			float nnt = into ? 1.0f / Eta : Eta / 1.0f;
			float ddn = glm::dot(InVec, N_l);
			float cos2t = 1.0f - nnt * nnt * (1.0f - ddn * ddn);

			if(cos2t < 0.0f){
				LocalDir = glm::reflect(InVec, N * (-1.0f));
			}
			else{
				glm::vec3 t_Dir = InVec * nnt - N * (into ? 1.0f : -1.0f) * (ddn * nnt + sqrtf(cos2t));
				t_Dir = glm::normalize(t_Dir);
				float a = Eta - 1.0f;
				float b = Eta + 1.0f;
				float R0 = a * a / (b * b);
				float c = 1.0f - (into ? (-ddn) : glm::dot(t_Dir, N)); 
				float Re = R0 + (1.0f - R0) * c * c * c * c * c;
				float Tr = 1.0f - Re;
				float P = 0.25f + 0.5f * Re;
				float RP = Re / P;
				float TP = Tr / (1.0f - P);
				float RndNum = RandomNumber();

				if(RndNum < P){
					if(into)
						LocalDir = ReflectVec;
					else
						LocalDir = glm::reflect(InVec, N * (-1.0f));
				}
				else{
					LocalDir = t_Dir;
				}

			}


			break;
		}
		case Phong:{
			break;
		}
		case Misc:{
			printf("What the hell!\n");
			break;
		}
	}

	return LocalDir;
}