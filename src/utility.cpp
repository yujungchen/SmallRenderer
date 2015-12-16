#include "utility.h"




BBox Union(BBox b, Point p) {
	BBox ret = b;
	ret.pMin.x = std::min(b.pMin.x, p.x);
	ret.pMin.y = std::min(b.pMin.y, p.y);
	ret.pMin.z = std::min(b.pMin.z, p.z);
	ret.pMax.x = std::max(b.pMax.x, p.x);
	ret.pMax.y = std::max(b.pMax.y, p.y);
	ret.pMax.z = std::max(b.pMax.z, p.z);
	return ret;
}
BBox Union(BBox b, BBox b2) {
	BBox ret;
	ret.pMin.x = std::min(b.pMin.x, b2.pMin.x);
	ret.pMin.y = std::min(b.pMin.y, b2.pMin.y);
	ret.pMin.z = std::min(b.pMin.z, b2.pMin.z);
	ret.pMax.x = std::max(b.pMax.x, b2.pMax.x);
	ret.pMax.y = std::max(b.pMax.y, b2.pMax.y);
	ret.pMax.z = std::max(b.pMax.z, b2.pMax.z);
	return ret;
};


void frameToPPM(char *f, unsigned char *data, int width, int height) {
  FILE *out = fopen(f, "wb");
  if (!out) {
	  fprintf(stderr, "***Error: Unable to capture frame.  fopen() failed!!\n");
	  return;
  }
  
  fprintf(out, "P6\n");  
  fprintf(out, "%d %d\n", width, height);
  fprintf(out, "%d\n", 255); 
  
  for (int y = height-1; y >= 0; y--)
	  fwrite(data+(3*y*width), 1, 3*width, out);
  fprintf(out, "\n");
  fclose(out);
}


Point CosineHemisphereSampler(float u1, float u2){
	float r = sqrt(u1);
	float theta = 2 * M_PI * u2;
 
	float x = r * cos(theta);
	float y = r * sin(theta);
 
	return Point(x, y, sqrt(std::max(0.0f, 1 - u1)));
}

float EvalPSNR(unsigned char* Img, unsigned char* GoldenImg, int ImgWidth, int ImgHeight){
	
	float PSNR = 0.0f;
	float MSE_R = 0.0f;
	float MSE_G = 0.0f;
	float MSE_B = 0.0f;
	for(int j = 0 ; j < ImgHeight ; j++){
		for(int i = 0 ; i < ImgWidth ; i++){
			unsigned char TmpR, TmpG, TmpB;
			unsigned char GoldenR, GoldenG, GoldenB;
			int Offset = j * ImgWidth + i;
			
			TmpR = Img[Offset * 3 + 0];
			TmpG = Img[Offset * 3 + 1];
			TmpB = Img[Offset * 3 + 2];
			
			int Offset1 = (511 - j) * ImgWidth + i;
			GoldenR = GoldenImg[Offset1 * 3 + 0];
			GoldenG = GoldenImg[Offset1 * 3 + 1];
			GoldenB = GoldenImg[Offset1 * 3 + 2];
			
		  
			float SqarR = (float)((int)TmpR - (int)GoldenR) * (float)((int)TmpR - (int)GoldenR);
			float SqarG = (float)((int)TmpG - (int)GoldenG) * (float)((int)TmpG - (int)GoldenG);
			float SqarB = (float)((int)TmpB - (int)GoldenB) * (float)((int)TmpB - (int)GoldenB);
			
			MSE_R = MSE_R + SqarR;
			MSE_G = MSE_G + SqarG;
			MSE_B = MSE_B + SqarB;
		  
		}
	  
	}
	
	MSE_R = MSE_R / (ImgHeight * ImgWidth);
	MSE_G = MSE_G / (ImgHeight * ImgWidth);
	MSE_B = MSE_B / (ImgHeight * ImgWidth);
	
	float MSE = (MSE_R + MSE_G + MSE_B) / 3.0;
	
	PSNR = 10.0 * log10( (255.0 * 255.0) / MSE);
	
	return PSNR;
}



