#include "core.h"
#include "glm.h"

extern GLMmodel* model;

class Primitive;

void UnitizeScene(float Scale) {
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
    float minX = FLT_MAX, minY = FLT_MAX, minZ = FLT_MAX;
    float maxX = FLT_MIN, maxY = FLT_MIN, maxZ = FLT_MIN;
    // Walk through vertex array, find max & min (x, y, z)
    for (int i = 1; i < model->numvertices+1; i++) {
		Point vertex = Point(model->vertices[i*3], model->vertices[i*3+1], model->vertices[i*3+2]);
        minX = std::min(minX, vertex.x);
        maxX = std::max(maxX, vertex.x);
        minY = std::min(minY, vertex.y);
		maxY = std::max(maxY, vertex.y);
        minZ = std::min(minZ, vertex.z);
        maxZ = std::max(maxZ, vertex.z);
    }
	
    // Find the center of the object & the width of the object
    float ctrX = 0.5*(minX+maxX);
    float ctrY = 0.5*(minY+maxY); 
    float ctrZ = 0.5*(minZ+maxZ);
    float deltaX = 0.5*(maxX-minX), deltaY = 0.5*(maxY-minY), deltaZ = 0.5*(maxZ-minZ);
    float delta = std::max(deltaX, std::max(deltaY, deltaZ))+0.0001;	// add little extend here!!

    // Walk through the vertex array and update the positions
     for (int i = 1; i < model->numvertices+1; i++) {
		Point vertex = Point(model->vertices[i*3], model->vertices[i*3+1], model->vertices[i*3+2]);
		model->vertices[i*3] = Scale * (vertex.x - ctrX) / delta;
		model->vertices[i*3+1] = Scale * (vertex.y - ctrY) / delta;
		model->vertices[i*3+2] = Scale * (vertex.z - ctrZ) / delta;
    }
}


bool BrouteForceTracing(Ray &ray, Intersection *insect) {
	bool hasHit = false;
	for(int i = 0; i < model->numtriangles; i++) {
		Primitive pri(i);
		if( pri.Intersect(ray, insect) )
			hasHit = true;
	}
	return hasHit;
}



Primitive::~Primitive() { }


void Primitive::GetKdKsNs(float u, float v, Vector &Kd, Vector &Ks, float &Ns) {
	GLMtriangle tri = model->triangles[primitiveId];
	GLMmaterial mtl = model->materials[tri.matId];
	if(mtl.diffuse_map[0] != '\0') {
		int width = mtl.width;
		int height = mtl.height;
		
		int index = Clamp(int(v*height+0.5), 0, height-1) * width + Clamp(int(u*width+0.5), 0, width-1);
		Vector interpolatedTex(mtl.texelData[index*3]/255.f, mtl.texelData[index*3+1]/255.f, mtl.texelData[index*3+2]/255.f);
		Kd = Vector(mtl.diffuse[0]*interpolatedTex.x, mtl.diffuse[1]*interpolatedTex.y, mtl.diffuse[2]*interpolatedTex.z);
		Ks = Vector(mtl.specular[0]*interpolatedTex.x, mtl.specular[1]*interpolatedTex.y, mtl.specular[2]*interpolatedTex.z);
		Ns = mtl.shininess;
	} else {
		Kd = Vector(mtl.diffuse[0], mtl.diffuse[1], mtl.diffuse[2]);
		Ks = Vector(mtl.specular[0], mtl.specular[1], mtl.specular[2]);
		Ns = mtl.shininess;
	}
}

Vector Primitive::GetKd(float u, float v) {
	GLMtriangle tri = model->triangles[primitiveId];
	GLMmaterial mtl = model->materials[tri.matId];
	if(mtl.diffuse_map[0] != '\0') {
		int width = mtl.width;
		int height = mtl.height;
		
		int index = Clamp(int(v*height+0.5), 0, height-1) * width + Clamp(int(u*width+0.5), 0, width-1);
		Vector interpolatedTex(mtl.texelData[index*3]/255.f, mtl.texelData[index*3+1]/255.f, mtl.texelData[index*3+2]/255.f);

		return Vector(mtl.diffuse[0]*interpolatedTex.x, mtl.diffuse[1]*interpolatedTex.y, mtl.diffuse[2]*interpolatedTex.z);
	} else {
		return Vector(mtl.diffuse[0], mtl.diffuse[1], mtl.diffuse[2]);
		//return Vector(0, 0, 0);
	}
}

BBox Primitive::GetWorldBound() {
	GLMtriangle tri = model->triangles[primitiveId];
	Point p1 = Point(model->vertices[tri.vindices[0]*3], model->vertices[tri.vindices[0]*3+1], model->vertices[tri.vindices[0]*3+2]);
	Point p2 = Point(model->vertices[tri.vindices[1]*3], model->vertices[tri.vindices[1]*3+1], model->vertices[tri.vindices[1]*3+2]);
	Point p3 = Point(model->vertices[tri.vindices[2]*3], model->vertices[tri.vindices[2]*3+1], model->vertices[tri.vindices[2]*3+2]);
	return Union(BBox(p1, p2), p3);
};


bool Primitive::Intersect(Ray &ray, Intersection *isect) {
	GLMtriangle tri = model->triangles[primitiveId];
	Point V0 = Point(model->vertices[tri.vindices[0]*3], model->vertices[tri.vindices[0]*3+1], model->vertices[tri.vindices[0]*3+2]);
    Point V1 = Point(model->vertices[tri.vindices[1]*3], model->vertices[tri.vindices[1]*3+1], model->vertices[tri.vindices[1]*3+2]);
    Point V2 = Point(model->vertices[tri.vindices[2]*3], model->vertices[tri.vindices[2]*3+1], model->vertices[tri.vindices[2]*3+2]);
	Vector E1 = V1 - V0;
	Vector E2 = V2 - V0;
	Vector P = Cross(ray.d, E2);
	
	float PdotE1 = Dot(P, E1);
	if(PdotE1 == 0.f)
		 return false;
	float inv_PdotE1 = 1.0 / PdotE1;
	
	Vector T = ray.o - V0;
	float u = Dot(P, T) * inv_PdotE1;
	if(u < 0.0 || u > 1.0)	
		return false;
		
	Vector Q = Cross(T, E1);
	float v = Dot(Q, ray.d) * inv_PdotE1;
	if(v < 0.0 || (u+v) > 1.0)	
		return false;
		
	float t = Dot(Q, E2) * inv_PdotE1;
	if(t < ray.mint || t > ray.maxt) 	
		return false;
	
	 // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    float uvs[3][2];
    GetUVs(uvs);



    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - u - v;
    float tu = b0*uvs[0][0] + u*uvs[1][0] + v*uvs[2][0];
    float tv = b0*uvs[0][1] + u*uvs[1][1] + v*uvs[2][1];

	if(tu > 1.f)	tu = tu - int(tu);
	if(tv > 1.f)	tv = tv - int(tv);
	if(tu < 0.f)	tu = tu - int(tu) + 1;
	if(tv < 0.f)	tv = tv - int(tv) + 1;
	isect->triId = primitiveId;
	isect->uvt[0] = tu;
	isect->uvt[1] = tv;
	isect->uvt[2] = t;
	ray.maxt = t;
	return true;
	

}


bool Primitive::IntersectP(Ray &ray) {
	GLMtriangle tri = model->triangles[primitiveId];
	Point V0 = Point(model->vertices[tri.vindices[0]*3], model->vertices[tri.vindices[0]*3+1], model->vertices[tri.vindices[0]*3+2]);
    Point V1 = Point(model->vertices[tri.vindices[1]*3], model->vertices[tri.vindices[1]*3+1], model->vertices[tri.vindices[1]*3+2]);
    Point V2 = Point(model->vertices[tri.vindices[2]*3], model->vertices[tri.vindices[2]*3+1], model->vertices[tri.vindices[2]*3+2]);
	Vector E1 = V1 - V0;
	Vector E2 = V2 - V0;
	Vector P = Cross(ray.d, E2);
	
	float PdotE1 = Dot(P, E1);
	if(PdotE1 == 0.)
		 return false;
	float inv_PdotE1 = 1.0 / PdotE1;
	
	Vector T = ray.o - V0;
	float u = Dot(P, T) * inv_PdotE1;
	if(u < 0.0 || u > 1.0)	
		return false;
		
	Vector Q = Cross(T, E1);
	float v = Dot(Q, ray.d) * inv_PdotE1;
	if(v < 0.0 || (u+v) > 1.0)	
		return false;
		
	float t = Dot(Q, E2) * inv_PdotE1;
	if(t < ray.mint || t > ray.maxt) 	
		return false;
	
	ray.maxt = t;
	return true;

}

void Primitive::GetUVs(float uv[3][2]) {
		GLMtriangle tri = model->triangles[primitiveId];
		if (tri.tindices[0] == 0 && tri.tindices[1] == 0 && tri.tindices[2] == 0) {
            uv[0][0] = 0.; uv[0][1] = 0.;
            uv[1][0] = 1.; uv[1][1] = 0.;
            uv[2][0] = 1.; uv[2][1] = 1.;
        }
        else {
            uv[0][0] = model->texcoords[tri.tindices[0]*2];
            uv[0][1] = model->texcoords[tri.tindices[0]*2+1];
            uv[1][0] = model->texcoords[tri.tindices[1]*2];
            uv[1][1] = model->texcoords[tri.tindices[1]*2+1];
            uv[2][0] = model->texcoords[tri.tindices[2]*2];
            uv[2][1] = model->texcoords[tri.tindices[2]*2+1];
        }
    }


    void Primitive::GetKdKsNsEta(float u, float v, Vector &Kd, Vector &Ks, float &Ns, float &Eta) {
	GLMtriangle tri = model->triangles[primitiveId];
	GLMmaterial mtl = model->materials[tri.matId];
	if(mtl.diffuse_map[0] != '\0') {
		int width = mtl.width;
		int height = mtl.height;
		
		int index = Clamp(int(v*height+0.5), 0, height-1) * width + Clamp(int(u*width+0.5), 0, width-1);
		Vector interpolatedTex(mtl.texelData[index*3]/255.f, mtl.texelData[index*3+1]/255.f, mtl.texelData[index*3+2]/255.f);
		Kd = Vector(mtl.diffuse[0]*interpolatedTex.x, mtl.diffuse[1]*interpolatedTex.y, mtl.diffuse[2]*interpolatedTex.z);
		Ks = Vector(mtl.specular[0]*interpolatedTex.x, mtl.specular[1]*interpolatedTex.y, mtl.specular[2]*interpolatedTex.z);
		Ns = mtl.shininess;
		Eta = mtl.eta;
	} else {
		Kd = Vector(mtl.diffuse[0], mtl.diffuse[1], mtl.diffuse[2]);
		Ks = Vector(mtl.specular[0], mtl.specular[1], mtl.specular[2]);
		Ns = mtl.shininess;
		Eta = mtl.eta;
	}
}

void Primitive::GetKdKsNsEtaEmission(float u, float v, Vector &Kd, Vector &Ks, Vector &Emission, float &Ns, float &Eta) {
	GLMtriangle tri = model->triangles[primitiveId];
	GLMmaterial mtl = model->materials[tri.matId];
	if(mtl.diffuse_map[0] != '\0') {
		int width = mtl.width;
		int height = mtl.height;
		
		int index = Clamp(int(v*height+0.5), 0, height-1) * width + Clamp(int(u*width+0.5), 0, width-1);
		Vector interpolatedTex(mtl.texelData[index*3]/255.f, mtl.texelData[index*3+1]/255.f, mtl.texelData[index*3+2]/255.f);
		Kd = Vector(mtl.diffuse[0]*interpolatedTex.x, mtl.diffuse[1]*interpolatedTex.y, mtl.diffuse[2]*interpolatedTex.z);
		Ks = Vector(mtl.specular[0]*interpolatedTex.x, mtl.specular[1]*interpolatedTex.y, mtl.specular[2]*interpolatedTex.z);
		Ns = mtl.shininess;
		Eta = mtl.eta;
	} else {
		Kd = Vector(mtl.diffuse[0], mtl.diffuse[1], mtl.diffuse[2]);
		Ks = Vector(mtl.specular[0], mtl.specular[1], mtl.specular[2]);
		Emission = Vector(mtl.emissive[0], mtl.emissive[1], mtl.emissive[2]);
		Ns = mtl.shininess;
		Eta = mtl.eta;
	}
}