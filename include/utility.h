#include <math.h>
#include <algorithm>
#include <float.h>
#include <stdint.h>
#include "stdlib.h"
#include "stdio.h"
#include "glm.h"
#include <glm/glm.hpp>


#pragma once


#ifndef M_PI
#define M_PI       3.14159265358979323846f
#endif

#ifndef INV_PI
#define INV_PI     0.31830988618379067154f
#endif

#ifndef INV_TWOPI
#define INV_TWOPI  0.15915494309189533577f
#endif

#ifndef INFINITY
#define INFINITY 1000000
#endif




// Global Inline Functions
inline float Lerp(float t, float v1, float v2) {
    return (1.f - t) * v1 + t * v2;
}


inline float Clamp(float val, float low, float high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Mod(int a, int b) {
    int n = int(a/b);
    a -= n*b;
    if (a < 0) a += b;
    return a;
}


inline float Radians(float deg) {
    return ((float)M_PI/180.f) * deg;
}


inline float Degrees(float rad) {
    return (180.f/(float)M_PI) * rad;
}


inline float Log2(float x) {
    static float invLog2 = 1.f / logf(2.f);
    return logf(x) * invLog2;
}


inline int Floor2Int(float val);
inline int Log2Int(float v) {
    return Floor2Int(Log2(v));
}


inline bool IsPowerOf2(int v) {
    return v && !(v & (v - 1));
}


inline uint32_t RoundUpPow2(uint32_t v) {
    v--;
    v |= v >> 1;    v |= v >> 2;
    v |= v >> 4;    v |= v >> 8;
    v |= v >> 16;
    return v+1;
}


inline int Floor2Int(float val) {
    return (int)floorf(val);
}


inline int Round2Int(float val) {
    return Floor2Int(val + 0.5f);
}


inline int Float2Int(float val) {
    return (int)val;
}


inline int Ceil2Int(float val) {
    return (int)ceilf(val);
}


class Vector;
class Point;
class Normal;
class Ray;
class BBox;

// Geometry Declarations
class Vector {
public:
    // Vector Public Methods
    Vector() { x = y = z = 0.f; }
    Vector(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz) {
    }
    explicit Vector(Point p);
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.

    
    Vector operator=(Vector v) {
        x = v.x; y = v.y; z = v.z;
        return *this;
    }
    Vector operator+(Vector v) {
        return Vector(x + v.x, y + v.y, z + v.z);
    }
    
    Vector operator+=(Vector v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector operator-(Vector v) {
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    
    Vector operator-=(Vector v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Vector operator*(float f) { return Vector(f*x, f*y, f*z); }
    
    Vector operator*=(float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Vector operator/(float f) {
        float inv = 1.f / f;
        return Vector(x * inv, y * inv, z * inv);
    }
    
    Vector operator/=(float f) {
        float inv = 1.f / f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    Vector operator-() { return Vector(-x, -y, -z); }
    float operator[](int i) {
        return (&x)[i];
    }
    
   
    float LengthSquared() { return x*x + y*y + z*z; }
    float Length() { return sqrtf(LengthSquared()); }
    explicit Vector(Normal n);

    bool operator==(Vector v) {
        return x == v.x && y == v.y && z == v.z;
    }
    bool operator!=(Vector v) {
        return x != v.x || y != v.y || z != v.z;
    }

    // Vector Public Data
    float x, y, z;
};


class Point {
public:
    // Point Public Methods
    Point() { x = y = z = 0.f; }
    Point(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz) {
    }
    
    Point operator=(Point p) {
        x = p.x; y = p.y; z = p.z;
        return *this;
    }
    Point operator+(Vector v) {
        return Point(x + v.x, y + v.y, z + v.z);
    }
    
    Point operator+=(Vector v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    Vector operator-(Point p) {
        return Vector(x - p.x, y - p.y, z - p.z);
    }
    
    Point operator-(Vector v) {
        return Point(x - v.x, y - v.y, z - v.z);
    }
    
    Point operator-=(Vector v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }
    Point operator+=(Point p) {
        x += p.x; y += p.y; z += p.z;
        return *this;
    }
    Point operator+(Point p) {
        return Point(x + p.x, y + p.y, z + p.z);
    }
    Point operator* (float f) {
        return Point(f*x, f*y, f*z);
    }
    Point operator*=(float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Point operator/ (float f) {
        float inv = 1.f/f;
        return Point(inv*x, inv*y, inv*z);
    }
    Point operator/=(float f) {
        float inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
     float operator[](int i) const {
        return (&x)[i];
    }
    
    float operator[](int i) {
        return (&x)[i];
    }
    

    bool operator==(Point p) {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(Point p) {
        return x != p.x || y != p.y || z != p.z;
    }

    // Point Public Data
    float x, y, z;
};


class Normal {
public:
    // Normal Public Methods
    Normal() { x = y = z = 0.f; }
    Normal(float xx, float yy, float zz)
        : x(xx), y(yy), z(zz) {
    }
    Normal operator-() {
        return Normal(-x, -y, -z);
    }
    Normal operator+ (Normal n) {
        return Normal(x + n.x, y + n.y, z + n.z);
    }
    
    Normal operator+=(Normal n) {
        x += n.x; y += n.y; z += n.z;
        return *this;
    }
    Normal operator- (Normal n) {
        return Normal(x - n.x, y - n.y, z - n.z);
    }
    
    Normal operator-=(Normal n) {
        x -= n.x; y -= n.y; z -= n.z;
        return *this;
    }
    Normal operator*(float f) {
        return Normal(f*x, f*y, f*z);
    }
    
    Normal operator*=(float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    Normal operator/(float f) {
        float inv = 1.f/f;
        return Normal(x * inv, y * inv, z * inv);
    }
    
    Normal operator/=(float f) {
        float inv = 1.f/f;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }
    float LengthSquared() { return x*x + y*y + z*z; }
    float Length()        { return sqrtf(LengthSquared()); }
    
    
    Normal operator=(Normal n) {
        x = n.x; y = n.y; z = n.z;
        return *this;
    }
    explicit Normal(Vector v)
      : x(v.x), y(v.y), z(v.z) {
    }
    float operator[](int i) {
        return (&x)[i];
    }
    

    bool operator==(Normal n) {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(Normal n) {
        return x != n.x || y != n.y || z != n.z;
    }

    // Normal Public Data
    float x, y, z;
};


class Ray {
public:
    // Ray Public Methods
    Ray() : mint(0.f), maxt(INFINITY), depth(0) { }
    Ray(Point origin, Vector direction,
        float start, float end = INFINITY, int d = 0)
        : o(origin), d(direction), mint(start), maxt(end), depth(d) { }
    Ray(Point origin, Vector direction, Ray parent,
        float start, float end = INFINITY)
        : o(origin), d(direction), mint(start), maxt(end),
          depth(parent.depth+1) { }
    Point operator()(float t) { return o + d * t; }
  

    // Ray Public Data
    Point o;
    Vector d;
    mutable float mint, maxt;
    int depth;
};


class BBox {
public:
     // BBox Public Methods
    BBox() {
        pMin = Point( INFINITY,  INFINITY,  INFINITY);
        pMax = Point(-INFINITY, -INFINITY, -INFINITY);
    }
    BBox(Point p) : pMin(p), pMax(p) { }
    BBox(Point p1, Point p2) {
        pMin = Point(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
        pMax = Point(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
    }
    friend BBox Union(BBox b, Point p);
    friend BBox Union(BBox b, BBox b2);
    bool Overlaps(BBox &b) {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }
    bool Inside(Point pt) {
        return (pt.x >= pMin.x && pt.x <= pMax.x &&
                pt.y >= pMin.y && pt.y <= pMax.y &&
                pt.z >= pMin.z && pt.z <= pMax.z);
    }
    void Expand(float delta) {
        pMin -= Vector(delta, delta, delta);
        pMax += Vector(delta, delta, delta);
    }
    float SurfaceArea() {
        Vector d = pMax - pMin;
        return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    float Volume() {
        Vector d = pMax - pMin;
        return d.x * d.y * d.z;
    }
    int MaximumExtent() {
        Vector diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }
    Point &operator[](int i);
    Point Lerp(float tx, float ty, float tz) {
        return Point(::Lerp(tx, pMin.x, pMax.x), ::Lerp(ty, pMin.y, pMax.y),
                     ::Lerp(tz, pMin.z, pMax.z));
    }
    Vector Offset(Point p) {
        return Vector((p.x - pMin.x) / (pMax.x - pMin.x),
                      (p.y - pMin.y) / (pMax.y - pMin.y),
                      (p.z - pMin.z) / (pMax.z - pMin.z));
    }
    void BoundingSphere(Point *c, float *rad);
    bool IntersectP(Ray ray, float *hitt0 = NULL, float *hitt1 = NULL);

    bool operator==(BBox b) {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(BBox b) {
        return b.pMin != pMin || b.pMax != pMax;
    }

    // BBox Public Data
    Point pMin, pMax;
};



// Geometry Inline Functions
inline Vector::Vector(Point p)
    : x(p.x), y(p.y), z(p.z) {
}


inline Vector operator*(float f, Vector v) { return v*f; }
inline float Dot(Vector v1, Vector v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


inline float AbsDot(Vector v1, Vector v2) {
    return fabsf(Dot(v1, v2));
}


inline Vector Cross(Vector v1, Vector v2) {
    float v1x = v1.x, v1y = v1.y, v1z = v1.z;
    float v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector Cross(Vector v1, Normal v2) {
    float v1x = v1.x, v1y = v1.y, v1z = v1.z;
    float v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector Cross(Normal v1, Vector v2) {
    float v1x = v1.x, v1y = v1.y, v1z = v1.z;
    float v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector((v1y * v2z) - (v1z * v2y),
                  (v1z * v2x) - (v1x * v2z),
                  (v1x * v2y) - (v1y * v2x));
}


inline Vector Normalize(Vector v) { return v / v.Length(); }
inline void CoordinateSystem(Vector v1, Vector *v2, Vector *v3) {
    if (fabsf(v1.x) > fabsf(v1.y)) {
        float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
        *v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
    }
    else {
        float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
        *v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
    }
    *v3 = Cross(v1, *v2);
}


inline float Distance(Point p1, Point p2) {
    return (p1 - p2).Length();
}


inline float DistanceSquared(Point p1, Point p2) {
    return (p1 - p2).LengthSquared();
}


inline Point operator*(float f, Point p) {
    return p*f;
}


inline Normal operator*(float f, Normal n) {
    return Normal(f*n.x, f*n.y, f*n.z);
}


inline Normal Normalize(Normal n) {
    return n / n.Length();
}


inline Vector::Vector(Normal n)
  : x(n.x), y(n.y), z(n.z) {
}


inline float Dot(Normal n1, Vector v2) {
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}


inline float Dot(Vector v1, Normal n2) {
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}


inline float Dot(Normal n1, Normal n2) {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}


inline float AbsDot(Normal n1, Vector v2) {
    return fabsf(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}


inline float AbsDot(Vector v1, Normal n2) {
    return fabsf(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}


inline float AbsDot(Normal n1, Normal n2) {
    return fabsf(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}



inline Point &BBox::operator[](int i) {
    return (&pMin)[i];
}

inline float RAND() {
	return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

void frameToPPM(char *f, unsigned char *data, int width, int height);

Point CosineHemisphereSampler(float u1, float u2);

float EvalPSNR(unsigned char* Img, unsigned char* GoldenImg, int ImgWidth, int ImgHeight);

typedef struct _EnergyInfo {
    double PkgEnergy;
    double PP0Energy;
    double PP1Energy;
    double DRAMEnergy;
} EnergyInfo;


typedef struct _PixelPos {
    int x;
    int y;
} PixelPos;


inline glm::vec3 BarycentricInterpolation(glm::vec3 P0, glm::vec3 P1, glm::vec3 P2, glm::vec3 Coef){
    glm::vec3 Result = glm::vec3(0.0, 0.0, 0.0);
    Result = Coef.x * P0 + Coef.y * P1 + Coef.z * P2;
    return Result;
}


inline void LocalBasis(glm::vec3 Normal, glm::vec3 *BiNormal, glm::vec3 *BiTangent){
    Normal = glm::normalize(Normal);

    glm::vec3 TempBiNormal;
    glm::vec3 TempBiTangent;

    if(abs(Normal.x) > abs(Normal.z)){
        TempBiNormal.x = -Normal.y;
        TempBiNormal.y = Normal.x;
        TempBiNormal.z = 0.0;
    }
    else{
        TempBiNormal.x = 0.0;
        TempBiNormal.y = -Normal.z;
        TempBiNormal.z = Normal.y;
    }
    TempBiNormal = glm::normalize(TempBiNormal);

    
    TempBiTangent = glm::cross(TempBiNormal, Normal);
    *BiNormal = TempBiTangent;
    *BiTangent = -1.0f * TempBiNormal;
}

inline Normal BaryInterpolationN(Normal n0, Normal n1, Normal n2, float u, float v){
    return u * n0 + v * n1 + (1 - u - v) * n2;
}

inline void TBN(Vector Normal, Vector* BiNormal, Vector* BiTangent){
    Normal = Normalize(Normal);

    Vector TempBiNormal;
    Vector TempBiTangent;

    if(abs(Normal.x) > abs(Normal.z)){
        TempBiNormal.x = -Normal.y;
        TempBiNormal.y = Normal.x;
        TempBiNormal.z = 0.0;
    }
    else{
        TempBiNormal.x = 0.0;
        TempBiNormal.y = -Normal.z;
        TempBiNormal.z = Normal.y;
    }
    TempBiNormal = Normalize(TempBiNormal);

    *BiNormal = TempBiNormal;
    TempBiTangent = Cross(TempBiNormal, Normal);
    *BiTangent = TempBiTangent;
}

inline bool isLight(glm::vec3 &Emission){
    bool isL = false;

    if(Emission.x > 0.0f || Emission.y > 0.0f || Emission.z > 0.0f)
        isL = true;

    return isL;
}

inline bool isZero(glm::vec3 Var){
    bool Zero = false;

    if(Var.x == 0.0f && Var.y == 0.0f && Var.z == 0.0f)
        Zero = true;

    return Zero;
}