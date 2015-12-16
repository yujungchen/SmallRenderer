#pragma once

typedef enum {
	Diffuse = 0, 
	Glossy,
	Phong,
	Glass,
	Misc
} MaterialType;

typedef enum {
	MarcoSurface = 0,
	Torrance = 1
} MicroFacetType;

typedef enum {
	MarcoDistribution = 0,
	BlinnPhong = 1
} DistributionType;