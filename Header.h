#pragma once
#include <iostream>
#include <fstream>
#include "json.hpp"
using namespace std;
using json = nlohmann::json;
static const std::string base64_chars =
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789+/";
inline bool is_base64(unsigned char c);
void base64_decode(const std::string& encoded_string, std::string& decoded_string);
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <thread>
#include <omp.h>
#include <typeinfo>

#define COUNT_OF_NODES 4


struct float6
{
	float x, y, z, w, u, v;
};
typedef struct float6 float6;

struct double6
{
	double x, y, z, w, u, v;
};
typedef struct double6 double6;
struct double3
{
	double x, y, z;
	void operator=(double3& b)
	{
		x = b.x; y = b.y; z = b.z;
	}
};
typedef struct double3 double3;
struct double2
{
	double x, y;
	void operator=(double2& b)
	{
		x = b.x; y = b.y;
	}
	void operator=(double3& b)
	{
		x = b.x; y = b.y;
	}
};
typedef struct double2 double2;
struct int3
{
	int n[3];
};
typedef struct int3 int3;
struct int4
{
	int n[4];
};
typedef struct int4 int4;

#if COUNT_OF_NODES == 3
using vc = int3;
#else
using vc = int4;
#endif

double det2(double a11, double a12, double a21, double a22);
double det3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);

class Matrix;
class Force;
class Restraint;
class Element;
class Material;
int Read(vector<double3>& pt_list, vector<vc>& hex_list, Force& F, vector<Restraint>& R, double& E, double& Nu, int& number_of_nodes_of_elem);
Matrix MadeGlobalStiffnessMatrix(int element_count, int node_count, vector<Matrix> matrices, vector<vc> nums);
int Task();
//#define FLOAT_TYPE

#ifdef FLOAT_TYPE
typedef float FTYPE;
typedef float2 FTYPE2;
typedef float3 FTYPE3;
typedef float4 FTYPE4;
typedef float6 FTYPE6;
#else
typedef double FTYPE;
//typedef double2 FTYPE2;
typedef double3 FTYPE3;
//typedef double4 FTYPE4;
typedef double6 FTYPE6;
#endif