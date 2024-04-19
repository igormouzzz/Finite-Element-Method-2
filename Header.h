#pragma once
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
#define eps 0.01


#define COUNT_OF_NODES 3

struct float6
{
	float x, y, z, w, u, v;
};
typedef struct float6 float6;

struct double6
{
	double n[6];
	void operator=(double6& b)
	{
		for (int i = 0; i < 6; i++) n[i] = b.n[i];
	}
};
typedef struct double6 double6;
struct double3_
{
	double x, y, z;
	void operator=(double3_& b)
	{
		x = b.x; y = b.y; z = b.z;
	}
};
typedef struct double3_ double3_;
struct double2_
{
	double x, y;
	void operator=(double2_& b)
	{
		x = b.x; y = b.y;
	}
	void operator=(double3_& b)
	{
		x = b.x; y = b.y;
	}
};
typedef struct double2_ double2_;
struct int3_
{
	int n[3];
};
typedef struct int3_ int3_;
struct int4_
{
	int n[4];
};
typedef struct int4_ int4_;

#if COUNT_OF_NODES == 3
using vc = int3_;
#else
using vc = int4_;
#endif

double det2(double a11, double a12, double a21, double a22);
double det3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);

class Matrix;
class Force;
class Restraint;
class Element;
class Material;
class DivisionToLocals;
int Read(vector<double3_>& pt_list, vector<vc>& hex_list, Force& F, vector<Restraint>& R, double& E, double& Nu, int& number_of_nodes_of_elem, vector<vector<int>>& list_nodes_with_elem_nums);
Matrix MadeGlobalStiffnessMatrix(int element_count, int node_count, vector<Matrix> matrices, vector<vc> nums);
int Task();
int Task2();
int Task3();
int Example();

double norm_square(vector<double>& b);
double norm_l1(vector<double>& b);
double scalar_product(vector<double>& a, vector<double>& b);
//#define FLOAT_TYPE

#ifdef FLOAT_TYPE
typedef float FTYPE;
typedef float2 FTYPE2;
typedef float3 FTYPE3;
typedef float4 FTYPE4;
typedef float6 FTYPE6;
#else
typedef double FTYPE;
//typedef double2_ FTYPE2;
typedef double3_ FTYPE3;
//typedef double4 FTYPE4;
typedef double6 FTYPE6;
#endif