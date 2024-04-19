#pragma once
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
using namespace std;

#include "Header.h"

#include <thread>
#include <omp.h>
#include <vector>

#include <sycl/sycl.hpp>
#if FPGA_HARDWARE || FPGA_EMULATOR || FPGA_SIMULATOR
#include <sycl/ext/intel/fpga_extensions.hpp>
#endif

using row = vector<double>;

class Matrix
{
protected:
	vector<row> a;
	int N, M;
	void SetZero();
	void Get_matr(vector<row> matr, int n, vector<row> temp_matr, int indRow, int indCol);
	double Det0(vector<row> matr, int n);
public:
	Matrix();
	Matrix(int N, int M);
	Matrix(vector<row> rows);
	Matrix(const Matrix& b);
	double GetForIndicies(int i, int j) { return a[i][j]; }
	size_t GetN();
	void Set(int i, int j, double value);
	Matrix& operator=(const Matrix& b);
	Matrix operator+(const Matrix& b);
	Matrix operator-(const Matrix& b);
	Matrix operator*(double k);
	Matrix operator*(const Matrix& b);
	vector<double> operator*(vector<double>& b);
	const vector<double> operator*(const vector<double>& b);
	Matrix T();
	void Symmetric();
	double SumInRow(int k);
	void SumOfRows();
	double Det();
	Matrix Inversed();
	Matrix Inv2();
	void Clean();
	~Matrix() { Clean(); }
	friend class Element;
	friend class TriangularElement;
	friend class QuadElement;
	friend class Material;
	friend class DivisionToLocalsTri;
	friend class DivisionToLocalsTri2;
	friend Matrix MadeGlobalStiffnessMatrix(int element_count, int node_count, vector<Matrix> matrices, vector<vc> nums);
	friend ostream& operator<<(ostream& cout, const Matrix& b);

	vector<double> Gauss(vector<double>& b);
	vector<double> CG(vector<double>& b);
	vector<double> CG2(vector<double>& b, double tolerance);
	vector<double> CG3(vector<double>& b);

	void ToCSR(ofstream& f);
	double SumOfComponentsForProduct(vector<double>& b);
	void UnitRowAndColumn(int i, int j, int n);
	void ZeroRowAndColumn(int i, int j);
};