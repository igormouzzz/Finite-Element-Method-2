#pragma once
#include "Matrix.h"

class Material
{
private:
	double E;
	double Nu;
	Matrix A;
public:
	Material();
	Material(double e, double nu);
	double GetE();
	double GetNu();
	Matrix GetA();
	Matrix CreateMatrixA();
};
