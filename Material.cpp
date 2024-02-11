#include "Material.h"

Material::Material() {}
Material::Material(double e, double nu)
{
	E = e; Nu = nu;
	A = CreateMatrixA();
}
double Material::GetE() { return E; }
double Material::GetNu() { return Nu; }
Matrix Material::GetA() { return A; }
Matrix Material::CreateMatrixA()
{
	Matrix A(3, 3);
	double C = E / (1 - Nu * Nu);
	A.a[0][0] = A.a[1][1] = 1;
	A.a[0][1] = A.a[1][0] = Nu;
	A.a[2][2] = (1 - Nu) / 2;
	A = A * C;
	//cout << A << endl;
	return A;
}