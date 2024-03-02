#include "Element.h"

QuadElement::QuadElement() {}
QuadElement::QuadElement(vector<double3> a, Material m)
{
	for (int i = 0; i < a.size(); i++)
	{
		coord[i].x = a[i].x;
		coord[i].y = a[i].y;
	}
	mat = m;
}
QuadElement::QuadElement(vector<double2> a, Material m)
{
	for (int i = 0; i < a.size(); i++)
	{
		coord[i].x = a[i].x;
		coord[i].y = a[i].y;
	}
	mat = m;
}
QuadElement QuadElement::operator=(QuadElement b)
{
	for (int i = 0; i < 4; i++)
	{
		this->coord[i] = b.coord[i];
	}
	return *this;
}
double* QuadElement::FF(double ksi, double eta)
{
	double n[4];
	n[0] = 0.25 * (1 - ksi - eta + ksi * eta);
	n[1] = 0.25 * (1 + ksi - eta - ksi * eta);
	n[2] = 0.25 * (1 + ksi + eta + ksi * eta);
	n[3] = 0.25 * (1 - ksi + eta - ksi * eta);
	return n;
}
Matrix QuadElement::DerFF(double ksi, double eta)
{
	Matrix gradFF(2, 4);
	gradFF.a[0][0] = 0.25 * (-1 + eta);	//N1 ksi
	gradFF.a[1][0] = 0.25 * (-1 + ksi); //N1 eta
	gradFF.a[0][1] = 0.25 * (1 - eta); //N2 ksi
	gradFF.a[1][1] = 0.25 * (-1 - ksi); //N2 eta
	gradFF.a[0][2] = 0.25 * (1 + eta); //N3 ksi
	gradFF.a[1][2] = 0.25 * (1 + ksi); //N3 eta
	gradFF.a[0][3] = 0.25 * (-1 - eta); //N4 ksi
	gradFF.a[1][3] = 0.25 * (1 - ksi); //N4 eta
	return gradFF;
}
Matrix QuadElement::Jac(double ksi, double eta)
{
	Matrix M = DerFF(ksi, eta);
	//cout << "M = " << endl << M << endl;
	Matrix MatOfCoord(4, 2);
	MatOfCoord.a[0][0] = coord[0].x; MatOfCoord.a[0][1] = coord[0].y;
	MatOfCoord.a[1][0] = coord[1].x; MatOfCoord.a[1][1] = coord[1].y;
	MatOfCoord.a[2][0] = coord[2].x; MatOfCoord.a[2][1] = coord[2].y;
	MatOfCoord.a[3][0] = coord[3].x; MatOfCoord.a[3][1] = coord[3].y;
	//cout << "MatOfCoord = " << endl << MatOfCoord << endl;
	return M * MatOfCoord;
}
Matrix QuadElement::DerFFGlob(double ksi, double eta)
{
	Matrix derFF = DerFF(ksi, eta);
	Matrix J = Jac(ksi, eta);
	//cout << "J = " << endl << J_inv << endl;
	Matrix J_inv = J.Inversed();
	//cout << "J_inv = " << endl << J_inv << endl;
	return J_inv * derFF;
}
Matrix QuadElement::CreateMatrixB(double ksi, double eta)
{
	Matrix B(3, 8);
	Matrix W = DerFFGlob(ksi, eta);
	//cout << "W = " << endl << W << endl;
	//cout << "---" << endl;
	B.a[0][0] = B.a[2][1] = W.a[0][0];	//Ô1 x
	B.a[1][1] = B.a[2][0] = W.a[1][0];	//Ô1 y
	B.a[0][2] = B.a[2][3] = W.a[0][1];	//Ô2 x
	B.a[1][3] = B.a[2][2] = W.a[1][1];	//Ô2 y
	B.a[0][4] = B.a[2][5] = W.a[0][2];	//Ô3 x
	B.a[1][5] = B.a[2][4] = W.a[1][2];	//Ô3 y
	B.a[0][6] = B.a[2][7] = W.a[0][3];	//Ô4 x
	B.a[1][7] = B.a[2][6] = W.a[1][3];	//Ô4 y

	return B;
}
double2 QuadElement::GetCentre()
{
	double2 c, sum;
	sum.x = sum.y = 0;
	for (int i = 0; i < 4; i++)
	{
		sum.x += coord[i].x;
		sum.y += coord[i].y;
	}
	c.x = sum.x / 4; c.y = sum.y / 4;
	return c;
}
Matrix QuadElement::CreateMatrixK()
{
	Matrix A = mat.GetA();
	vector<Matrix> B(4);
	B[0] = CreateMatrixB(-a, -a); B[1] = CreateMatrixB(a, -a); B[2] = CreateMatrixB(a, a); B[3] = CreateMatrixB(-a, a);
	Matrix K(8, 8);
	for (int i = 0; i < 4; i++)
	{
		K = K + B[i].T() * A * B[i];
	}
	return K;
}
void QuadElement::CreateElement(vector<double3>& pt_list, vector<int3>& hex_list, vector<int3>& nums, int i)
{
	vector<QuadElement> elems2(pt_list.size());
	vector<double3> coord_of_element_nodes(4);
	for (int j = 0; j < 4; j++)
	{
		coord_of_element_nodes[j] = pt_list[hex_list[i].n[j]];
		nums[i].n[j] = hex_list[i].n[j];
	}
	*this = QuadElement(coord_of_element_nodes, mat);
}

void QuadElement::Print()
{
	for (int i = 0; i < 4; i++)
	{
		cout << coord[i].x << " " << coord[i].y << endl;
	}
	cout << endl;
}