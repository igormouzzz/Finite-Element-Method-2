#include "Element.h"

TriangularElement::TriangularElement() {}
TriangularElement::TriangularElement(vector<double3> a, Material m)
{
	for (int i = 0; i < 3; i++)
	{
		coord[i].x = a[i].x;
		coord[i].y = a[i].y;
	}
	mat = m;
}
TriangularElement::TriangularElement(vector<double2> a, Material m)
{
	for (int i = 0; i < 3; i++)
	{
		coord[i].x = a[i].x;
		coord[i].y = a[i].y;
	}
	mat = m;
}
TriangularElement TriangularElement::operator=(TriangularElement b)
{
	for (int i = 0; i < 3; i++)
	{
		this->coord[i] = b.coord[i];
	}
	return *this;
}
Matrix TriangularElement::DerFF(double x, double y)
{
	Matrix W(2, 3);
	W.a[0][0] = coord[1].y - coord[2].y; 
	W.a[1][0] = coord[2].x - coord[1].x;
	W.a[0][1] = coord[2].y - coord[0].y;
	W.a[1][1] = coord[0].x - coord[2].x;
	W.a[0][2] = coord[0].y - coord[1].y;
	W.a[1][2] = coord[1].x - coord[0].x;
	return W;
}
Matrix TriangularElement::CreateMatrixB(double x, double y)
{
	double Delta;
	Matrix B(3, 6);
	double C;
	Matrix W = DerFF(x,y);
	B.a[0][0] = B.a[2][1] = W.a[0][0];
	B.a[0][2] = B.a[2][3] = W.a[0][1];
	B.a[0][4] = B.a[2][5] = W.a[0][2];
	B.a[1][1] = B.a[2][0] = W.a[1][0];
	B.a[1][3] = B.a[2][2] = W.a[1][1];
	B.a[1][5] = B.a[2][4] = W.a[1][2];
	Delta = 0.5 * det3(1, coord[0].x, coord[0].y, 1, coord[1].x, coord[1].y, 1, coord[2].x, coord[2].y);
	C = 0.5 / Delta;
	B = B * C;
	return B;
}
Matrix TriangularElement::CreateMatrixK()
{
	double Delta = 0.5 * det3(1, coord[0].x, coord[0].y, 1, coord[1].x, coord[1].y, 1, coord[2].x, coord[2].y);
	Matrix A = mat.GetA();
	Matrix B = CreateMatrixB();
	return B.T() * A * B * Delta;
}
void TriangularElement::CreateElement(vector<double3>& pt_list, vector<int3>& hex_list, vector<int3>& nums, int i)
{
	vector<TriangularElement> elems2(pt_list.size());
	vector<double3> coord_of_element_nodes(3);
	for (int j = 0; j < 3; j++)
	{
		coord_of_element_nodes[j] = pt_list[hex_list[i].n[j]];
		nums[i].n[j] = hex_list[i].n[j];
	}
	*this = TriangularElement(coord_of_element_nodes, mat);
}

void TriangularElement::Print()
{
	for (int i = 0; i < 3; i++)
	{
		cout << coord[i].x << " " << coord[i].y << endl;
	}
	cout << endl;
}