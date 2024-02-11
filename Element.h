#include "Header.h"
#include "Material.h"
#include "Strains.h"
#include "Stresses.h"

class Element
{
protected:
	Material mat;
	virtual Matrix DerFF(double x = 0, double y = 0) = 0;
	virtual Matrix CreateMatrixB(double x = 0, double y = 0) = 0;
public:
	virtual void SetMaterial(Material m) = 0;
	virtual Matrix CreateMatrixK() = 0;
	virtual void CreateElement(vector<double3>& pt_list, vector<int3>& hex_list, vector<int3>& nums, int i) = 0;
	static Strains FindStrains(const vector<double>& X, const vector<double3> list_of_nodes_with_coords, const vector<int3>& list_elements_with_nodes, Element& elem, int index);
	static Stresses FindStresses(const vector<double>& X, const vector<double3> list_of_nodes_with_coords, const vector<int3>& list_elements_with_nodes, Element& elem, int index);
	virtual void Print() = 0;
	friend ostream& operator<<(ostream& cout, const Element& b);
};

class TriangularElement : public Element
{
protected:
	double2 coord[3];
	Matrix DerFF(double x, double y);
	Matrix CreateMatrixB(double x = 0, double y = 0);
public:
	TriangularElement();
	TriangularElement(vector<double3> a, Material m);
	TriangularElement(vector<double2> a, Material m);
	TriangularElement operator=(TriangularElement b);
	void SetMaterial(Material m) { mat = m; }
	Matrix CreateMatrixK();
	void CreateElement(vector<double3>& pt_list, vector<int3>& hex_list, vector<int3>& nums, int i);
	void Print();
};


class QuadElement : public Element
{
protected:
	double2 coord[4];
	const double a = 0.5773502691;
	const double H = 1.0;
	double* FF(double ksi, double eta);
	Matrix DerFF(double ksi, double eta);
	Matrix Jac(double ksi, double eta);
	Matrix DerFFGlob(double ksi, double eta);
	Matrix CreateMatrixB(double ksi, double eta);
	double2 GetCentre();
public:
	QuadElement();
	QuadElement(vector<double3> a, Material m);
	QuadElement(vector<double2> a, Material m);
	QuadElement operator=(QuadElement b);
	void SetMaterial(Material m) { mat = m; }
	Matrix CreateMatrixK();
	void CreateElement(vector<double3>& pt_list, vector<int3>& hex_list, vector<int3>& nums, int i);
	void Print();
};