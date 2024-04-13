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
	virtual void CreateElement(vector<double3_>& pt_list, vector<vc>& hex_list, vector<vc>& nums, int i) = 0;
	static Strains FindStrains(const vector<double>& X, const vector<double3_> list_of_nodes_with_coords, const vector<vc>& list_elements_with_nodes, Element& elem, int index);
	static Stresses FindStresses(const vector<double>& X, const vector<double3_> list_of_nodes_with_coords, const vector<vc>& list_elements_with_nodes, Element& elem, int index);
	virtual void Print() = 0;
	virtual void PrintToFile(ofstream& f) = 0;
	friend ostream& operator<<(ostream& cout, const Element& b);
};

class TriangularElement : public Element
{
protected:
	double2_ coord[3];
	Matrix DerFF(double x, double y);
	Matrix CreateMatrixB(double x = 0, double y = 0);
public:
	TriangularElement();
	TriangularElement(vector<double3_> a, Material m);
	TriangularElement(vector<double2_> a, Material m);
	TriangularElement operator=(TriangularElement b);
	void SetMaterial(Material m) { mat = m; }
	Matrix CreateMatrixK();
	void CreateElement(vector<double3_>& pt_list, vector<vc>& hex_list, vector<vc>& nums, int i);
	double2_ GetCoord(int i) { return coord[i]; }
	void Print();
	void PrintToFile(ofstream& f);
};


class QuadElement : public Element
{
protected:
	double2_ coord[4];
	const double a = 0.5773502691;
	const double H = 1.0;
	double* FF(double ksi, double eta);
	Matrix DerFF(double ksi, double eta);
	Matrix Jac(double ksi, double eta);
	Matrix DerFFGlob(double ksi, double eta);
	Matrix CreateMatrixB(double ksi, double eta);
	double2_ GetCentre();
public:
	QuadElement();
	QuadElement(vector<double3_> a, Material m);
	QuadElement(vector<double2_> a, Material m);
	QuadElement operator=(QuadElement b);
	void SetMaterial(Material m) { mat = m; }
	Matrix CreateMatrixK();
	void CreateElement(vector<double3_>& pt_list, vector<vc>& hex_list, vector<vc>& nums, int i);
	void Print();
	void PrintToFile(ofstream& f);
};