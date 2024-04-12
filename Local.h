#include "Header.h"
#include "Matrix.h"
using vector_loc = vector<double>;

class DivisionToLocals
{
public:
	virtual void PrintVector() = 0;
	friend ostream& operator<<(ostream& cout, DivisionToLocals& b);
	friend class Restraint;
};

class DivisionToLocalsTri : public DivisionToLocals
{
protected:
	vector<vector_loc> v;	//M*v
	vector<int> n_adjelem;
	vector<vc> list_elements_with_nodes;
	vector<Matrix> M;		// -> x_loc		M*x_loc
	const int size = 6;
public:
	DivisionToLocalsTri(vector<double>& W, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies);
	vector<vector_loc> Multiply(vector<vector_loc>& x_loc);
	vector<double> MakeGlobalVector(vector<vector_loc>& x_loc);
	vector<vector_loc> MakeLocalVectors(vector<double>& b);
	vector<vector_loc> GetV() { return v; }
	vector<Matrix> GetMatrices() { return M; }
	void PrintVector();

	vector<vector_loc> SolveGauss();
	vector<double> SolveGaussGlobal();

	vector<vector_loc> SolveCG();
	vector<double> SolveCGGlobal();

	vector<double> CG4(vector<double>& b);

	vector<double> Test(vector<double>& b)
	{
		vector<vector_loc> b_loc = MakeLocalVectors(b);
		vector<double> b2 = MakeGlobalVector(b_loc);
		return b;
	}
	
	Matrix GetMatrixForIndex(int j) { return M[j]; }

	double norm_square_of_locals(vector<vector_loc>& b);
	double scalar_product_of_locals(vector<vector_loc>& a, vector<vector_loc>& b);
	friend class Restraint;
};