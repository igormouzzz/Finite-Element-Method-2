#include "Header.h"
#include "Matrix.h"
using vector_loc = vector<double>;

class DivisionToLocals2
{
public:
	virtual void PrintVector() = 0;
	friend ostream& operator<<(ostream& cout, DivisionToLocals& b);
	friend class Restraint;
};

class DivisionToLocalsTri2 : public DivisionToLocals2
{
protected:
	vector<double> v;	//M*v
	vector<int> n_adjelem;
	vector<vc> list_elements_with_nodes;
	vector<double> matr;		// -> x_loc		M*x_loc
	const int size = 6;
public:
	DivisionToLocalsTri2(vector<double>& W, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies);
	void MakeGlobalVector(vector<double>& x_loc, vector<double>& b);
	void MakeLocalVectors(vector<double>& b, vector<double>& b_loc);
	void PrintVector();

	vector<double> CG4(vector<double>& b);

	friend class Restraint;
};