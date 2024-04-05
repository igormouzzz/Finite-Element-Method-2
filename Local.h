#include "Header.h"
#include "Matrix.h"
using vector_loc = vector<double>;

class DivisionToLocals
{
public:
	virtual void PrintVector() = 0;
	friend ostream& operator<<(ostream& cout, DivisionToLocals& b);
};

class DivisionToLocalsTri : public DivisionToLocals
{
protected:
	vector<vector_loc> v;
	vector<Matrix> M;
	const int size = 6;
public:
	DivisionToLocalsTri(vector<double> W, vector<Matrix> matricies);
	DivisionToLocalsTri(DivisionToLocalsTri& b);
	vector<vector_loc> GetV() { return v; }
	void PrintVector();
	vector<vector_loc> SolveGauss();
	vector<vector_loc> SolveCG();
};
