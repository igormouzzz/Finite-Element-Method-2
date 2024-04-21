#include "Header.h"
#include "Matrix.h"
using vector_loc = vector<double>;

#include <sycl/sycl.hpp>
#if FPGA_HARDWARE || FPGA_EMULATOR || FPGA_SIMULATOR
#include <sycl/ext/intel/fpga_extensions.hpp>
#endif

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
	vector<vector_loc> v;	//M*v
	vector<int> n_adjelem;
	vector<vc> list_elements_with_nodes;
	vector<vector<vector<double>>> matr;
	const int size = 6;
public:
	DivisionToLocalsTri2(vector<double>& W, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies);
	void MakeGlobalVector(vector<vector_loc>& x_loc, vector<double>& b);
	void MakeLocalVectors(vector<double>& b, vector<vector_loc>& b_loc);
	void PrintVector();

	vector<double> CG4(vector<double>& b);

	friend class Restraint;
};