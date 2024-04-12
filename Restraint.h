//#include "Header.h"
//#include "Matrix.h"
#include "Local.h"

class Restraint
{
private:
	double3 F;
	vector<int> nodes;
	vector<int> flag;
public:
	Restraint();
	Restraint(double3 f, vector<int> nums, vector<int> flags);
	Restraint(double3 f, vc nums, vector<int> flags);
	Restraint& operator=(Restraint b);
	double3 GetF();
	vector<int> GetNumbersOfNodes();
	vector<int> GetFlag();
	static void ApplyRestraints(Matrix& K, Restraint R);
	friend class Local;
	static void ApplyRestraintsLocal(DivisionToLocalsTri& L, Restraint R, vector<vc>& list_elements_with_nodes, vector<vector<int>>& list_nodes_with_elem_nums);
	static void ApplyRestraintsLocal2(DivisionToLocalsTri& L, Restraint R, vector<vc>& list_elements_with_nodes, vector<vector<int>>& list_nodes_with_elem_nums);
};