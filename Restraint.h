#include "Header.h"
#include "Matrix.h"

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
};