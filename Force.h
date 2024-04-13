#include "Header.h"

class Force
{
private:
	double3_ F;
	vector<int> nodes;
public:
	Force();
	Force(double3_ f, vector<int> nums);
	Force(double3_ f, int3_ nums);
	Force& operator=(Force b);
	double3_ GetF();
	void SetF(double x, double y, double z) { F.x = x; F.y = y; F.z = z; }
	vector<int> GetNumbersOfNodes();
};