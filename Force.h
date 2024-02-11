#include "Header.h"

class Force
{
private:
	double3 F;
	vector<int> nodes;
public:
	Force();
	Force(double3 f, vector<int> nums);
	Force(double3 f, int3 nums);
	Force& operator=(Force b);
	double3 GetF();
	vector<int> GetNumbersOfNodes();
};