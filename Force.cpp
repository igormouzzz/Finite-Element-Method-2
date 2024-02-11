#include "Force.h"

Force::Force() {}
Force::Force(double3 f, vector<int> nums)
{
	F = f; nodes = nums;
}
Force::Force(double3 f, int3 nums)
{
	F = f;
	nodes.resize(3);
	for (int j = 0; j < 3; j++) nodes[j] = nums.n[j];
}
Force& Force::operator=(Force b)
{
	F = b.F;
	nodes = b.nodes;
	return *this;
}
double3 Force::GetF() { return F; }
vector<int> Force::GetNumbersOfNodes() { return nodes; }