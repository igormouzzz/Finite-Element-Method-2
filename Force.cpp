#include "Force.h"

Force::Force() {}
Force::Force(double3_ f, vector<int> nums)
{
	F = f; nodes = nums;
}
Force::Force(double3_ f, int3_ nums)
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
double3_ Force::GetF() { return F; }
vector<int> Force::GetNumbersOfNodes() { return nodes; }