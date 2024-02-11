#include "Restraint.h"

Restraint::Restraint() {}
Restraint::Restraint(double3 f, vector<int> nums, vector<int> flags)
{
	F = f; nodes = nums; flag = flags;
}
Restraint::Restraint(double3 f, int3 nums, vector<int> flags)
{
	F = f;
	nodes.resize(3);
	for (int j = 0; j < 3; j++) nodes[j] = nums.n[j];
	flag = flags;
}
Restraint& Restraint::operator=(Restraint b)
{
	F = b.F;
	nodes = b.nodes;
	flag = b.flag;
	return *this;
}
double3 Restraint::GetF() { return F; }
vector<int> Restraint::GetNumbersOfNodes() { return nodes; }
vector<int> Restraint::GetFlag() { return flag; }

void Restraint::ApplyRestraints(Matrix& K, Restraint R)
{
	vector<int> nodes = R.GetNumbersOfNodes();
	for (int i = 0; i < nodes.size(); i++)
	{
		for (int j = 0; j < K.GetN(); j++)
		{
			if (R.GetFlag()[0])
			{
				if (2 * nodes[i] - 2 == j)
				{
					K.Set(2 * nodes[i] - 2, j, 1.0);
				}
				else
				{

					K.Set(2 * nodes[i] - 2, j, 0.0); K.Set(j, 2 * nodes[i] - 2, 0.0);
				}
			}
			if (R.GetFlag()[1])
			{
				if (2 * nodes[i] - 1 == j)
				{
					K.Set(2 * nodes[i] - 1, j, 1.0);
				}
				else
				{
					K.Set(2 * nodes[i] - 1, j, 0.0); K.Set(j, 2 * nodes[i] - 1, 0.0);
				}
			}
		}
	}
}