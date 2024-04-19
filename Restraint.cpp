#include "Restraint.h"

Restraint::Restraint() {}
Restraint::Restraint(double3_ f, vector<int> nums, vector<int> flags)
{
	F = f; nodes = nums; flag = flags;
}
Restraint::Restraint(double3_ f, vc nums, vector<int> flags)
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
double3_ Restraint::GetF() { return F; }
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
					K.Set(2 * nodes[i] - 2, j, 1.0);	//out of range
				}
				else
				{
					K.Set(2 * nodes[i] - 2, j, 0.0);	//out of range
					K.Set(j, 2 * nodes[i] - 2, 0.0);	//out of range
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
					K.Set(2 * nodes[i] - 1, j, 0.0);
					K.Set(j, 2 * nodes[i] - 1, 0.0);
				}
			}
		}
	}
}

void Restraint::ApplyRestraintsLocal(DivisionToLocalsTri& L, Restraint R, vector<vc>& list_elements_with_nodes, vector<vector<int>>& list_nodes_with_elem_nums)
{
	vector<int> nodes = R.GetNumbersOfNodes();
	for (int i = 0; i < nodes.size(); i++)		//идём по закреплённым узлам
	{
		for (int j = 0; j < list_nodes_with_elem_nums[nodes[i]-1].size(); j++)	//идём по элементам, которым принадлежат закреп узлы
		{
			for (int k = 0; k < 3; k++)			//идём по локальным номерам узлов каждого рассматриваемого элемента
			{
				if (list_elements_with_nodes[list_nodes_with_elem_nums[nodes[i] - 1][j]].n[k] == nodes[i]-1)	//если номер локального узла совпадает с рассматриваемым закреплённым
				{
					if (R.GetFlag()[0])
					{
						L.M[list_nodes_with_elem_nums[nodes[i] - 1][j]].UnitRowAndColumn(2*k, 2*k, list_nodes_with_elem_nums[nodes[i] - 1].size());			//функция "вычёркивания" строки и столбца (меняем числа на 1 и 0)
					}
					if (R.GetFlag()[1])
					{
						L.M[list_nodes_with_elem_nums[nodes[i] - 1][j]].UnitRowAndColumn(2*k+1, 2*k+1, list_nodes_with_elem_nums[nodes[i] - 1].size());
					}
				}
			}
		}
	}
}

void Restraint::ApplyRestraintsLocal2(DivisionToLocalsTri2& L, Restraint R, vector<vc>& list_elements_with_nodes, vector<vector<int>>& list_nodes_with_elem_nums)
{
	vector<int> nodes = R.GetNumbersOfNodes();
	for (int i = 0; i < nodes.size(); i++)		//идём по закреплённым узлам
	{
		for (int k = 0; k < 3; k++)			//идём по локальным номерам узлов каждого рассматриваемого элемента
		{
			if (list_elements_with_nodes[list_nodes_with_elem_nums[nodes[i] - 1][0]].n[k] == nodes[i] - 1)	//если номер локального узла совпадает с рассматриваемым закреплённым
			{
				if (R.GetFlag()[0])
				{
					//L.M[list_nodes_with_elem_nums[nodes[i] - 1][0]].UnitRowAndColumn(2 * k, 2 * k, 1);			//функция "вычёркивания" строки и столбца (меняем числа на 1 и 0)
					L.matr[36 * list_nodes_with_elem_nums[nodes[i] - 1][0] + 6 * 2 * k + 2 * k] = 1;
				}
				if (R.GetFlag()[1])
				{
					//L.M[list_nodes_with_elem_nums[nodes[i] - 1][0]].UnitRowAndColumn(2 * k + 1, 2 * k + 1, 1);
					L.matr[36 * list_nodes_with_elem_nums[nodes[i] - 1][0] + 6 * (2*k+1) + 2 * k+1] = 1;
				}
			}
		}
		for (int j = 1; j < list_nodes_with_elem_nums[nodes[i] - 1].size(); j++)	//идём по элементам, которым принадлежат закреп узлы
		{
			for (int k = 0; k < 3; k++)			//идём по локальным номерам узлов каждого рассматриваемого элемента
			{
				if (list_elements_with_nodes[list_nodes_with_elem_nums[nodes[i] - 1][j]].n[k] == nodes[i] - 1)	//если номер локального узла совпадает с рассматриваемым закреплённым
				{
					if (R.GetFlag()[0])
					{
						//L.M[list_nodes_with_elem_nums[nodes[i] - 1][j]].ZeroRowAndColumn(2 * k, 2 * k);			//функция "вычёркивания" строки и столбца (меняем числа на 1 и 0)
						L.matr[36 * list_nodes_with_elem_nums[nodes[i] - 1][j] + 6 * 2 * k + 2 * k] = 0;
					}
					if (R.GetFlag()[1])
					{
						//L.M[list_nodes_with_elem_nums[nodes[i] - 1][j]].ZeroRowAndColumn(2 * k + 1, 2 * k + 1);
						L.matr[36 * list_nodes_with_elem_nums[nodes[i] - 1][j] + 6 * (2 * k+1) + 2 * k+1] = 0;
					}
				}
			}
		}
	}
}