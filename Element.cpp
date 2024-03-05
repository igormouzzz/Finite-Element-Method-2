#include "Element.h"

Strains Element::FindStrains(const vector<double>& X, const vector<double3> list_of_nodes_with_coords, const vector<vc>& list_elements_with_nodes, Element& elem, int index)
{
	Matrix B = elem.CreateMatrixB();
	vector<double> u(2 * COUNT_OF_NODES);

	for (int i = 0; i < COUNT_OF_NODES; i++)
	{
		u[2 * i] = X[2 * list_elements_with_nodes[index].n[i]];
		u[2 * i + 1] = X[2 * list_elements_with_nodes[index].n[i] + 1];
	}

	return Strains(B * u);
}
Stresses Element::FindStresses(const vector<double>& X, const vector<double3> list_of_nodes_with_coords, const vector<vc>& list_elements_with_nodes, Element& elem, int index)
{
	Matrix A = elem.mat.GetA();

	vector<double> strains = FindStrains(X, list_of_nodes_with_coords, list_elements_with_nodes, elem, index).GetStrains();

	return Stresses(A * strains);
}

ostream& operator<<(ostream& cout, Element& b)
{
	b.Print();
	return cout;
}