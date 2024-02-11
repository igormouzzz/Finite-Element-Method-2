#include "Header.h"
#include "Element.h"
#include "Force.h"
#include "Restraint.h"

int Task()
{
	struct timespec ts1, ts2;
	char timestamp1[100], timestamp2[100];
	int t = 0;
	double sec = 0;
	timespec_get(&ts1, TIME_UTC);
	vector<double3> list_of_nodes_with_coords;
	vector<int3> list_elements_with_nodes;
	double E, Nu;
	Force F; vector<Restraint> R;
	int number_of_nodes_of_elem = 0;
	Read(list_of_nodes_with_coords, list_elements_with_nodes, F, R, E, Nu, number_of_nodes_of_elem);
	cout << endl << endl;

	Material mat(E, Nu);

	int element_count = int(list_elements_with_nodes.size());
	int node_count = int(list_of_nodes_with_coords.size());
	
	vector<TriangularElement> elements(element_count);
	//vector<QuadElement> elements(element_count);
	vector<Matrix> matrices(element_count);
	vector<int3> nums(element_count);
	vector<Strains> Epsilon(element_count);
	vector<Stresses> Sigma(element_count);

	for (int i = 0; i < elements.size(); i++)
	{
		elements[i].SetMaterial(mat);
		elements[i].CreateElement(list_of_nodes_with_coords, list_elements_with_nodes, nums, i);
		//cout << nums[i].n[0] << nums[i].n[1] << nums[i].n[2] << endl;
		//cout << i << ":\t" << endl;
		//cout << elements[i] << endl;
	}

	for (int i = 0; i < elements.size(); i++)
	{
		matrices[i] = elements[i].CreateMatrixK();
		//cout << matrices[i] << endl;
	}

	Matrix K = MadeGlobalStiffnessMatrix(element_count, node_count, matrices, nums);

	for (int k = 0; k < R.size(); k++)
	{
		Restraint::ApplyRestraints(K, R[k]);
	}
	//cout << K << endl;


	vector<double> b(K.GetN());
	for (int i = 0; i < F.GetNumbersOfNodes().size(); i++)
	{
		//cout << F.GetNumbersOfNodes()[i] << " ";
		b[2 * F.GetNumbersOfNodes()[i] - 2] = F.GetF().x;
		b[2 * F.GetNumbersOfNodes()[i] - 1] = F.GetF().y;
	}
	cout << K.GetN() << endl;
	vector<double> X = K.Gauss(b);
	//vector<double> X = K.Iter(b);

	for (int i = 0; i < elements.size(); i++)
	{
		Epsilon[i] = Element::FindStrains(X, list_of_nodes_with_coords, list_elements_with_nodes, elements[i], i);
		Sigma[i] = Element::FindStresses(X, list_of_nodes_with_coords, list_elements_with_nodes, elements[i], i);
	}
	ofstream f3("strains.txt");
	ofstream f4("stresses.txt");
	cout << Epsilon.size() << endl;
	for (int i = 0; i < Epsilon.size(); i++)
	{
		f3 << Epsilon[i];
		f4 << Sigma[i];
	}

	timespec_get(&ts2, TIME_UTC);
	strftime(timestamp1, 100, "%Y-%m-%d %H:%M:%S:", gmtime(&ts1.tv_sec));
	strftime(timestamp2, 100, "%Y-%m-%d %H:%M:%S:", gmtime(&ts2.tv_sec));
	t = ts2.tv_nsec - ts1.tv_nsec;
	cout << timestamp1 << ts1.tv_nsec << endl;
	cout << timestamp2 << ts2.tv_nsec << endl;
	sec = (double(ts2.tv_sec) + double(ts2.tv_nsec) / 1000000000) - (double(ts1.tv_sec) + double(ts1.tv_nsec) / 1000000000);
	cout << sec << endl;

	ofstream f5("nodes.txt");
	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
		f5 << list_elements_with_nodes[i].n[0] << " " << list_elements_with_nodes[i].n[1] << " " << list_elements_with_nodes[i].n[2] << endl;
	}
	f5.close();
	ofstream f6("coord.txt");
	for (int i = 0; i < list_of_nodes_with_coords.size(); i++)
	{
		f6 << list_of_nodes_with_coords[i].x << " " << list_of_nodes_with_coords[i].y << endl;
	}
	f6.close();

	return 0;
}