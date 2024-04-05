#include "Header.h"
#include "Element.h"
#include "Force.h"
#include "Restraint.h"
#include "Local.h"

int Example()
{
	srand(time(NULL));
	vector<double> a1 = { 32,5,0,0,0 };
	vector<double> a2 = { 5,6,0,0,0 };
	vector<double> a3 = { 0,0,29,15,0 };
	vector<double> a4 = { 0,0,15,10,0 };
	vector<double> a5 = { 0,0,0,0,1 };
	vector<row> M = { a1,a2,a3,a4,a5 };
	vector<double> b = { 1,1,0,0,9 };
	Matrix A(M);
	cout << A << endl;
	vector<double> x1 = A.Gauss(b);
	for (int i = 0; i < x1.size(); i++) cout << x1[i] << " ";
	cout << endl;
	vector<double> x2 = A.CG3(b, pow(10, -4));
	for (int i = 0; i < x2.size(); i++) cout << x2[i] << " ";
	cout << endl;

	return 0;
}

int Task()
{
	struct timespec ts1, ts2;
	char timestamp1[100], timestamp2[100];
	int t = 0;
	double sec = 0;
	timespec_get(&ts1, TIME_UTC);
	vector<double3> list_of_nodes_with_coords;
	vector<vc> list_elements_with_nodes;
	double E, Nu;
	Force F;
	vector<Restraint> R;
	int number_of_nodes_of_elem = 0;
	Read(list_of_nodes_with_coords, list_elements_with_nodes, F, R, E, Nu, number_of_nodes_of_elem);
	cout << endl << endl;

	

	Material mat(E, Nu);

	int element_count = int(list_elements_with_nodes.size());
	int node_count = int(list_of_nodes_with_coords.size());

#if COUNT_OF_NODES == 3
	vector<TriangularElement> elements(element_count);
#else
	vector<QuadElement> elements(element_count);
#endif
	vector<Matrix> matrices(element_count);
	vector<vc> nums(element_count);

	for (int i = 0; i < elements.size(); i++)
	{
		elements[i].SetMaterial(mat);
		elements[i].CreateElement(list_of_nodes_with_coords, list_elements_with_nodes, nums, i);
	}

	for (int i = 0; i < elements.size(); i++)
	{
		matrices[i] = elements[i].CreateMatrixK();
	}
	
	Matrix K = MadeGlobalStiffnessMatrix(element_count, node_count, matrices, nums);
	cout << K.GetN() << endl;

	for (int k = 0; k < R.size(); k++)
	{
		Restraint::ApplyRestraints(K, R[k]);
	}

	vector<double> b(K.GetN());
	for (int i = 0; i < F.GetNumbersOfNodes().size(); i++)
	{
		//cout << F.GetNumbersOfNodes()[i] << " ";
		b[2 * F.GetNumbersOfNodes()[i] - 2] = F.GetF().x;
		b[2 * F.GetNumbersOfNodes()[i] - 1] = F.GetF().y;
	}
	
	vector<double> X = K.CG3(b, 0.0001);
	/*for (int i = 0; i < X.size(); i++)
	{
		cout << X[i] << " ";
	}
	cout << endl;*/

	timespec_get(&ts2, TIME_UTC);
	strftime(timestamp1, 100, "%Y-%m-%d %H:%M:%S:", gmtime(&ts1.tv_sec));
	strftime(timestamp2, 100, "%Y-%m-%d %H:%M:%S:", gmtime(&ts2.tv_sec));
	t = ts2.tv_nsec - ts1.tv_nsec;
	cout << timestamp1 << ts1.tv_nsec << endl;
	cout << timestamp2 << ts2.tv_nsec << endl;
	sec = (double(ts2.tv_sec) + double(ts2.tv_nsec) / 1000000000) - (double(ts1.tv_sec) + double(ts1.tv_nsec) / 1000000000);
	cout << sec << endl;

	ofstream x("X.txt");
	for (int i = 0; i < X.size(); i++)
	{
		x << X[i] << " ";
	}
	x << endl;

	vector<Strains> Epsilon(element_count);
	vector<Stresses> Sigma(element_count);

	for (int i = 0; i < elements.size(); i++)
	{
		Epsilon[i] = Element::FindStrains(X, list_of_nodes_with_coords, list_elements_with_nodes, elements[i], i);
		Sigma[i] = Element::FindStresses(X, list_of_nodes_with_coords, list_elements_with_nodes, elements[i], i);
	}
	ofstream f3("strains.txt");
	ofstream f4("stresses.txt");
	//cout << Epsilon.size() << endl;
	for (int i = 0; i < Epsilon.size(); i++)
	{
		f3 << Epsilon[i];
		f4 << Sigma[i];
	}

	ofstream f5("nodes.txt");
	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
#if COUNT_OF_NODES == 3
		f5 << list_elements_with_nodes[i].n[0] << " " << list_elements_with_nodes[i].n[1] << " " << list_elements_with_nodes[i].n[2] << endl;
#else
		f5 << list_elements_with_nodes[i].n[0] << " " << list_elements_with_nodes[i].n[1] << " " << list_elements_with_nodes[i].n[2] << " " << list_elements_with_nodes[i].n[3] << endl;
#endif
	}
	f5.close();
	ofstream f6("coord.txt");
	for (int i = 0; i < list_of_nodes_with_coords.size(); i++)
	{
		f6 << list_of_nodes_with_coords[i].x << " " << list_of_nodes_with_coords[i].y << endl;
	}
	f6.close();

	ofstream f7("elements.txt");
	for (int i = 0; i < elements.size(); i++)
	{
		elements[i].PrintToFile(f7);
	}
	f7.close();
	
	return 0;
}

int Task2()
{
	struct timespec ts1, ts2;
	char timestamp1[100], timestamp2[100];
	int t = 0;
	double sec = 0;
	timespec_get(&ts1, TIME_UTC);
	vector<double3> list_of_nodes_with_coords;
	vector<vc> list_elements_with_nodes;
	double E, Nu;
	Force F;
	vector<Restraint> R;
	int number_of_nodes_of_elem = 0;
	Read(list_of_nodes_with_coords, list_elements_with_nodes, F, R, E, Nu, number_of_nodes_of_elem);
	cout << endl << endl;



	Material mat(E, Nu);

	int element_count = int(list_elements_with_nodes.size());
	int node_count = int(list_of_nodes_with_coords.size());

#if COUNT_OF_NODES == 3
	vector<TriangularElement> elements(element_count);
#else
	vector<QuadElement> elements(element_count);
#endif
	vector<Matrix> matrices(element_count);
	vector<vc> nums(element_count);

	for (int i = 0; i < elements.size(); i++)
	{
		elements[i].SetMaterial(mat);
		elements[i].CreateElement(list_of_nodes_with_coords, list_elements_with_nodes, nums, i);
	}

	for (int i = 0; i < elements.size(); i++)
	{
		matrices[i] = elements[i].CreateMatrixK();
	}

	/*
	for (int i = 0; i < F.GetNumbersOfNodes().size(); i++)
	{
		cout << F.GetNumbersOfNodes()[i] << " ";
	}
	cout << endl;
	*/
	vector<int> F_nodes = F.GetNumbersOfNodes();
	vector<double> ElemForces(2 * COUNT_OF_NODES * elements.size());
	for (int i = 0; i < list_elements_with_nodes.size();)
	{
		for (int j = 0; j < COUNT_OF_NODES; j++)
		{
			for (int k = 0; k < F_nodes.size(); k++)
			{
				if (list_elements_with_nodes[i].n[j] == F_nodes[k] - 1)
				{
					ElemForces[i + 2 * j] = F.GetF().x;
					ElemForces[i + 2 * j + 1] = F.GetF().y;
				}
			}
		}
		i += 2 * COUNT_OF_NODES;
	}

#if COUNT_OF_NODES == 3
	DivisionToLocalsTri Local(ElemForces, matrices);
#else
	DivisionToLocalsFour Local(ElemForces, matrices);
#endif

	//Local.PrintVector();
	
	vector<vector<double>> X_local = Local.SolveCG();

	for (int i = 0; i < X_local.size(); i++)
	{
		for (int j = 0; j < 2 * COUNT_OF_NODES; j++)
		{
			cout << X_local[i][j] << " ";
		}
		cout << endl;
	}
	
	timespec_get(&ts2, TIME_UTC);
	strftime(timestamp1, 100, "%Y-%m-%d %H:%M:%S:", gmtime(&ts1.tv_sec));
	strftime(timestamp2, 100, "%Y-%m-%d %H:%M:%S:", gmtime(&ts2.tv_sec));
	t = ts2.tv_nsec - ts1.tv_nsec;
	cout << timestamp1 << ts1.tv_nsec << endl;
	cout << timestamp2 << ts2.tv_nsec << endl;
	sec = (double(ts2.tv_sec) + double(ts2.tv_nsec) / 1000000000) - (double(ts1.tv_sec) + double(ts1.tv_nsec) / 1000000000);
	cout << sec << endl;

	/*
	vector<Strains> Epsilon(element_count);
	vector<Stresses> Sigma(element_count); 
	
	for (int i = 0; i < elements.size(); i++)
	{
		Epsilon[i] = Element::FindStrains(X, list_of_nodes_with_coords, list_elements_with_nodes, elements[i], i);
		Sigma[i] = Element::FindStresses(X, list_of_nodes_with_coords, list_elements_with_nodes, elements[i], i);
	}
	ofstream f3("strains.txt");
	ofstream f4("stresses.txt");
	//cout << Epsilon.size() << endl;
	for (int i = 0; i < Epsilon.size(); i++)
	{
		f3 << Epsilon[i];
		f4 << Sigma[i];
	}

	ofstream f5("nodes.txt");
	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
#if COUNT_OF_NODES == 3
		f5 << list_elements_with_nodes[i].n[0] << " " << list_elements_with_nodes[i].n[1] << " " << list_elements_with_nodes[i].n[2] << endl;
#else
		f5 << list_elements_with_nodes[i].n[0] << " " << list_elements_with_nodes[i].n[1] << " " << list_elements_with_nodes[i].n[2] << " " << list_elements_with_nodes[i].n[3] << endl;
#endif
	}
	f5.close();
	ofstream f6("coord.txt");
	for (int i = 0; i < list_of_nodes_with_coords.size(); i++)
	{
		f6 << list_of_nodes_with_coords[i].x << " " << list_of_nodes_with_coords[i].y << endl;
	}
	f6.close();

	ofstream f7("elements.txt");
	for (int i = 0; i < elements.size(); i++)
	{
		elements[i].PrintToFile(f7);
	}
	f7.close();
	*/
	return 0;
}