#include "Local.h"

DivisionToLocalsTri::DivisionToLocalsTri(vector<double>& b, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies)
{
	n_adjelem = n_adj;
	list_elements_with_nodes = list_elements_with_nodes2;
	M = matricies;
	v.resize(b.size());
	for (int i = 0; i < b.size(); i++)
	{
		v[i].resize(6);
	}
	MakeLocalVectors(b, v);
}

void DivisionToLocalsTri::Multiply(vector<vector_loc>& x_loc, vector<vector_loc>& b_loc)
{
	for (int i = 0; i < M.size(); i++)
	{
		b_loc[i] = M[i] * x_loc[i];
	}
}

void DivisionToLocalsTri::MakeGlobalVector(vector<vector_loc>& x_loc, vector<double>& X)
{
	for (int j = 0; j < x_loc.size(); j++)		//b_loc
	{
		for (int k = 0; k < 3; k++)
		{
			X[2*list_elements_with_nodes[j].n[k]] += x_loc[j][2*k];
			X[2*list_elements_with_nodes[j].n[k]+1] += x_loc[j][2*k+1];
		}
	}
}

void DivisionToLocalsTri::MakeLocalVectors(vector<double>& b, vector<vector_loc>& b_loc)
{
	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			b_loc[i][2*j] = b[2 * list_elements_with_nodes[i].n[j]];
			b_loc[i][2*j+1] = b[2 * list_elements_with_nodes[i].n[j]+1];
		}
	}
}

void DivisionToLocalsTri::PrintVector()
{
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

ostream& operator<<(ostream& cout, DivisionToLocals& b)
{
	b.PrintVector();
	return cout;
}

double DivisionToLocalsTri::norm_square_of_locals(vector<vector_loc>& b)
{
	double norm = 0;
	vector<bool> flag(n_adjelem.size(), false);
	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < COUNT_OF_NODES; j++)
		{
			if (!flag[list_elements_with_nodes[i].n[j]])	//если не посчитали этот узел
			{
				norm =  norm + b[i][2 * j] * b[i][2 * j] + b[i][2 * j + 1] * b[i][2 * j + 1];
				flag[list_elements_with_nodes[i].n[j]] = true;
			}
			else											//иначе пропускаем
			{
				continue;
			}
		}
	}
	return norm;
}

double DivisionToLocalsTri::scalar_product_of_locals(vector<vector_loc>& a, vector<vector_loc>& b)
{
	double sp = 0;
	vector<bool> flag(n_adjelem.size(), false);
	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < COUNT_OF_NODES; j++)
		{
			if (!flag[list_elements_with_nodes[i].n[j]])	//если не посчитали этот узел
			{
				sp = sp + a[i][2 * j] * b[i][2 * j] + a[i][2 * j + 1] * b[i][2 * j + 1];
				flag[list_elements_with_nodes[i].n[j]] = true;
			}
			else											//иначе пропускаем
			{
				continue;
			}
		}
	}
	return sp;
}

vector<double> DivisionToLocalsTri::CG4(vector<double>& b)
{
	int SIZE = M.size();
	vector<double> matr(36 * SIZE);
	vector<double> xx_loc(6 * SIZE);
	vector<double> bb_loc(6 * SIZE);

	for (int k = 0; k < SIZE; k++)
	{
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				matr[36 * k + 6 * i + j] = M[k].a[i][j];
			}
		}
	}

	double b_norm = norm_square(b);
	//cout << b_norm << endl;
	if (b_norm == 0) return vector<double>(b.size());
	int n = b.size();
	vector<double> x(n, 0.0);  // Initial guess for the solution
	vector<double> r = b;      // Residual vector
	vector<double> p = b;      // Search direction vector

	

	unsigned int iteration = 0;
	unsigned int max_iter = 1000;
	
	vector<vector_loc> p_loc(list_elements_with_nodes.size());
	vector<vector_loc> Ap_loc(list_elements_with_nodes.size());
	for (int i = 0; i < p_loc.size(); i++) { p_loc[i].resize(6); Ap_loc[i].resize(6); }
	vector<double> Ap(2*n_adjelem.size());
	double alpha, beta;

	sycl::queue q;

	sycl::buffer <double, 1> dM(matr.data(), 36 * SIZE);
	sycl::buffer <double, 1> dx_loc(xx_loc.data(), 6 * SIZE);
	sycl::buffer <double, 1> db_loc(bb_loc.data(), 6 * SIZE);

	do
	{
		//cout << iteration << endl;
		MakeLocalVectors(p, p_loc);
		Multiply(p_loc, Ap_loc);		
		MakeGlobalVector(Ap_loc, Ap);
		alpha = scalar_product(r, r) / scalar_product(p, Ap);
		//cout << "alpha = " << alpha << endl;
		for (size_t i = 0; i < n; ++i)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}
		beta = scalar_product(r, r) / scalar_product(p, p);
		for (size_t i = 0; i < n; ++i)
		{
			p[i] = r[i] + beta * p[i];
		}
		iteration++;

		fill(Ap.begin(), Ap.end(), 0);

		//cout << norm_square(r) / b_norm << endl;
	} while (norm_square(r) / b_norm > eps * eps && iteration < max_iter);

	cout << "iterations: " << iteration << endl;
	//cout << "tolerance: " << eps << endl;

	return x;
}