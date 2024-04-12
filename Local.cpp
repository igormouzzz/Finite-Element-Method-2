#include "Local.h"

DivisionToLocalsTri::DivisionToLocalsTri(vector<double>& b, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies)
{
	n_adjelem = n_adj;
	list_elements_with_nodes = list_elements_with_nodes2;
	M = matricies;

	v = MakeLocalVectors(b);
}

vector<vector_loc> DivisionToLocalsTri::Multiply(vector<vector_loc>& x_loc)
{
	vector<vector_loc> b_loc(M.size());
#pragma omp parallel for
	for (int i = 0; i < M.size(); i++)
	{
		b_loc[i] = M[i] * x_loc[i];
	}
	return b_loc;
}

vector<double> DivisionToLocalsTri::MakeGlobalVector(vector<vector_loc>& x_loc)
{
	vector<double> X(2*n_adjelem.size());
	for (int j = 0; j < x_loc.size(); j++)		//b_loc
	{
		for (int k = 0; k < 3; k++)
		{
			X[2*list_elements_with_nodes[j].n[k]] += x_loc[j][2*k];
			X[2*list_elements_with_nodes[j].n[k]+1] += x_loc[j][2*k+1];
		}
	}
	
	for (int i = 0; i < n_adjelem.size(); i++)
	{
		//X[2*i] /= n_adjelem[i];
		//X[2*i+1] /= n_adjelem[i];
	}
	return X;
}

vector<vector_loc> DivisionToLocalsTri::MakeLocalVectors(vector<double>& b)
{
	vector<vector_loc> b_loc(list_elements_with_nodes.size());
	
	for (int i = 0; i < b_loc.size(); i++) b_loc[i].resize(6);

	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			b_loc[i][2 * j] = b[2 * list_elements_with_nodes[i].n[j]];// *double(n_adjelem[list_elements_with_nodes[i].n[j]]);	//?
			b_loc[i][2 * j + 1] = b[2 * list_elements_with_nodes[i].n[j] + 1];// *double(n_adjelem[list_elements_with_nodes[i].n[j]]);
		}
	}
	return b_loc;
}

vector<vector_loc> DivisionToLocalsTri::SolveGauss()
{
	vector<vector_loc> X_local(v.size());
	for (int i = 0; i < X_local.size(); i++)
	{
		X_local[i] = M[i].Gauss(v[i]);
	}
	return X_local;
}

vector<double> DivisionToLocalsTri::SolveGaussGlobal()
{
	vector<vector_loc> X = SolveGauss();
	return MakeGlobalVector(X);
}

vector<vector_loc> DivisionToLocalsTri::SolveCG()
{
	vector<vector_loc> X_local(v.size());
	for (int i = 0; i < X_local.size(); i++)
	{
		X_local[i] = M[i].CG3(v[i]);
	}
	return X_local;
}

vector<double> DivisionToLocalsTri::SolveCGGlobal()
{
	vector<vector_loc> X = SolveCG();
	return MakeGlobalVector(X);
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
	double b_norm = norm_square(b);
	if (b_norm == 0) return vector<double>(b.size());
	int n = b.size();
	vector<double> x(n, 0.0);  // Initial guess for the solution
	vector<double> r = b;      // Residual vector
	vector<double> p = b;      // Search direction vector

	unsigned int iteration = 0;
	unsigned int max_iter = 1000;
	
	do
	{
		vector<vector_loc> p_loc = MakeLocalVectors(p);
		vector<vector_loc> Ap_loc = Multiply(p_loc);
		vector<double> Ap = MakeGlobalVector(Ap_loc);
		/*for (int i = 0; i < Ap.size(); i++)
		{
			cout << Ap[i] << " ";
		}
		cout << endl;*/
		double alpha = scalar_product(r, r) / scalar_product(p, Ap);
		for (size_t i = 0; i < n; ++i)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}
		double beta = scalar_product(r, r) / scalar_product(p, p);
		for (size_t i = 0; i < n; ++i)
		{
			p[i] = r[i] + beta * p[i];
		}
		iteration++;
	} while (norm_square(r) / b_norm > epsilon * epsilon && iteration < max_iter);

	cout << "iterations: " << iteration << endl;
	//cout << "tolerance: " << epsilon << endl;

	return x;
}