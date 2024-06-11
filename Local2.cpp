#include "Local2.h"

DivisionToLocalsTri2::DivisionToLocalsTri2(vector<double>& b, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies)
{
	n_adjelem = n_adj;
	list_elements_with_nodes = list_elements_with_nodes2;
	//matr.resize(matricies.size());
	int size = list_elements_with_nodes.size();
	matr.resize(size);
	/*for (int k = 0; k < size; k++)
	{
		matr[k].resize(6);
		for (int i = 0; i < 6; i++)
		{
			matr[k][i].resize(6);
		}
	}*/
	for (int k = 0; k < matricies.size(); ++k)
	{
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				matr[k][i][j] = matricies[k].a[i][j];
			}
		}
	}
	v.resize(6*list_elements_with_nodes.size());
	MakeLocalVectors(b, v);
}

void DivisionToLocalsTri2::MakeGlobalVector(vector<vector_loc>& x_loc, vector<double>& X)
{
	for (int j = 0; j < x_loc.size(); j++)		//b_loc
	{
		for (int k = 0; k < 3; k++)
		{
			X[2 * list_elements_with_nodes[j].n[k]] += x_loc[j][2*k];
			X[2 * list_elements_with_nodes[j].n[k] + 1] += x_loc[j][2*k+1];
		}
	}
}

void DivisionToLocalsTri2::MakeLocalVectors(vector<double>& b, vector<vector_loc>& b_loc)
{
	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			b_loc[i][2*j] = b[2 * list_elements_with_nodes[i].n[j]];
			b_loc[i][2*j+1] = b[2 * list_elements_with_nodes[i].n[j] + 1];
		}
	}
}

void DivisionToLocalsTri2::PrintVector()
{
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

ostream& operator<<(ostream& cout, DivisionToLocals2& b)
{
	b.PrintVector();
	return cout;
}

/*vector<double> DivisionToLocalsTri2::CG4(vector<double>& b)
{
	struct timespec ts1, ts2;

	int SIZE = list_elements_with_nodes.size();
	double b_norm = norm_square(b);
	if (b_norm == 0) return vector<double>(b.size());
	int n = b.size();
	vector<double> x(n, 0.0);  // Initial guess for the solution
	vector<double> r = b;      // Residual vector
	vector<double> p = b;      // Search direction vector

	unsigned int iteration = 0;
	unsigned int max_iter = 2*n_adjelem.size();

	//vector<double> p_loc(6*SIZE);
	//vector<double> Ap_loc(6 * SIZE);
	vector<vector_loc> p_loc(SIZE);
	vector<vector_loc> Ap_loc(SIZE);
	vector<double> Ap(2 * n_adjelem.size());
	double alpha, beta;

	sycl::queue q(sycl::default_selector_v);

	timespec_get(&ts1, TIME_UTC);
	do
	{
		//cout << iteration << endl;
		MakeLocalVectors(p, p_loc);

		Multiply(q, matr, p_loc, Ap_loc);

		MakeGlobalVector(Ap_loc, Ap);
		alpha = scalar_product(r, r) / scalar_product(p, Ap);
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

	} while (norm_square(r) / b_norm > eps * eps && iteration < max_iter);
	timespec_get(&ts2, TIME_UTC);

	cout << "iterations: " << iteration << endl;
	cout << endl;
	int t = ts2.tv_nsec - ts1.tv_nsec;
	double sec = (double(ts2.tv_sec) + double(ts2.tv_nsec) / 1000000000) - (double(ts1.tv_sec) + double(ts1.tv_nsec) / 1000000000);
	cout << sec << endl;

	return x;
}*/