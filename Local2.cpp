#include "Local2.h"

DivisionToLocalsTri2::DivisionToLocalsTri2(vector<double>& b, vector<int>& n_adj, vector<vc>& list_elements_with_nodes2, vector<Matrix>& matricies)
{
	n_adjelem = n_adj;
	list_elements_with_nodes = list_elements_with_nodes2;
	matr.resize(36*matricies.size());
	for (int k = 0; k < matricies.size(); ++k)
	{
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				matr[36 * k + 6 * i + j] = matricies[k].a[i][j];
			}
		}
	}
	v.resize(6*list_elements_with_nodes.size());
	MakeLocalVectors(b, v);
}

void DivisionToLocalsTri2::MakeGlobalVector(vector<double>& x_loc, vector<double>& X)
{
	for (int j = 0; j < x_loc.size(); j++)		//b_loc
	{
		for (int k = 0; k < 3; k++)
		{
			X[2 * list_elements_with_nodes[j].n[k]] += x_loc[6*j + 2*k];
			X[2 * list_elements_with_nodes[j].n[k] + 1] += x_loc[6*j + 2*k+1];
		}
	}
}

void DivisionToLocalsTri2::MakeLocalVectors(vector<double>& b, vector<double>& b_loc)
{
	for (int i = 0; i < list_elements_with_nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			b_loc[6*i + 2*j] = b[2 * list_elements_with_nodes[i].n[j]];
			b_loc[6*i + 2*j+1] = b[2 * list_elements_with_nodes[i].n[j] + 1];
		}
	}
}

void DivisionToLocalsTri2::PrintVector()
{
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cout << v[6*i + j] << " ";
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


vector<double> DivisionToLocalsTri2::CG4(vector<double>& b)
{
	int SIZE = list_elements_with_nodes.size();
	vector<double> xx_loc(6 * SIZE);
	vector<double> bb_loc(6 * SIZE);

	double b_norm = norm_square(b);
	//cout << b_norm << endl;
	if (b_norm == 0) return vector<double>(b.size());
	int n = b.size();
	vector<double> x(n, 0.0);  // Initial guess for the solution
	vector<double> r = b;      // Residual vector
	vector<double> p = b;      // Search direction vector



	unsigned int iteration = 0;
	unsigned int max_iter = 1000;

	vector<double> p_loc(6*SIZE);
	vector<double> Ap_loc(6 * SIZE);
	vector<double> Ap(2 * n_adjelem.size());
	double alpha, beta;

	sycl::queue q;

	sycl::buffer <double, 1> dM(matr.data(), 36 * SIZE);
	sycl::buffer <double, 1> dx_loc(p_loc.data(), 6 * SIZE);
	sycl::buffer <double, 1> db_loc(Ap_loc.data(), 6 * SIZE);

	do
	{
		//cout << iteration << endl;
		MakeLocalVectors(p, p_loc);

		//sycl::buffer <double, 1> dx_loc(xx_loc.data(), 6 * SIZE);
		//sycl::buffer <double, 1> db_loc(bb_loc.data(), 6 * SIZE);

		q.submit([&](sycl::handler& h) {

			sycl::accessor pM(dM, h, sycl::read_only);
			sycl::accessor px_loc(dx_loc, h, sycl::read_only);
			sycl::accessor pb_loc(db_loc, h, sycl::write_only);

			h.parallel_for(SIZE, [=](auto k)
				{
					for (int i = 0; i < 6; i++)
					{
						double val = 0;
						for (int j = 0; j < 6; j++)
						{
							val += pM[36 * k + 6 * i + j] * px_loc[6 * k + j];
						}
						pb_loc[6 * k + i] = val;
					}

				});
			});

		q.wait();

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