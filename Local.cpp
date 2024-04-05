#include "Local.h"

DivisionToLocalsTri::DivisionToLocalsTri(vector<double> W, vector<Matrix> matricies)
{
	int n = W.size() / size;
	v.resize(n);
	for (int i = 0; i < n; i++)
	{
		v[i].resize(size);
		for (int j = 0; j < size; j++)
		{
			v[i][j] = W[6 * i + j];
		}
	}
	M = matricies;
}

DivisionToLocalsTri::DivisionToLocalsTri(DivisionToLocalsTri& b)
{
	int n = b.v.size();
	v.resize(n);
	for (int i = 0; i < n; i++)
	{
		v[i].resize(size);
		for (int j = 0; j < size; j++)
		{
			v[i][j] = b.v[i][j];
		}
	}
	M = b.M;
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

vector<vector_loc> DivisionToLocalsTri::SolveCG()
{
	vector<vector_loc> X_local(v.size());
	for (int i = 0; i < X_local.size(); i++)
	{
		X_local[i] = M[i].CG3(v[i], 0.01);
	}
	return X_local;
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