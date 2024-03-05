#include "Matrix.h"

ostream& operator<<(ostream& cout, const Matrix& b)
{
	for (int i = 0; i < b.N; i++)
	{
		for (int j = 0; j < b.M; j++)
		{
			cout << b.a[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl << endl;
	return cout;
}

void Matrix::SetZero() { N = M = 0; }

void Matrix::Clean()
{
	a.clear();
	SetZero();
}

size_t Matrix::GetN() { return N; }

void Matrix::Set(int i, int j, double value) { a[i][j] = value; }

Matrix::Matrix() { SetZero(); }

Matrix::Matrix(int N, int M)
{
	this->N = N;
	this->M = M;
	this->a.resize(N);
	for (int i = 0; i < N; i++)
	{
		this->a[i].resize(M);
	}
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			//this->a[i][j] = rand()% 20 - 10;
			this->a[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix& b)
{
	this->N = b.N;
	this->M = b.M;
	this->a.resize(N);
	for (int i = 0; i < N; i++)
	{
		this->a[i].resize(M);
	}
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			this->a[i][j] = b.a[i][j];
		}
	}
}

Matrix& Matrix:: operator=(const Matrix& b)
{
	if (this != &b)
	{
		Clean();
		this->N = b.N;
		this->M = b.M;
		this->a.resize(N);
		for (int i = 0; i < N; i++)
		{
			this->a[i].resize(M);
		}
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				this->a[i][j] = b.a[i][j];
			}
		}
	}
	return *this;
}

Matrix Matrix::operator+(const Matrix& b)
{
	if ((N != b.N) || (M != b.M))
	{
		//cout << "No" << endl;
		return Matrix();
	}
	else
	{
		Matrix S(N, M);
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				S.a[i][j] = a[i][j] + b.a[i][j];
			}
		}

		return S;
	}
}
Matrix Matrix::operator-(const Matrix& b)
{
	if ((N != b.N) || (M != b.M))
	{
		cout << "No" << endl;
		return Matrix();
	}
	else
	{
		Matrix S(N, M);
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				S.a[i][j] = a[i][j] - b.a[i][j];
			}
		}

		return S;
	}
}
Matrix Matrix::operator*(double k)
{
	Matrix S(N, M);
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			S.a[i][j] = k * a[i][j];
		}
	}
	return S;
}
Matrix Matrix::operator*(const Matrix& b)
{
	if (M != b.N)
	{
		//cout << "No" << endl;
		return Matrix();
	}
	else
	{
		Matrix S(N, b.M);
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < b.M; j++)
			{
				S.a[i][j] = a[i][0] * b.a[0][j];
				for (int k = 1; k < M; k++)
				{
					S.a[i][j] += a[i][k] * b.a[k][j];
				}
			}
		}
		return S;
	}
}
vector<double> Matrix::operator*(vector<double>& b)
{
	if (M != b.size())
	{
		cout << "No" << endl;
		throw - 100;
	}
	else
	{
		vector<double> s(N);
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				s[i] += a[i][j] * b[j];
			}
		}
		return s;
	}
}
Matrix Matrix::T()
{
	Matrix S(M, N);
#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			S.a[i][j] = a[j][i];
		}
	}
	return S;
}
Matrix Matrix::Inv2()
{
	Matrix Y(2, 2);
	double detinv = 1 / det2(a[0][0], a[1][0], a[0][1], a[1][1]);
	Y.a[0][0] = a[1][1];
	Y.a[1][1] = a[0][0];
	Y.a[1][0] = -a[1][0];
	Y.a[0][1] = -a[0][1];
	return Y * detinv;
}
Matrix Matrix::Inversed()
{
	Matrix Y(N, 2 * N);
	double temp;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Y.a[i][j] = a[i][j];
		}
		Y.a[i][N + i] = 1;
	}
	for (int i = N - 1; i > 0; i--)
	{
		if (Y.a[i - 1][0] < Y.a[i][0])
		{
			row temp = Y.a[i];
			Y.a[i] = Y.a[i - 1];
			Y.a[i - 1] = temp;
		}
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (j != i)
			{
				temp = Y.a[j][i] / Y.a[i][i];
				for (int k = 0; k < 2 * N; k++)
				{
					Y.a[j][k] -= Y.a[i][k] * temp;
				}
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		temp = Y.a[i][i];
		for (int j = 0; j < 2 * N; j++)
		{
			Y.a[i][j] = Y.a[i][j] / temp;
		}
	}
	Matrix Inv(N, N);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Inv.a[i][j] = Y.a[i][N + j];
		}
	}
	return Inv;
}
void Matrix::Symmetric()
{
	if (N != M)
	{
		cout << "Matrix is not symmetric!" << endl;
		throw - 1;
	}
	else
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < i; j++)
			{
				a[i][j] = a[j][i];
			}
		}
	}
}
double Matrix::SumInRow(int k)
{
	if (k < 0 || k >= N)
	{
		throw - 1;
	}
	double s = 0;
	for (int i = 0; i < M; i++)
	{
		s += a[k][i];
	}
	return s;
}
void Matrix::SumOfRows()
{
	for (int i = 0; i < N; i++)
	{
		cout << i << ": " << SumInRow(i) << endl;
	}
}

void Matrix::Get_matr(vector<row> matr, int n, vector<row> temp_matr, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matr[ki][kj] = matr[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}
double Matrix::Det0(vector<row> matr, int n)
{
	double temp = 0;   //временная переменная для хранения определителя
	int k = 1;      //степень
	if (n < 1) {
		cout << "Не верный размер матрицы!!!" << endl;
		return 0;
	}
	else if (n == 1)
		temp = matr[0][0];
	else if (n == 2)
		temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
	else
	{
		for (int i = 0; i < n; i++) {
			int m = n - 1;
			//double** temp_matr = new double* [m];
			vector<row> temp_matr(m);
			for (int j = 0; j < m; j++)
				temp_matr[j].resize(m);
			Get_matr(matr, n, temp_matr, 0, i);
			temp = temp + k * matr[0][i] * Det0(temp_matr, m);
			k = -k;
		}
	}
	return temp;
}

double Matrix::Det()
{
	if (N != M)
	{
		cout << "Wrong size!" << endl;
		throw - 1;
	}
	return Det0(a, N);
}

vector<double> Matrix::Gauss(vector<double>& b)
{
	for (int i = 0; i < N - 1; i++) {
		// Поиск максимального элемента в столбце
		int maxi = i;
		double max_value = std::abs(a[i][i]);
		for (int j = i + 1; j < N; j++) {
			if (std::abs(a[j][i]) > max_value) {
				maxi = j;
				max_value = std::abs(a[j][i]);
			}
		}

		// Перестановка строк, чтобы максимальный элемент был на диагонали
		if (maxi != i) {
			std::swap(a[i], a[maxi]);
			std::swap(b[i], b[maxi]);
		}

		// Преобразование текущей строки и строк ниже
		for (int j = i + 1; j < N; j++) {
			double factor = a[j][i] / a[i][i];
			for (int k = i; k < M; k++) {
				a[j][k] -= factor * a[i][k];
			}
			b[j] -= factor * b[i];
		}
	}

	// Обратный ход
	std::vector<double> x(N);
	for (int i = N - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < M; j++) {
			sum += a[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / a[i][i];
	}

	return x;
}