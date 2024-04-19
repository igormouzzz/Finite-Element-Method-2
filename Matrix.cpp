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

Matrix::Matrix(vector<row> rows)
{
	this->N = rows.size();
	this->M = rows[0].size();
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
			this->a[i][j] = rows[i][j];
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
//#pragma omp parallel for
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
//#pragma omp parallel for
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
//#pragma omp parallel for
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
//#pragma omp parallel for
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
//#pragma omp parallel for
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
const vector<double> Matrix::operator*(const vector<double>& b)
{
	const vector<double> s = (*this) * b;
	return s;
}
Matrix Matrix::T()
{
	Matrix S(M, N);
//#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			S.a[i][j] = a[j][i];
		}
	}
	return S;
}
void Matrix::UnitRowAndColumn(int k, int l, int n)
{
	for (int i = 0; i < N; i++)
	{
		if (i == k)
		{
			a[i][l] = 1.0/double(n);
		}
		else
		{
			a[i][l] = 0;
		}
	}
	for (int j = 0; j < M; j++)
	{
		if (j == l)
		{
			a[k][j] = 1.0/double(n);
		}
		else
		{
			a[k][j] = 0;
		}
	}
}
void Matrix::ZeroRowAndColumn(int k, int l)
{
	for (int i = 0; i < N; i++)
	{
		if (i == k)
		{
			a[i][l] = 0;
		}
		else
		{
			a[i][l] = 0;
		}
	}
	for (int j = 0; j < M; j++)
	{
		if (j == l)
		{
			a[k][j] = 0;
		}
		else
		{
			a[k][j] = 0;
		}
	}
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

double Matrix::SumOfComponentsForProduct(vector<double>& b)
{
	vector<double> x = (*this) * b;
	double sum = 0;
	for (int i = 0; i < x.size(); i++) sum += x[i];
	return sum;
}

void Matrix::ToCSR(ofstream& f)
{
	unsigned int non_zero = 0;
	unsigned int non_zero_in_rows = 0;
	unsigned int non_zero_in_columns = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			if (abs(a[i][j]) > 1e-5)
			{
				non_zero++;
			}
		}
	}
	//f << N << " " << non_zero << endl;
	
	//f << 0;
	for (int i = 0; i < N; i++)
	{
		//non_zero_in_rows = 0;
		for (int j = 0; j < M; j++)
		{
			if (abs(a[i][j]) > 1e-5)
			{
				non_zero_in_rows++;
			}
		}
		//f << " " << non_zero_in_rows;
	}
	//f << endl;

	//f << 0;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < N; i++)
		{
			if (abs(a[i][j]) > 1e-5)
			{
				//f << " " << j;
			}
		}
	}
	//f << endl;
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			if (abs(a[i][j]) > 1e-5)
			{
				//f << a[i][j] << " ";
			}
		}
	}
	//f << endl;
}

vector<double> Matrix::Gauss(vector<double>& b)
{
	if (norm_square(b) < 1e-30) return vector<double>(b.size());

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

vector<double> Matrix::CG(vector<double>& b)
{
	vector<double> Xk(M), Zk(M), Rk(M), Sz(M);
	double alpha, beta, mf;
	double Spr, Spr1, Spz;

	int kl = 0;
	int i, j;
	double max_iter = 1000;
	/* Вычисляем сумму квадратов элементов вектора b*/
	for (mf = 0, i = 0; i < M; i++) 
	{
		mf += b[i] * b[i];
	}


	/* Задаем начальное приближение корней. В Хk хранятся значения корней
	 * к-й итерации. */
	for (i = 0; i < M; i++) 
	{
		Xk[i] = 0;
	}

	/* Задаем начальное значение r0 и z0. */
	for (i = 0; i < M; i++) {
		for (Sz[i] = 0, j = 0; j < M; j++)
		{
			Sz[i] += a[i][j] * Xk[j];
		}
		Rk[i] = b[i] - Sz[i];
		Zk[i] = Rk[i];
	}
	int Iteration = 0;
	do {
		Iteration++;
		/* Вычисляем числитель и знаменатель для коэффициента
		 * alpha = (rk-1,rk-1)/(Azk-1,zk-1) */
		Spz = 0;
		Spr = 0;
		for (i = 0; i < M; i++) {
			for (Sz[i] = 0, j = 0; j < M; j++) 
			{
				Sz[i] += a[i][j] * Zk[j];
			}
			Spz += Sz[i] * Zk[i];
			Spr += Rk[i] * Rk[i];
		}
		alpha = Spr / Spz;             /*  alpha    */


		/* Вычисляем вектор решения: xk = xk-1+ alpha * zk-1,
			вектор невязки: rk = rk-1 - alpha * A * zk-1 и числитель для betaa равный (rk,rk) */
		Spr1 = 0;
		for (i = 0; i < M; i++) 
		{
			Xk[i] += alpha * Zk[i];
			Rk[i] -= alpha * Sz[i];
			Spr1 += Rk[i] * Rk[i];
			//cout << "Iter #" << kl;
			//cout << " " << "X[" << i << "] = " << Xk[i] << endl;
		}
		//cout << endl;
		kl++;

		/* Вычисляем  beta  */
		beta = Spr1 / Spr;

		/* Вычисляем вектор спуска: zk = rk+ beta * zk-1 */
		for (i = 0; i < M; i++)
			Zk[i] = Rk[i] + beta * Zk[i];
	}
	/* Проверяем условие выхода из итерационного цикла  */
	while (Spr1 / mf > eps * eps && Iteration < max_iter);

	cout << "kol-vo iter: " << kl << endl;

	return Xk;
}

double norm_square(vector<double>& b)
{
	double norm = 0;
	#pragma omp parallel for reduction(+:norm)
	for (int i = 0; i < b.size(); i++) norm += b[i] * b[i];
	return norm;
}

double norm_l1(vector<double>& b)
{
	double norm = 0;
	#pragma omp parallel for reduction(+:norm)
	for (int i = 0; i < b.size(); i++) norm += abs(b[i]);
	return norm;
}

double scalar_product(vector<double>& a, vector<double>& b)
{
	double sp = 0;
	#pragma omp parallel for reduction(+:sp)
	for (int i = 0; i < a.size(); i++) sp += a[i] * b[i];
	return sp;
}

vector<double> Matrix::CG2(vector<double>& b, double tolerance)
{
	double alpha, beta;
	vector<double> Xk(N), Xk_1(N), Rk(N), Rk_1(N), Zk(N), Zk_1(N), AXk_1(N);
	
	unsigned int iteration = 0;
	unsigned int max_iter = 1000;

	vector<double> tmp(N);

	for (int i = 0; i < N; i++)
	{
		Xk_1[i] = 0;
	}

	AXk_1 = (*this) * Xk_1;
	for (int i = 0; i < N; i++)
	{
		Rk_1[i] = b[i] - AXk_1[i];
	}
	Zk_1 = Rk_1;

	do
	{
		tmp = (*this) * Zk_1;
		alpha = norm_square(Rk_1) / scalar_product(tmp, Zk_1);

		for (int i = 0; i < N; i++)
		{
			Xk[i] = Xk_1[i] + alpha * Zk_1[i];
			Rk[i] = Rk_1[i] - alpha * tmp[i];
		}

		beta = norm_square(Rk) / norm_square(Rk_1);

		for (int i = 0; i < N; i++)
		{
			Zk[i] = Rk[i] + beta * Zk_1[i];
		}
		iteration++;

		cout << norm_square(Rk) / norm_square(b) << endl;

		Xk_1 = Xk; Rk_1 = Rk; Zk_1 = Zk; AXk_1 = (*this) * Xk;

	} while (norm_square(Rk_1)/ norm_square(b) >= tolerance * tolerance && iteration < max_iter);

	cout << "iterations: " << iteration << endl;
	cout << "tolerance: " << tolerance << endl;

	return Xk;
}

vector<double> Matrix::CG3(vector<double>& b)
{
	double b_norm = norm_square(b);
	//double b_norm = norm_l1(b);

	if (b_norm == 0) return vector<double>(b.size());

	size_t n = N;

	std::vector<double> x(n, 0.0);  // Initial guess for the solution
	std::vector<double> r = b;      // Residual vector
	std::vector<double> p = r;      // Search direction vector

	unsigned int iteration = 0;
	unsigned int max_iter = 1000;

	do
	{
		vector<double> Ap = (*this)*p;
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
	} while (norm_square(r) / b_norm > eps * eps && iteration < max_iter);

	cout << "iterations: " << iteration << endl;
	//cout << "tolerance: " << eps << endl;

	return x;
}