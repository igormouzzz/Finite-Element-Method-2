#include "Header.h"

/*void initialize_array(int n, vector<double>& a, vector<double>& b, vector<double>& c_result)
{
	for (int i = 0; i < n; i++) {
		a[i] = 1.2;
		b[i] = 2.8;
		c_result[i] = a[i] + b[i];
	}
}

void initialize_array2(int n, matrix&m, vector<vector_loc>& b, vector<vector_loc>& c_result)
{
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				m[k][i][j] = 1;
			}
			b[k][i] = k+i;
		}
	}
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < 6; ++i)
		{
			double val = 0;
			for (int j = 0; j < 6; ++j)
			{
				val += m[k][i][j] * b[k][j];
			}
			c_result[k][i] = val;
		}
	}
}

void vector_add(sycl::queue& Q, int n, vector<double>& a, vector<double>& b, vector<double>& c)
{
	sycl::buffer a_buffer(a);
	sycl::buffer b_buffer(b);
	sycl::buffer c_buffer(c);
	auto task_add = Q.submit([&](sycl::handler& cgh) 
	{
		sycl::accessor a_accessor(a_buffer, cgh, sycl::read_only);
		sycl::accessor b_accessor(b_buffer, cgh, sycl::read_only);
		sycl::accessor c_accessor(c_buffer, cgh, sycl::write_only, sycl::no_init);
		cgh.parallel_for(sycl::range<1>(n), [=](sycl::id<1> idx) 
			{
			for(int j = 0; j < 36; j++)
				c_accessor[idx] = a_accessor[idx] + b_accessor[idx];
			});
		});
	task_add.wait();
}

void multiply2(sycl::queue& Q, int n, matrix& m, vector<vector_loc>& b, vector<vector_loc>& c)
{
	sycl::buffer m_buffer(m);
	sycl::buffer b_buffer(b);
	sycl::buffer c_buffer(c);

	auto task_add = Q.submit([&](sycl::handler& cgh)
		{
			sycl::accessor m_accessor(m_buffer, cgh, sycl::read_only);
			sycl::accessor b_accessor(b_buffer, cgh, sycl::read_only);
			sycl::accessor c_accessor(c_buffer, cgh, sycl::write_only, sycl::no_init);
			cgh.parallel_for(sycl::range<1>(n), [=](sycl::id<1> k)
				{
					for (int p = 0; p < 500; p++)
					{
						for (int i = 0; i < 6; i++)
						{
							double val = 0;
							for (int j = 0; j < 6; ++j)
							{
								val += m_accessor[k][i][j] * b_accessor[k][j];
							}
							c_accessor[k][i] = val;
						}
					}
				});
		});
	task_add.wait();
}

int Test()
{
	sycl::queue Q;
	int n = 10000000;
	vector<double> a(n);
	vector<double> b(n);
	vector<double> c_result(n);
	vector<double> c(n);

	initialize_array(n, a, b, c_result);
	vector_add(Q, n, a, b, c);

	return 0;
}

int Test2()
{
	sycl::queue Q;
	int n = 1000000;
	matrix m(n);
	vector<vector_loc> b(n);
	vector<vector_loc> c_result(n);
	vector<vector_loc> c(n);
	for (int k = 0; k < n; k++)
	{
		b[k].resize(6); c_result[k].resize(6); c[k].resize(6);
	}

	initialize_array2(n, m, b, c_result);
	multiply2(Q, n, m, b, c);

}
	return 0;
*/
int main(int argc, char* argv[])
{
	//Example();
	//Task();
	Task2();
	//Task3();

	//Test();
	//Test2();

	return 0;
}