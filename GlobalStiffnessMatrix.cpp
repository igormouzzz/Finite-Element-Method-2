#include "Header.h"
#include "Matrix.h"

Matrix MadeGlobalStiffnessMatrix(int element_count, int node_count, vector<Matrix> matrices, vector<vc> nums)
{
	int size = 2 * node_count;
	Matrix K(size, size);
	int i = 0; int j = 0; int k = 0; int l = 0;

	for (int s = 0; s < element_count; s++)
	{
		for (int n = 0; n < COUNT_OF_NODES; n++)
		{
			for (int m = 0; m < COUNT_OF_NODES; m++)
			{
				K.a[2 * nums[s].n[n]][2 * nums[s].n[m]] = K.a[2 * nums[s].n[n]][2 * nums[s].n[m]] + matrices[s].a[2 * n][2 * m];
				K.a[2 * nums[s].n[n]][2 * nums[s].n[m] + 1] = K.a[2 * nums[s].n[n]][2 * nums[s].n[m] + 1] + matrices[s].a[2 * n][2 * m + 1];
				K.a[2 * nums[s].n[n] + 1][2 * nums[s].n[m]] = K.a[2 * nums[s].n[n] + 1][2 * nums[s].n[m]] + matrices[s].a[2 * n + 1][2 * m];
				K.a[2 * nums[s].n[n] + 1][2 * nums[s].n[m] + 1] = K.a[2 * nums[s].n[n] + 1][2 * nums[s].n[m] + 1] + matrices[s].a[2 * n + 1][2 * m + 1];
			}
		}
	}
	return K;
}