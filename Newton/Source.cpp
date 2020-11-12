#include <cmath>
#include <iostream>

bool Gauss(double** a1, int n, double* x)
{
	int f = 0;
	double** a = new double* [n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n + 1];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n + 1; j++)
			a[i][j] = a1[i][j];

	for (int k = 0; k < n; k++)
	{
		double pivot = a[k][k];

		int q = k;

		for (int i = k; i < n; i++)
			if (abs(pivot) < abs(a[i][k]))
			{
				pivot = a[i][k];
				q = i;
			}

		if (abs(pivot) < 1e-12)
		{
			f = 1;
			break;
		}

		if (q != k)
			std::swap(a[q], a[k]);

		double z = 0;

		z = a[k][k];

		for (int j = k; j < n + 1; j++)
			a[k][j] /= z;

		for (int i = k + 1; i < n; i++)
		{
			z = a[i][k];
			for (int j = k; j < n + 1; j++)
				a[i][j] -= z * a[k][j];
		}
	}

	double sum = 0;

	if (f == 0)
	{
		x[n - 1] = a[n - 1][n];

		for (int k = n - 2; k >= 0; k--)
		{
			sum = 0;

			for (int j = k + 1; j < n; j++)
				sum += a[k][j] * x[j];

			x[k] = a[k][n] - sum;
		}
	}
	else
	{
		for (int i = 0; i < n; i++)
			delete[]a[i];
		delete[]a;

		std::cout << "error :" << f << std::endl;
		return 0;
	}

	for (int i = 0; i < n; i++)
		delete[]a[i];
	delete[]a;

	return 1;
}