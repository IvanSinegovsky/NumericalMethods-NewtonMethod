#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f1(double x1, double x2)
{
	return (x1 * x1 * x1 + x2 * x2 * x2 - x1 * 6 + 3);
}
double f2(double x1, double x2)
{
	return (x1 * x1 * x1 - x2 * x2 * x2 - 6 * x2 + 2);
}
double df11(double x1, double x2)
{
	return 3 * x1 * x1 - 6;
}
double df12(double x1, double x2)
{
	return 3 * x2;
}
double df21(double x1, double x2)
{
	return 3 * x1 * x1;
}
double df22(double x1, double x2)
{
	return -3 * x2 * x2 - 6;
}

double** ArrayEdit(int n)
{
	double** v = new double* [n];

	for (int i = 0; i < n; i++)
		v[i] = new double[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			v[i][j] = 0;

	return v;
}

bool Gauss(double** a1, int n, double* x);

double IterativeComputationMethod(double* x, int n, int i, int j)
{
	double m = 0.01; //относительное приращение
	double* x1 = new double[n];

	for (int k = 0; k < n; k++)
		x1[k] = x[k];

	x1[j] += x1[j] * m;
	if (i == 0)
		return (f1(x1[0], x1[1]) - f1(x[0], x[1])) / m / x[j];
	if (i == 1)
		return (f2(x1[0], x1[1]) - f2(x[0], x[1])) / m / x[j];
}

void Newton(double* x, int n)
{
	x[0] = 1;
	x[1] = 1;
	double d = 1e-9;
	int NIT = 100; //максимальное количество операций
	double** a = new double* [n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n + 1];
	double* dx = new double[n]; //delta x
	double* F = new double[n]; //значение функции при х
	double** J = ArrayEdit(n);

	cout << "     l      d1      d2" << endl << endl;

	for (int l = 1; abs(f1(x[0], x[1]) + f2(x[0], x[1])) > d || abs(dx[1] + dx[0]) > d && l < NIT; l++)
	{
		F[0] = -f1(x[0], x[1]);
		F[1] = -f2(x[0], x[1]);
		//точное вычисление матрицы Якоби
		J[0][0] = df11(x[0], x[1]);
		J[0][1] = df12(x[0], x[1]);
		J[1][0] = df21(x[0], x[1]);
		J[1][1] = df22(x[0], x[1]);

		//приближенное вычисление матрицы Якоби
		/*for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				J[i][j] = IterativeComputationMethod(x, n, i, j);*/

		for (int i = 0; i < n; i++) //заполнение матрицы для вызова метода Гаусса
			for (int j = 0; j < n; j++)
				a[i][j] = J[i][j];

		for (int i = 0; i < n; i++) //заполнение последнего столбца b
			a[i][n] = F[i];

		if (Gauss(a, n, dx))
			for (int i = 0; i < n; i++)
				x[i] += dx[i];
		else
			break;
		if(f1(x[0], x[1]) > f2(x[0], x[1]))
			cout << setw(6) << l << " " << setw(6) << f1(x[0], x[1])
				<< " " << setw(6) << abs(dx[1]) << endl << endl;
		else
			cout << setw(6) << l << " " << setw(6) << f2(x[0], x[1])
			<< " " << setw(6) << abs(dx[1]) << endl << endl;
	}

	for (int i = 0; i < n; i++)
		delete[]a[i];
	delete[]a;
	delete[]dx;
	delete[]F;
	for (int i = 0; i < n; i++)
		delete[]J[i];
	delete[]J;
}

int main()
{
	const int n = 2;
	double* x = new double[n];
	Newton(x, n);
	cout << "X :" << endl;
	for (int i = 0; i < n; i++)
		cout << x[i] << endl;
}