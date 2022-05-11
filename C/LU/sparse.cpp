#include <iostream>
#include <mm_malloc.h>
#include <ctime>
#include <chrono>

const int N = 10;

void sparselu(int squareSize, int *Ap, int *Ai, double *Ax, int *Lp, int *Li, double *Lx, int *Up, int *Ui, double *Ux)
{
	int countL = 1, countU = 0;

	for (int i = 0; i < squareSize; i++)
	{
		Lp[i] = i;
		Li[i] = i;
		Lx[i] = 1;
	}

	for (int i = 0; i < squareSize; i++)
	{
		double b[N] = {0}, x[N] = {0};

		if (Ap[i] != Ap[i + 1])
		{
			for (int j = Ap[i]; j < Ap[i + 1]; j++)
				b[Ai[j]] = Ax[j];
		}
		else
		{
			// Lp[i] = 0;
			Up[i] = 0;
			continue;
		}
		for (int j = 0; j < squareSize; j++)
		{
			double sum = 0;
			for (int k = 0; k <= j; k++)
			{

				for (int t = Lp[k]; t < Lp[k + 1]; t++)
					if (Li[t] == j)
						sum += Lx[Li[t]] * x[k];
			}
			x[j] = b[j] - sum;
		}

		for (int j = 0; j <= i; j++)
		{
			if (std::abs(x[j]) < 1e-6)
				continue;
			else
			{
				Ui[countU] = j;
				Ux[countU++] = x[j];
			}
		}
		Up[i + 1] = countU;

		for (int j = i + 1; j < squareSize; j++)
		{
			if (std::abs(x[j]) < 1e-6)
				continue;
			else
			{
				Li[countL] = j;
				Lx[countL++] = x[j] / Ux[Up[i + 1] - 1];
			}
		}
		Lp[i + 1] = countL;
		if (i < squareSize - 1)
		{
			Li[countL] = i + 1;
			Lx[countL++] = 1;
		}
	}
}

int main(int argc, char *argv[])
{
	int n = 10;
	int Ap[] = {0, 2, 4, 7, 9, 13, 14, 18, 21, 23, 24};
	int Ai[] = {0, 7, 1, 5, 0, 2, 7, 3, 8, 3, 4, 5, 8, 5, 0, 6, 7, 8, 4, 5, 7, 2, 8, 9};
	double Ax[] = {8, 8, 2, 4, 10, 3, 10, 1, 5, 2, 1, 4, 10, 2, 3, 5, 3,
				   15, 4, 16, 3, 7, 2, 9};
	int *Lp, *Li, *Up, *Ui;
	double *Lx, *Ux;

	Lp = (int *)_mm_malloc(n * (n + 1) / 2 * sizeof(int), 32);
	Li = (int *)_mm_malloc(n * (n + 1) / 2 * sizeof(int), 32);
	Up = (int *)_mm_malloc(n * (n + 1) / 2 * sizeof(int), 32);
	Ui = (int *)_mm_malloc(n * (n + 1) / 2 * sizeof(int), 32);
	Lx = (double *)_mm_malloc(n * (n + 1) / 2 * sizeof(double), 32);
	Ux = (double *)_mm_malloc(n * (n + 1) / 2 * sizeof(double), 32);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	sparselu(n, Ap, Ai, Ax, Lp, Li, Lx, Up, Ui, Ux);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	for (int i; i < 9; i++)
		printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);

	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
	return 0;
}