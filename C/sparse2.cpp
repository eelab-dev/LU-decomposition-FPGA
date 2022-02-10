#include <iostream>
#include <mm_malloc.h>
#include <ctime>
#include <chrono>

const int N = 10;

void sparselu(int squareSize, int *Ap, int *Ai, double *Ax, int *Lp, int *Li, double *Lx, int *Up, int *Ui, double *Ux, int lnz, int unz)
{
	for (int i = 0; i < squareSize; i++)
	{
		double b[N] = {0}, x[N] = {0};

		for (int j = Ap[i]; j < Ap[i + 1]; j++)
			b[Ai[j]] = Ax[j];

		for (int j = 0; j < squareSize; j++)
		{
			double sum = 0;

			for (int k = 0; k <= j; k++)
			{
				for (int t = Lp[k]; t < Lp[k + 1]; t++)
				{
					if (Li[t] == j)
						sum += Lx[t] * x[k];
					else if (Li[t] > j)
						t = Lp[k + 1];
				}
			}
			x[j] = b[j] - sum;
		}

		for (int j = Up[i]; j < Up[i + 1]; j++)
		{
			if ((x[Ui[j]] > 1e-6) || (x[Ui[j]] < -1e6))
				Ux[j] = x[Ui[j]];
			else
				Ux[j] = 0;
		}

		for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
		{
			if ((x[Li[j]] > 1e-6) || (x[Li[j]] < -1e6))
				Lx[j] = x[Li[j]] / Ux[Up[i + 1] - 1];
			else
				Ux[j] = 0;
		}
	}
}

void analyse(int *Ap, int *Ai, int *Up, int *Ui, int *Lp, int *Li, double *Lx, int *lnz, int *unz, int n)
{
	for (int i = 0, count = 0; i < n; i++)
	{
		for (int j = Ap[i]; j < Ap[i + 1]; j++)
		{
			if (Ai[j] < i)
			{
				Ui[*unz] = Ai[count];
				(*unz)++;
			}
			else if (Ai[j] == i)
			{
				Li[*lnz] = i;
				Lx[(*lnz)++] = 1;
				Ui[(*unz)++] = i;
			}
			else
			{
				Li[*lnz] = Ai[count];
				(*lnz)++;
			}
			count++;
		}
		Lp[i + 1] = *lnz;
		Up[i + 1] = *unz;
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

	Lp = (int *)_mm_malloc((n + 1) * sizeof(int), 32);
	int data_size = std::size(Ai);
	Li = (int *)_mm_malloc(data_size * sizeof(int), 32);
	Up = (int *)_mm_malloc(data_size * sizeof(int), 32);
	Ui = (int *)_mm_malloc((n + 1) * sizeof(int), 32);
	Lx = (double *)_mm_malloc(data_size * sizeof(double), 32);
	Ux = (double *)_mm_malloc(data_size * sizeof(double), 32);
	std::cout << "Malloc Finished\n";

	std::cout << data_size << std::endl;

	int lnz = 0, unz = 0;

	analyse(Ap, Ai, Up, Ui, Lp, Li, Lx, &lnz, &unz, n);

	// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	// sparselu(n, Ap, Ai, Ax, Lp, Li, Lx, Up, Ui, Ux, lnz, unz);
	// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	for (int i = 0; i < n + 1; i++)
		printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);
	for (int i = n + 1; i < lnz; i++)
		printf("Li[%d]=%d\tLx[%d]=%lf\n", i, Li[i], i, Lx[i]);
	for (int i = n + 1; i < unz; i++)
		printf("Ui[%d]=%d\tUx[%d]=%lf\n", i, Ui[i], i, Ux[i]);

	// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
	return 0;
}