#include <mm_malloc.h>
#include "cholmod.h"
#include <chrono>
#include <iostream>
#include <vector>

#define N 90

void sparselu(int squareSize, int *Ap, int *Ai, double *Ax, int *Lp, int *Li, double *Lx, int *Up, int *Ui, double *Ux, int lnz, int unz)
{
	for (int i = 0; i < squareSize; i++)
	{
		double x[N] = {0};

		for (int j = Ap[i]; j < Ap[i + 1]; j++)
			x[Ai[j]] = Ax[j];

		for (int j = 0; j < squareSize; j++)
		{
			for (int k = 0; k <= j; k++)
			{
				for (int t = Lp[k] + 1; t < Lp[k + 1]; t++)
				{
					if (Li[t] == j)
						x[j] -= Lx[t] * x[k];
					else if (Li[t] > j)
						t = Lp[k + 1];
				}
			}
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
	cholmod_sparse *A;
	cholmod_common ch;
	cholmod_start(&ch);
	// FILE *f = fopen("../Matrix_Sample/impcol_a.mtx", "r");
	A = cholmod_read_sparse(stdin, &ch);

	if (A)
	{
		int runtime = 10;
		std::chrono::steady_clock::time_point begin, end;
		long total = 0;
		int *Ap, *Ai, *Lp, *Li, *Up, *Ui;
		double *Lx, *Ux, *Ax;
		int lnz, unz, n;

		for (int i = 0; i < 1; i++)
		{
			int data_size = A->nzmax;
			std::cout << "data_size=" << data_size << std::endl;
			Ap = (int *)(A->p);
			Ai = (int *)(A->i);
			Ax = (double *)(A->x);
			std::vector<int> Lp(A->nrow + 1);
			// Lp = (int *)_mm_malloc((A->nrow + 1) * sizeof(int), 32);
			Li = (int *)_mm_malloc(data_size * sizeof(int), 32);
			Up = (int *)_mm_malloc(data_size * sizeof(int), 32);
			Ui = (int *)_mm_malloc((A->nrow + 1) * sizeof(int), 32);
			Lx = (double *)_mm_malloc(data_size * sizeof(double), 32);
			Ux = (double *)_mm_malloc(data_size * sizeof(double), 32);
			printf("Malloc Finished\n");

			lnz = 0;
			unz = 0;
			n = A->nrow;
			analyse(Ap, Ai, Up, Ui, Lp.data(), Li, Lx, &lnz, &unz, n);
			printf("Analysis Finished,%d\n", ((int *)(A->p))[n]);

			begin = std::chrono::steady_clock::now();
			sparselu(n, Ap, Ai, Ax, Lp.data(), Li, Lx, Up, Ui, Ux, lnz, unz);
			end = std::chrono::steady_clock::now();
			total += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		}

		// for (int i = 0; i < A->nrow + 1; i++)
		// 	printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);
		// for (int i = A->nrow + 1; i < lnz; i++)
		// 	printf("Li[%d]=%d\tLx[%d]=%lf\n", i, Li[i], i, Lx[i]);
		// for (int i = A->nrow + 1; i < unz; i++)
		// 	printf("Ui[%d]=%d\tUx[%d]=%lf\n", i, Ui[i], i, Ux[i]);

		printf("Time:%lfns\n", total / 10.0);
		cholmod_free_sparse(&A, &ch);
	}

	cholmod_finish(&ch);
	return 0;
}