#include <stdio.h>
#include <mm_malloc.h>
#include "cholmod.h"
#include <time.h>

#define N 39

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
	cholmod_sparse *A;
	cholmod_common ch;
	cholmod_start(&ch);
	// FILE *f = fopen("../Matrix_Sample/impcol_a.mtx", "r");
	A = cholmod_read_sparse(stdin, &ch);

	if (A)
	{
		int runtime = 10;
		clock_t begin, end, total = 0;
		int *Lp, *Li, *Up, *Ui;
		double *Lx, *Ux;

		// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		for (int i = 0; i < runtime; i++)
		{
			Lp = (int *)_mm_malloc((A->nrow + 1) * sizeof(int), 32);
			int data_size = A->nzmax;
			Li = (int *)_mm_malloc(data_size * sizeof(int), 32);
			Up = (int *)_mm_malloc(data_size * sizeof(int), 32);
			Ui = (int *)_mm_malloc((A->nrow + 1) * sizeof(int), 32);
			Lx = (double *)_mm_malloc(data_size * sizeof(double), 32);
			Ux = (double *)_mm_malloc(data_size * sizeof(double), 32);
			printf("Malloc Finished\n");

			int lnz = 0, unz = 0, n = A->nrow;
			analyse(A->p, A->i, Up, Ui, Lp, Li, Lx, &lnz, &unz, n);
			printf("Analysis Finished,%d\n", ((int *)(A->p))[n]);

			begin = clock();
			sparselu(n, A->p, A->i, A->x, Lp, Li, Lx, Up, Ui, Ux, lnz, unz);
			end = clock();
			total += end - begin;
		}
		// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

		// for (int i = 0; i < A->nrow + 1; i++)
		// 	printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);
		// for (int i = A->nrow + 1; i < lnz; i++)
		// 	printf("Li[%d]=%d\tLx[%d]=%lf\n", i, Li[i], i, Lx[i]);
		// for (int i = A->nrow + 1; i < unz; i++)
		// 	printf("Ui[%d]=%d\tUx[%d]=%lf\n", i, Ui[i], i, Ux[i]);

		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
		printf("Time:%lfus\n", total / 10.0);
		cholmod_free_sparse(&A, &ch);
	}

	cholmod_finish(&ch);
	return 0;
}