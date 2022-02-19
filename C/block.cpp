#include <iostream>
#include <fstream>
#include <mm_malloc.h>
#include "timer.h"

using namespace std;

#define BlockSize 5
#define BlockGrainSize 8 //BlockGrainSize >= BlockSize; Good if BlockGrainSize < 2*BlockSize

double *ReadMatrixFromFile(char *fileName, int &N)
{
	FILE *pFile = fopen(fileName, "r");
	if (!pFile)
	{
		cout << "Input file not found\n";
		return 0;
	}

	int h, w;
	fscanf(pFile, "%d %d", &h, &w);
	if (h != w)
	{
		fclose(pFile);
		printf("Error!\n");
		return 0;
	}
	else
	{
		N = h;
	}

	double *A = (double *)_mm_malloc(N * N * sizeof(double), 32);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fscanf(pFile, "%lf", &A[i * N + j]);
		}
	}

	fclose(pFile);

	return A;
}

void WriteLMatrixToFile(char *fileName, const double *L, int N)
{
	FILE *pFile = fopen(fileName, "w");
	if (!pFile)
	{
		printf("File %s not found!\n", fileName);
		return;
	}

	fprintf(pFile, "%d %d\n", N, N);
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			fprintf(pFile, "%.16f ", L[i * (i + 1) / 2 + j]);
		}
		for (int j = i + 1; j < N - 1; j++)
		{
			fprintf(pFile, "%.16f ", 0.0);
		}
		fprintf(pFile, "%.16f\n", 0.0);
	}
	for (int j = 0; j < N - 1; j++)
	{
		fprintf(pFile, "%.16f ", L[(N - 1) * N / 2 + j]);
	}
	fprintf(pFile, "%.16f", 1.0);

	fclose(pFile);
}

void WriteUMatrixToFile(char *fileName, const double *U, int N)
{
	FILE *pFile = fopen(fileName, "w");
	if (!pFile)
	{
		printf("File %s not found!\n", fileName);
		return;
	}

	fprintf(pFile, "%d %d\n", N, N);
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < i; j++)
		{
			fprintf(pFile, "%.16f ", 0.0);
		}
		for (int j = i; j < N - 1; j++)
		{
			fprintf(pFile, "%.16f ", U[j * (j + 1) / 2 + i]);
		}
		fprintf(pFile, "%.16f\n", U[(N - 1) * N / 2 + i]);
	}
	for (int j = 0; j < N - 1; j++)
	{
		fprintf(pFile, "%.16f ", 0.0);
	}
	fprintf(pFile, "%.16f", U[(N - 1) * N / 2 + N - 1]);

	fclose(pFile);
}

void WriteTimeToFile(char *fileName, double time)
{
	FILE *pFile = fopen(fileName, "w");
	if (!pFile)
	{
		printf("File %s not found!\n", fileName);
		return;
	}

	fprintf(pFile, "%.16f", time);

	fclose(pFile);
}

void DiagonalSubmatrixLUDecompose(int bias, int squareSize, int matrixSize, int *Ap, int *Ai, double *Ax, int *Lp, int *Li, double *Lx, int *Up, int *Ui, double *Ux)
{
	int N = bias + squareSize;

	int i, j, k, countL = 1, countU = 0;

	double sum;

	for (int i = 0; i < squareSize; i++)
	{
		Lp[i] = i;
		Li[i] = i;
		Lx[i] = 1;
	}

	for (i = 0; i < squareSize; i++)
	{
		double b[squareSize], x[squareSize];
		for (int j = 0; j < squareSize; j++)
		{
			b[j] = 0;
			x[j] = 0;
		}

		std::cout << i << std::endl;
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
			for (int k = 0; k <= j; k++)
			{
				double sum = 0;
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

void SolveUpperEquation(int bias, int squareSize, int matrixSize, int *Ap, int *Ai, double *Ax, int *Lp, int *Li, double *Lx, int *Up, int *Ui, double *Ux)
{
	int i, j, k, countL, countU;

	for (int i = bias + squareSize; i < bias + 2 * matrixSize; i++)
	{
		for (int j = Ap[i]; j < Ap[i + 1]; j++)
		{
			if (Ai[j] > bias + 2 * matrixSize)
				break;
			if (Ai[j] == bias)
			{
				Up[i - squareSize + 1] += 1;
				Ui[countU] = bias;
				Ux[countU++] = Ax[j];
			}
		}

		for (int j = Ap[i]; j < Ap[i + 1]; j++) // column direction
		{
			double sum = 0;
			for (int k = Ai[j]; k < i - squareSize; k++)
			{
				for (int t; t < k; t++)
				{
				}
			}
		}
	}

	for (j = bias + squareSize; j < matrixSize; j++)
	{
		U[j * (j + 1) / 2 + bias] = A[bias * matrixSize + j];
		for (i = bias + 1; i < bias + squareSize; i++)
		{
			sum = 0.0;
			for (k = bias; k < i; k++)
			{
				sum += U[j * (j + 1) / 2 + k] * L[i * (i + 1) / 2 + k];
			}
			U[j * (j + 1) / 2 + i] = A[i * matrixSize + j] - sum;
		}
	}
}

void SolveLeftEquation(int bias, int squareSize, int matrixSize, const double *A, double *L, const double *U)
{
	int i, j, k;
	double sum;

	for (i = bias + squareSize; i < matrixSize; i++)
	{
		L[i * (i + 1) / 2 + bias] = A[i * matrixSize + bias] / U[bias * (bias + 3) / 2];
		for (j = bias + 1; j < bias + squareSize; j++)
		{
			sum = 0;
			for (k = bias; k < j; k++)
			{
				sum += L[i * (i + 1) / 2 + k] * U[j * (j + 1) / 2 + k];
			}
			L[i * (i + 1) / 2 + j] = (A[i * matrixSize + j] - sum) / U[j * (j + 3) / 2];
		}
	}
}

void UpdateDiagonalSubmatrix(int bias, int squareSize, int matrixSize, double *A, double *L, double *U)
{
	int i, j, k;
	double sum;

	for (i = bias + squareSize; i < matrixSize; i++)
	{
		for (j = bias + squareSize; j < matrixSize; j++)
		{
			sum = 0;
			for (k = bias; k < bias + squareSize; k++)
			{
				sum += L[i * (i + 1) / 2 + k] * U[j * (j + 1) / 2 + k];
			}
			A[i * matrixSize + j] -= sum;
		}
	}
}

// void LUDecompose(int matrixSize, double *A, double *L, double *U)
// {
// 	for (int bias = 0; bias < matrixSize; bias += BlockSize)
// 	{
// 		if (matrixSize - bias <= BlockGrainSize)
// 		{
// 			DiagonalSubmatrixLUDecompose(bias, matrixSize - bias, matrixSize, A, L, U);
// 			std::cout << "Entered" << std::endl;
// 			break;
// 		}
// 		DiagonalSubmatrixLUDecompose(bias, BlockSize, matrixSize, A, L, U);
// 		SolveUpperEquation(bias, BlockSize, matrixSize, A, L, U);
// 		SolveLeftEquation(bias, BlockSize, matrixSize, A, L, U);
// 		UpdateDiagonalSubmatrix(bias, BlockSize, matrixSize, A, L, U);
// 	}
// }

int main(int argc, char *argv[])
{
	int n = 3;
	int Ap[] = {0, 2, 3, 5};
	int Ai[] = {0, 2, 1, 1, 2};
	double Ax[] = {2, 6, 5, 2, 1};
	int *Lp, *Li, *Up, *Ui;
	double *Lx, *Ux;

	Lp = (int *)_mm_malloc(9 * sizeof(int), 32);
	Li = (int *)_mm_malloc(9 * sizeof(int), 32);
	Up = (int *)_mm_malloc(9 * sizeof(int), 32);
	Ui = (int *)_mm_malloc(9 * sizeof(int), 32);
	Lx = (double *)_mm_malloc(9 * sizeof(double), 32);
	Ux = (double *)_mm_malloc(9 * sizeof(double), 32);

	// if (argc < 5)
	// {
	// 	printf("Error!\n");
	// 	return 0;
	// }

	DiagonalSubmatrixLUDecompose(0, 3, 3, Ap, Ai, Ax, Lp, Li, Lx, Up, Ui, Ux);

	std::cout << "Finished" << std::endl;

	for (int i; i < 9; i++)
		printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);

	// double *A, *L, *U;
	// int N;

	// A = ReadMatrixFromFile(argv[1], N);
	// if (!A)
	// {
	// 	return 0;
	// }

	// std::cout << N << std::endl;
	// int triangularMatrixSize = N * (N + 1) / 2;
	// L = (double *)_mm_malloc(triangularMatrixSize * sizeof(double), 32);
	// U = (double *)_mm_malloc(triangularMatrixSize * sizeof(double), 32);

	// Timer timer;
	// timer.start();
	// LUDecompose(N, A, L, U);
	// timer.stop();

	// WriteLMatrixToFile(argv[2], L, N);
	// WriteUMatrixToFile(argv[3], U, N);
	// WriteTimeToFile(argv[4], timer.getElapsed());

	// _mm_free(A);
	// _mm_free(L);
	// _mm_free(U);

	return 0;
}