#include "stdio.h"
#include "assert.h"
#include "klu_solve.h"

extern "C"
{
	void lu(int *INT_LU, double *DOUBLE_LU, int N, int NNZ, int NBLOCKS, int NRHS, double *B)
	{
		if (N > MAX_SIZE || NNZ > MAX_NNZ)
			return;

		int r[MAX_SIZE], q[MAX_SIZE], Pnum[MAX_SIZE], Lip[MAX_SIZE], Llen[MAX_SIZE], LUsize[MAX_SIZE], Uip[MAX_SIZE], Ulen[MAX_SIZE], Offp[MAX_SIZE], Offi[MAX_NNZ];
		double Rs[MAX_SIZE], LUbx[MAX_NNZ], UDiag[MAX_SIZE], Offx[MAX_NNZ];

		for (int j = 0; j < N; j++)
		{
			q[j] = INT_LU[j];
			Rs[j] = DOUBLE_LU[j];
		}

		for (int j = 0; j < N; j++)
		{
			Pnum[j] = INT_LU[N + j];
			UDiag[j] = DOUBLE_LU[N + j];
		}

		for (int j = 0; j < N; j++)
			Lip[j] = INT_LU[2 * N + j];

		for (int j = 0; j < N; j++)
			Llen[j] = INT_LU[3 * N + j];

		for (int j = 0; j < N; j++)
			LUsize[j] = INT_LU[4 * N + j];

		for (int j = 0; j < N; j++)
			Uip[j] = INT_LU[5 * N + j];

		for (int j = 0; j < N; j++)
			Ulen[j] = INT_LU[6 * N + j];

		for (int j = 0; j < N + 1; j++)
			r[j] = INT_LU[7 * N + j];

		for (int j = 0; j < N + 1; j++)
			Offp[j] = INT_LU[8 * N + 1 + j];

		for (int j = 0; j < NNZ; j++)
		{
			Offi[j] = INT_LU[9 * N + 2 + j];
			LUbx[j] = DOUBLE_LU[2 * N + j];
		}

		for (int j = 0; j < NNZ; j++)
			Offx[j] = DOUBLE_LU[2 * N + NNZ + j];

		klu_solve(r, q, Pnum, Rs, Lip, Llen, LUbx, LUsize, Uip, Ulen, UDiag, Offp, Offi, Offx, N, NRHS, NBLOCKS, B);
	}
}
