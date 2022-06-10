#include "stdio.h"
#include "assert.h"
#include "klu_solve.h"

extern "C"
{
	void lu(int *R, int *Q, int *PNUM, double *RS, int *LIP, int *LLEN, double *LUBX, int *LUSIZE, int *UIP, int *ULEN, double *UDIAG, int *OFFP, int *OFFI, double *OFFX, int N, int NNZ, int NBLOCKS, int NRHS, double *B)
	{
		if (N > MAX_SIZE)
			return;

		int r[MAX_SIZE], q[MAX_SIZE], Pnum[MAX_SIZE], Lip[MAX_SIZE], Llen[MAX_SIZE], LUsize[MAX_SIZE], Uip[MAX_SIZE], Ulen[MAX_SIZE], Offp[MAX_SIZE], Offi[MAX_NNZ];
		double Rs[MAX_SIZE], LUbx[MAX_NNZ], UDiag[MAX_SIZE], Offx[MAX_NNZ];

		for (int i = 0; i < N; i++)
		{
			r[i] = R[i];
			q[i] = Q[i];
			Pnum[i] = PNUM[i];
			Lip[i] = LIP[i];
			Llen[i] = LLEN[i];
			LUsize[i] = LUSIZE[i];
			Uip[i] = UIP[i];
			Ulen[i] = ULEN[i];
			Offp[i] = OFFP[i];
			Rs[i] = RS[i];
			UDiag[i] = UDIAG[i];
		}
		Offp[N] = OFFP[N];
		for (int i = 0; i < NNZ; i++)
		{
			Offi[i] = OFFI[i];
			LUbx[i] = LUBX[i];
			Offx[i] = OFFX[i];
		}

		klu_solve(r, q, Pnum, Rs, Lip, Llen, LUbx, LUsize, Uip, Ulen, UDiag, Offp, Offi, Offx, N, NRHS, NBLOCKS, B);
	}
}
