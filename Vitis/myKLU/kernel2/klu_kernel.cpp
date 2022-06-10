#include "stdio.h"
#include "assert.h"
#include "klu_solve.h"

extern "C"
{
	void lu(int *R, int *Q, int *PNUM, double *RS, int *LIP, int *LLEN, double *LUBX, int *LUSIZE, int *UIP, int *ULEN, double *UDIAG, int *OFFP, int *OFFI, double *OFFX, int N, int NBLOCKS, int NRHS, double *B)
	{
		if (N > MAX_SIZE)
			return;

		klu_solve(R, Q, PNUM, RS, LIP, LLEN, LUBX, LUSIZE, UIP, ULEN, UDIAG, OFFP, OFFI, OFFX, N, NRHS, NBLOCKS, B);
	}
}
