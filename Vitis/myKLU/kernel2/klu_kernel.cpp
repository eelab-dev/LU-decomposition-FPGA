#include "stdio.h"
#include "assert.h"
#include "klu_factor.h"
#include "klu_solve.h"

void read_Ap(int *AP, int *P, int *Q, int *R, int *Ap, klu_symbolic *Symbolic, int N)
{
Read_loop_1:
	for (int i = 0; i < N; i++)
	{
		Ap[i] = AP[i];
		Symbolic->P[i] = P[i];
		Symbolic->Q[i] = Q[i];
		Symbolic->R[i] = R[i];
	}
	Ap[N] = AP[N];
	Symbolic->R[N] = R[N];
}

void read_Ax(int *AI, double *AX, double *LNZ, double *B, int *Ai, double *Ax, double *b, klu_symbolic *Symbolic, int N, int NZ, int nrhs)
{
Read_loop_2:
	for (int i = 0; i < N; i++)
		Symbolic->Lnz[i] = LNZ[i];

	for (int i = 0; i < NZ; i++)
	{
		Ai[i] = AI[i];
		Ax[i] = AX[i];
	}

	// for (int i = 0; i < N * nrhs; i++)
	// 	b[i] = B[i];
}

extern "C"
{
	void lu(int *AP, int *AI, double *AX, int *P, int *Q, int *R, double *LNZ, int N, int NBLOCKS, int MAXBLOCK, int NZOFF, int NZ, int nrhs, double *B)
	{
		if (N > MAX_SIZE || NZ > MAX_NNZ)
			return;

		klu_common Common;
		klu_numeric Numeric;
		klu_symbolic Symbolic;
		klu_defaults(&Common);

		int Ap[MAX_SIZE], Ai[MAX_NNZ];
		double Ax[MAX_NNZ], b[MAX_NNZ];
		read_Ap(AP, P, Q, R, Ap, &Symbolic, N);
		read_Ax(AI, AX, LNZ, B, Ai, Ax, b, &Symbolic, N, NZ, nrhs);

		//        #pragma HLS array_partition variable = Ap factor = 16
		//        #pragma HLS array_partition variable = Ai factor = 16
		//        #pragma HLS array_partition variable = Symbolic.P factor = 16
		//        #pragma HLS array_partition variable = Symbolic.Q factor = 16
		//        #pragma HLS array_partition variable = Symbolic.R factor = 16
		//        #pragma HLS array_partition variable = Ax factor = 16
		//		  #pragma HLS array_partition variable = b type=cyclic factor=16 dim=1
		//        #pragma HLS array_partition variable = Symbolic.Lnz factor = 16
		//        #pragma HLS array_partition variable = Numeric.Pnum type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Offp type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Offi type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Lip type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Uip type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Llen type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Ulen type = cyclic factor = 16
		//      ` #pragma HLS array_partition variable = Numeric.Pinv type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.LUsize type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Offx type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Udiag type = cyclic factor = 16
		// 		  #pragma HLS array_partition variable = Numeric.LUbx type = cyclic factor = 16
		//        #pragma HLS array_partition variable = Numeric.Rs type = cyclic factor = 16
#pragma HLS array_partition variable = Numeric.Xwork type = cyclic factor = 16 dim = 1
#pragma HLS array_partition variable = Numeric.xusolve type = cyclic factor = 16 dim = 1

		Symbolic.n = N;
		Symbolic.nz = NZ;
		Symbolic.nzoff = NZOFF;
		Symbolic.nblocks = NBLOCKS;
		Symbolic.maxblock = MAXBLOCK;

		Numeric.n = Symbolic.n;
		Numeric.nblocks = Symbolic.nblocks;
		Numeric.nzoff = Symbolic.nzoff;
		klu_factor(Ap, Ai, Ax, &Symbolic, &Numeric, &Common);
		klu_solve(&Symbolic, &Numeric, N, nrhs, B, &Common);

		// 	Write_loop:
		// 		for (int i = 0; i < N * nrhs; i++)
		// 		{
		// #pragma HLS pipeline
		// 			B[i] = b[i];
		// 		}
	}
}
