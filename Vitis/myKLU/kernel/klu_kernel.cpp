#include "stdio.h"
#include "assert.h"
#include "klu_factor.h"
#include "klu_solve.h"

void read_Ap(int *AP, int *P, int *Q, int *R, double *LNZ, int *Ap, klu_symbolic *Symbolic, int N)
{
Read_loop_1:
	for (int i = 0; i < N; i++)
	{
		Ap[i] = AP[i];
		Symbolic->P[i] = P[i];
		Symbolic->Q[i] = Q[i];
		Symbolic->R[i] = R[i];
		Symbolic->Lnz[i] = LNZ[i];
	}
	Ap[N] = AP[N];
	Symbolic->R[N] = R[N];
}

void read_Ax(int *AI, double *AX, int *Ai, double *Ax, int NZ)
{
Read_loop_2:
	for (int i = 0; i < NZ; i++)
	{
		Ai[i] = AI[i];
		Ax[i] = AX[i];
	}
}

void solve_b(klu_symbolic *Symbolic, klu_numeric *Numeric, int N, int nrhs, double *b, klu_common *Common)
{
	klu_solve(Symbolic, Numeric, N, nrhs, b, Common);
}

extern "C"
{
	void lu(int *AP, int *AI, double *AX, int *P, int *Q, int *R, double *LNZ, int N, int NBLOCKS, int MAXBLOCK, int NZOFF, int NZ, int nrhs, double *B)
	{
		if(N>MAX_SIZE || NZ>MAX_NNZ)
			return;

		klu_common Common;
		klu_numeric Numeric;
		klu_symbolic Symbolic;
		klu_defaults(&Common);

#pragma HLS INTERFACE m_axi depth = 2048 port = AP max_read_burst_length = 64 offset = slave bundle = gmem
#pragma HLS INTERFACE m_axi depth = 2048 port = AI max_read_burst_length = 64 offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi depth = 2048 port = AX max_read_burst_length = 64 offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi depth = 2048 port = P max_read_burst_length = 64 offset = slave bundle = gmem3
#pragma HLS INTERFACE m_axi depth = 2048 port = Q max_read_burst_length = 64 offset = slave bundle = gmem4
#pragma HLS INTERFACE m_axi depth = 2048 port = R max_read_burst_length = 64 offset = slave bundle = gmem5
#pragma HLS INTERFACE m_axi depth = 2048 port = LNZ max_read_burst_length = 64 offset = slave bundle = gmem6
#pragma HLS INTERFACE m_axi depth = 2048 port = B max_read_burst_length = 64 offset = slave bundle = gmem7

		int Ap[MAX_SIZE], Ai[MAX_NNZ];
		double Ax[MAX_NNZ], b[MAX_NNZ];
		read_Ap(AP, P, Q, R, LNZ, Ap, &Symbolic, N);
		read_Ax(AI, AX, Ai, Ax, NZ);
		for (int i = 0; i < N * nrhs; i++)
			b[i] = B[i];
		//        #pragma HLS array_partition variable = Ap factor = 16
		//        #pragma HLS array_partition variable = Ai factor = 16
		//        #pragma HLS array_partition variable = Symbolic.P factor = 16
		//        #pragma HLS array_partition variable = Symbolic.Q factor = 16
		//        #pragma HLS array_partition variable = Symbolic.R factor = 16
		//        #pragma HLS array_partition variable = Ax factor = 16
		//        #pragma HLS array_partition variable = b factor = 16
		//        #pragma HLS array_partition variable = Symbolic.Lnz factor = 16
		//        #pragma HLS array_partition variable = Numeric.Pnum
		//        #pragma HLS array_partition variable = Numeric.Offp
		//        #pragma HLS array_partition variable = Numeric.Offi
		//        #pragma HLS array_partition variable = Numeric.Lip
		//        #pragma HLS array_partition variable = Numeric.Uip
		//        #pragma HLS array_partition variable = Numeric.Llen
		//        #pragma HLS array_partition variable = Numeric.Ulen
		//        #pragma HLS array_partition variable = Numeric.Pinv
		//        #pragma HLS array_partition variable = Numeric.LUsize
		//        #pragma HLS array_partition variable = Numeric.Offx
		//        #pragma HLS array_partition variable = Numeric.Udiag
		//        #pragma HLS array_partition variable = Numeric.LUbx
		//        #pragma HLS array_partition variable = Numeric.Rs
//		        #pragma HLS array_partition variable = Numeric.Xwork

		Symbolic.n = N;
		Symbolic.nz = NZ;
		Symbolic.nzoff = NZOFF;
		Symbolic.nblocks = NBLOCKS;
		Symbolic.maxblock = MAXBLOCK;

		Numeric.n = Symbolic.n;
		Numeric.nblocks = Symbolic.nblocks;
		Numeric.nzoff = Symbolic.nzoff;
		klu_factor(Ap, Ai, Ax, &Symbolic, &Numeric, &Common);
		klu_solve(&Symbolic, &Numeric, N, nrhs, b, &Common);
	Write_loop:
		for (int i = 0; i < N * nrhs; i++)
		{
#pragma HLS pipeline
			B[i] = b[i];
		}
	}
}
