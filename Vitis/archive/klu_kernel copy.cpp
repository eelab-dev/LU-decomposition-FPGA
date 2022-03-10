/**
 * Copyright (C) 2019-2021 Xilinx, Inc
 *
 * Licensed under the Apache License, Version 2.0 (the "License"). You may
 * not use this file except in compliance with the License. A copy of the
 * License is located at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

/*******************************************************************************
Description:

    This example uses the load/compute/store coding style which is generally
    the most efficient for implementing kernels using HLS. The load and store
    functions are responsible for moving data in and out of the kernel as
    efficiently as possible. The core functionality is decomposed across one
    of more compute functions. Whenever possible, the compute function should
    pass data through HLS streams and should contain a single set of nested loops.

    HLS stream objects are used to pass data between producer and consumer
    functions. Stream read and write operations have a blocking behavior which
    allows consumers and producers to synchronize with each other automatically.

    The dataflow pragma instructs the compiler to enable task-level pipelining.
    This is required for to load/compute/store functions to execute in a parallel
    and pipelined manner.

    The kernel operates on vectors of NUM_WORDS integers modeled using the hls::vector
    data type. This datatype provides intuitive support for parallelism and
    fits well the vector-add computation. The vector length is set to NUM_WORDS
    since NUM_WORDS integers amount to a total of 64 bytes, which is the maximum size of
    a kernel port. It is a good practice to match the compute bandwidth to the I/O
    bandwidth. Here the kernel loads, computes and stores NUM_WORDS integer values per
    clock cycle and is implemented as below:
                                       _____________
                                      |             |<----- Input Vector 1 from Global Memory
                                      |  load_input |       __
                                      |_____________|----->|  |
                                       _____________       |  | in1_stream
Input Vector 2 from Global Memory --->|             |      |__|
                               __     |  load_input |        |
                              |  |<---|_____________|        |
                   in2_stream |  |     _____________         |
                              |__|--->|             |<--------
                                      | compute_add |      __
                                      |_____________|---->|  |
                                       ______________     |  | out_stream
                                      |              |<---|__|
                                      | store_result |
                                      |______________|-----> Output result to Global Memory

*******************************************************************************/

// Includes
#include <hls_vector.h>
#include <hls_stream.h>
#include "assert.h"
#include "klu_factor.h"
#include "klu_solve.h"

extern "C"
{
    void lu(int *AP, int *AI, double *AX, int *P, int *Q, int *R, double *LNZ, int N, int NBLOCKS, int MAXBLOCK, int NZOFF, int NZ, double *B)
    {
        klu_common Common;
        KLU_numeric Numeric;
        klu_symbolic Symbolic;
        klu_defaults(&Common);

#pragma HLS INTERFACE m_axi depth=2048 port=AP max_read_burst_length=32 offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi depth=2048 port=AI max_read_burst_length=32 offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi depth=2048 port=AX max_read_burst_length=32 offset=slave bundle=gmem2
#pragma HLS INTERFACE m_axi depth=2048 port=P max_read_burst_length=32 offset=slave bundle=gmem3
#pragma HLS INTERFACE m_axi depth=2048 port=Q max_read_burst_length=32 offset=slave bundle=gmem4
#pragma HLS INTERFACE m_axi depth=2048 port=R max_read_burst_length=32 offset=slave bundle=gmem5
#pragma HLS INTERFACE m_axi depth=2048 port=LNZ max_read_burst_length=32 offset=slave bundle=gmem6
#pragma HLS INTERFACE m_axi depth=2048 port=B max_read_burst_length=32 offset=slave bundle=gmem7


        int Ap[MAX_SIZE], Ai[MAX_SIZE], p[MAX_SIZE], q[MAX_SIZE], r[MAX_SIZE];
        double Ax[MAX_SIZE], b[MAX_SIZE], Lnz[MAX_SIZE];

        int Pnum[MAX_SIZE], Offp[MAX_SIZE], Offi[MAX_SIZE], Lip[MAX_SIZE], Uip[MAX_SIZE], Llen[MAX_SIZE], Ulen[MAX_SIZE], Pinv[MAX_SIZE], LUsize[MAX_SIZE];

        double Offx[MAX_SIZE], Udiag[MAX_SIZE], LUbx[MAX_SIZE], Rs[MAX_SIZE], Xwork[MAX_SIZE];

//#pragma HLS array_partition variable = Ap factor = 16
//#pragma HLS array_partition variable = Ai factor = 16
//#pragma HLS array_partition variable = p
//#pragma HLS array_partition variable = q
//#pragma HLS array_partition variable = r
//#pragma HLS array_partition variable = Ax
//#pragma HLS array_partition variable = b factor = 16
//#pragma HLS array_partition variable = Lnz
//#pragma HLS array_partition variable = Pnum
//#pragma HLS array_partition variable = Offp
//#pragma HLS array_partition variable = Offi
//#pragma HLS array_partition variable = Lip
//#pragma HLS array_partition variable = Uip
//#pragma HLS array_partition variable = Llen
//#pragma HLS array_partition variable = Ulen
//#pragma HLS array_partition variable = Pinv
//#pragma HLS array_partition variable = LUsize
//#pragma HLS array_partition variable = Offx
//#pragma HLS array_partition variable = Udiag
//#pragma HLS array_partition variable = LUbx
//#pragma HLS array_partition variable = Rs
//#pragma HLS array_partition variable = Xwork

        for (int i = 0; i < N; i++)
        {
#pragma HLS pipeline
            Ap[i] = AP[i];
            p[i] = P[i];
            q[i] = Q[i];
            r[i] = R[i];
            b[i] = B[i];
            Lnz[i] = LNZ[i];
            Ai[i] = AI[i];
            Ax[i] = AX[i];
        }
        Ap[N] = AP[N];
        r[N] = R[N];
        for (int i = N; i < NZ; i++)
        {
#pragma HLS pipeline
            Ai[i] = AI[i];
            Ax[i] = AX[i];
        }

        Symbolic.n = N;
        Symbolic.nz = NZ;
        Symbolic.P = p;
        Symbolic.Q = q;
        Symbolic.R = r;
        Symbolic.nzoff = NZOFF;
        Symbolic.nblocks = NBLOCKS;
        Symbolic.maxblock = MAXBLOCK;
        Symbolic.Lnz = Lnz;

        Numeric.n = Symbolic.n;
        Numeric.nblocks = Symbolic.nblocks;
        Numeric.nzoff = Symbolic.nzoff;
        Numeric.Pnum = Pnum;
        Numeric.Offp = Offp;
        Numeric.Offi = Offi;
        Numeric.Offx = Offx;
        Numeric.Lip = Lip;
        Numeric.Uip = Uip;
        Numeric.Llen = Llen;
        Numeric.Ulen = Ulen;
        Numeric.LUsize = LUsize;
        Numeric.LUbx = LUbx;
        Numeric.Udiag = Udiag;
        Numeric.Rs = Rs;
        Numeric.Pinv = Pinv;
        Numeric.Xwork = Xwork;

        klu_factor(Ap, Ai, Ax, &Symbolic, &Numeric, &Common);

        // klu_solve(Symbolic.Q, Symbolic.R, Numeric.Pnum, Numeric.Offp, Numeric.Offi, Numeric.Lip, Numeric.Uip, Numeric.Llen, Numeric.Ulen, Numeric.Offx, Numeric.Xwork, Numeric.Udiag, Numeric.Rs, Numeric.LUbx, Numeric.LUsize, Symbolic.nblocks, N, Numeric.lusize_sum, b, &Common);

        klu_solve(&Symbolic, &Numeric, N, b, &Common);

        for (int i = 0; i < N; i++)
        {
#pragma HLS pipeline
            B[i] = b[i];
        }
    }
}
