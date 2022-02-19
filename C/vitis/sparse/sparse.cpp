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

#define MAX_SIZE 30

extern "C"
{

    /*
    Vector Addition Kernel

    Arguments:
        in1  (input)  --> Input vector 1
        in2  (input)  --> Input vector 2
        out  (output) --> Output vector
        size (input)  --> Number of elements in vector
*/

    void sparselu(int squareSize, int *AP, int *AI, double *AX, int *LP, int *LI, double *LX, int *UP, int *UI, double *UX)
    {
        int Ap[MAX_SIZE], Ai[MAX_SIZE], Lp[MAX_SIZE], Li[MAX_SIZE], Up[MAX_SIZE], Ui[MAX_SIZE];
        double Ax[MAX_SIZE], Lx[MAX_SIZE], Ux[MAX_SIZE];

        // #pragma HLS ARRAY_PARTITION variable = Ap factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Ai factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Ax factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Li factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Lp factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Li factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Ux factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Ui factor = 4 type = block
        // #pragma HLS ARRAY_PARTITION variable = Ux factor = 4 type = block

        for (int i = 0; i < squareSize; i++)
            Ap[i] = AP[i];
        for (int i = 0; i < AP[squareSize]; i++)
        {
            Ai[i] = AI[i];
            Ax[i] = AX[i];
        }

        for (int i = 0; i < squareSize; i++)
        {
            double b[MAX_SIZE] = {0}, x[MAX_SIZE] = {0};

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
                double sum = 0;

                for (int k = 0; k <= j; k++)
                {
                    for (int t = Lp[k]; t < Lp[k + 1]; t++)
                    {
                        if (Li[t] == j)
                            sum += Lx[t] * x[k];
                        else if (Li[t] > j)
                            break;
                    }
                }
                x[j] = b[j] - sum;
            }

            for (int j = Up[i]; j < Up[i + 1]; j++)
            {
                if ((x[Ui[j]] > 1e-6) || (x[Ui[j]] < -1e6))
                    Ux[j] = x[Ui[j]];
            }

            for (int j = Lp[i] + 1; j < Lp[i + 1]; j++)
            {
                if ((x[Ui[j]] > 1e-6) || (x[Ui[j]] < -1e6))
                    Lx[j] = x[Li[j]] / Ux[Up[i + 1] - 1];
            }
        }

        for (int i = 0; i <= squareSize; i++)
        {
            LP[i] = Lp[i];
            UP[i] = Up[i];
        }
        for (int i = 0; i < Up[squareSize]; i++)
        {
            UI[i] = Ui[i];
            UX[i] = Ux[i];
            LI[i] = Li[i];
            LX[i] = Lx[i];
        }
    }
}
