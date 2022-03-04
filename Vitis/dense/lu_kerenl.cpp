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
#include <stdio.h>
#include <string.h>

#define MAX_SIZE 16

// TRIPCOUNT identifiers
const unsigned int c_dim = MAX_SIZE;

extern "C" {

void lu_kernel(int *aa, float *ll, float *uu, int size)
{
	float a[MAX_SIZE][MAX_SIZE], l[MAX_SIZE][MAX_SIZE], u[MAX_SIZE][MAX_SIZE];

#pragma HLS ARRAY_PARTITION variable=l factor=4 type=block dim=2
#pragma HLS ARRAY_PARTITION variable=u factor=4 type=block dim=2

	for (int itr = 0, i = 0, j = 0; itr < size * size; itr++, j++)
	{
	#pragma HLS LOOP_TRIPCOUNT min = c_dim* c_dim max = c_dim * c_dim
	        if (j == size) {
	            j = 0;
	            i++;
	        }
	        a[i][j] = aa[itr];
	}

	for (int i = 0; i < size; i++)
	    {
	        // Upper Triangular
	        for (int k = i; k < size; k++)
	        {
	            // Summation of L(i, j) * U(j, k)
	            float sum = 0;
	            for (int j = 0; j < i; j++)
	                sum += (l[i][j] * u[j][k]);

	            // Evaluating U(i, k)
	            u[i][k] = a[i][k] - sum;
	        }

	        // Lower Triangular
	        for (int k = i; k < size; k++)
	        {
	            if (i == k)
	                l[i][i] = 1; // Diagonal as 1
	            else
	            {
	                // Summation of L(k, j) * U(j, i)
	            	float sum = 0;
	                for (int j = 0; j < i; j++)
	                    sum += (l[k][j] * u[j][i]);

	                // Evaluating L(k, i)
	                l[k][i] = (a[k][i] - sum) / u[i][i];
	            }
	        }
	    }


	    for (int itr = 0, i = 0, j = 0; itr < size * size; itr++, j++) {
	#pragma HLS LOOP_TRIPCOUNT min = c_dim* c_dim max = c_dim * c_dim
	        if (j == size) {
	            j = 0;
	            i++;
	        }
	        ll[itr] = l[i][j];
	        uu[itr] = u[i][j];
	    }
	}
}
