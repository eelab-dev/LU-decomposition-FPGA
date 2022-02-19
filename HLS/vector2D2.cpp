// Unit Lower matrix

#include <hls_vector.h>
#include "ap_int.h"

const int SIZE = 30;

typedef ap_uint<10> sizet;
typedef hls::vector<hls::vector<float,SIZE>,SIZE> mat;

typedef struct{
	float l[SIZE][SIZE];
	float u[SIZE][SIZE];
} matrix_lu;

matrix_lu LUdecomposition(float a[SIZE][SIZE])//, float l[SIZE][SIZE], float u[SIZE][SIZE])
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix

	int n = SIZE;
	matrix_lu lu;

#pragma HLS ARRAY_PARTITION variable=lu.l factor=10 type=block dim=2
#pragma HLS ARRAY_PARTITION variable=lu.u factor=10 type=block dim=2

    for (sizet i = 0; i < n; i++)
    {
        // Upper Triangular
        for (sizet k = i; k < n; k++)
        {
            // Summation of L(i, j) * U(j, k)
            float sum = 0;
            for (sizet j = 0; j < i; j++)
            {
//				#pragma HLS unroll factor=2
                sum += (lu.l[i][j] * lu.u[j][k]);
            }

            // Evaluating U(i, k)
            lu.u[i][k] = a[i][k] - sum;
        }

        // Lower Triangular
        for (sizet k = i; k < n; k++)
        {
            if (i == k)
                lu.l[i][i] = 1; // Diagonal as 1
            else
            {
                // Summation of L(k, j) * U(j, i)
                float sum = 0;
                for (sizet j = 0; j < i; j++)
                {
//					#pragma HLS unroll factor=2
                    sum += (lu.l[k][j] * lu.u[j][i]);
                }

                // Evaluating L(k, i)
                lu.l[k][i] = (a[k][i] - sum) / lu.u[i][i];
            }
        }
    }

    return lu;
}
