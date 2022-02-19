/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) 2005-2012 by Timothy A. Davis,                       */
/* http://www.suitesparse.com. All Rights Reserved.                           */
/* See ../Doc/License.txt for License.                                        */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include "../myKLU/mmio.c"
#include <chrono>
#include "umfpack.h"
#include <numeric>

// int n = 5;
// int Ap[] = {0, 2, 5, 9, 10, 12};
// int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
// double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
// double b[] = {8., 45., -3., 3., 19.};
// double x[5];

int main(void)
{
    Mtx mtx;
    char filename[] = "../../Matrix_Sample/bcircuit.mtx";

    if (read_sparse(filename, &mtx))
        return 1;

    int n = mtx.n;
    int *Ap = mtx.Ap;
    int *Ai = mtx.Ai;
    double *Ax = mtx.Ax;
    double b[n];
    double x[n];

    const int runtime = 10;
    std::chrono::steady_clock::time_point begin[3], end[3];
    long total[3] = {0};

    double *null = (double *)NULL;
    for (int i = 0; i < runtime; i++)
    {

        std::iota(b, b + n, 0);
        void *Symbolic, *Numeric;

        begin[0] = std::chrono::steady_clock::now();
        (void)umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, null, null);
        end[0] = std::chrono::steady_clock::now();
        total[0] += std::chrono::duration_cast<std::chrono::microseconds>(end[0] - begin[0]).count();

        begin[1] = std::chrono::steady_clock::now();
        (void)umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
        end[1] = std::chrono::steady_clock::now();
        total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

        umfpack_di_free_symbolic(&Symbolic);

        begin[2] = std::chrono::steady_clock::now();
        (void)umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
        end[2] = std::chrono::steady_clock::now();
        total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();

        umfpack_di_free_numeric(&Numeric);
    }

    for (int i = 0; i < n; i++)
        printf("x [%d] = %g\n", i, x[i]);

    std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
    return (0);
}
