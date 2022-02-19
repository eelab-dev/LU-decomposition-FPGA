/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */

#include "klu.h"
#include <iostream>
#include "../myKLU/mmio.c"
#include <chrono>
#include <numeric>

// int n = 5;
// int Ap[] = {0, 2, 5, 9, 10, 12};
// int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
// double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
// double b[] = {8., 45., -3., 3., 19.};

int main(void)
{
    Mtx mtx;
    char filename[] = "../../Matrix_Sample/test.mtx";

    if (read_sparse(filename, &mtx))
        return 1;

    int n = mtx.n;
    int *Ap = mtx.Ap;
    int *Ai = mtx.Ai;
    double *Ax = mtx.Ax;
    double b[n];

    // for (int i = 0; i < 100; i++)
    //     printf("Ai[%d]=%d\tAx[%d]=%.20lf\n", i, mtx.Ai[i], i, mtx.Ax[i]);

    klu_common Common;
    klu_defaults(&Common);

    const int runtime = 10;
    std::chrono::steady_clock::time_point begin[3], end[3];
    long total[3] = {0};
    klu_symbolic *Symbolic;
    klu_numeric *Numeric;
    for (int i = 0; i < runtime; i++)
    {

        std::iota(b, b + n, 0);

        begin[0] = std::chrono::steady_clock::now();
        Symbolic = klu_analyze(n, Ap, Ai, &Common);
        end[0] = std::chrono::steady_clock::now();
        total[0] += std::chrono::duration_cast<std::chrono::microseconds>(end[0] - begin[0]).count();

        begin[1] = std::chrono::steady_clock::now();
        Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);
        end[1] = std::chrono::steady_clock::now();
        total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

        begin[2] = std::chrono::steady_clock::now();
        klu_solve(Symbolic, Numeric, n, 1, b, &Common);
        end[2] = std::chrono::steady_clock::now();
        total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();

        // klu_free_symbolic(&Symbolic, &Common);
        // klu_free_numeric(&Numeric, &Common);
    }

    for (int i = 0; i < n; i++)
        printf("x [%d] = %g\n", i, b[i]);

    std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
    return (0);
}
