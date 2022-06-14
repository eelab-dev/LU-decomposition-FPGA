/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */

#include "klu.h"
#include <iostream>
#include "../myKLU/include/mmio.h"
#include <chrono>
#include <numeric>

// int n = 5;
// int Ap[] = {0, 2, 5, 9, 10, 12};
// int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
// double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
// double b[] = {8., 45., -3., 3., 19.};

int main(void)
{
    char filename[] = "../Matrix_Sample/host.mtx";
    char bmatrix[] = "../Matrix_Sample/host_b.mtx";

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b;
    int n;
    if (read_sparse(filename, &n, Ap, Ai, Ax))
        return 1;

    klu_common Common;
    klu_defaults(&Common);

    const int runtime = 1000;
    std::chrono::steady_clock::time_point begin[3], end[3];
    long total[3] = {0};
    klu_symbolic *Symbolic;
    klu_numeric *Numeric;
    int nrhs;
    for (int i = 0; i < runtime; i++)
    {
        read_bmatrix(bmatrix, b, &nrhs);
        // std::iota(b.begin(), b.end(), 10);

        begin[0] = std::chrono::steady_clock::now();
        Symbolic = klu_analyze(n, Ap.data(), Ai.data(), &Common);
        end[0] = std::chrono::steady_clock::now();
        total[0] += std::chrono::duration_cast<std::chrono::microseconds>(end[0] - begin[0]).count();

        begin[1] = std::chrono::steady_clock::now();
        Numeric = klu_factor(Ap.data(), Ai.data(), Ax.data(), Symbolic, &Common);
        end[1] = std::chrono::steady_clock::now();
        total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

        begin[2] = std::chrono::steady_clock::now();
        klu_solve(Symbolic, Numeric, n, nrhs, b.data(), &Common);
        end[2] = std::chrono::steady_clock::now();
        total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();

        // klu_free_symbolic(&Symbolic, &Common);
        // klu_free_numeric(&Numeric, &Common);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs - 1; j++)
        {
            printf("x [%d,%d] = %g\t", i, j, b[i + n * j]);
            if (j > 64)
                break;
        }

        printf("x [%d,%d] = %g\n", i, nrhs - 1, b[i + n * (nrhs - 1)]);
        if (i > 10)
            break;
    }

    std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
    return (0);
}
