/* ========================================================================== */
/* === KLU_kernel =========================================================== */
/* ========================================================================== */

#include <iostream>
#include "klu_factor.h"
#include "klu_solve.h"
#include "mmio.h"
#include <chrono>
#include <numeric>

const int runtime = 10;

int main(void)
{
    std::string filename, bmatrix;
    std::cout << "Left matrix file path (default - " << "./host.mtx): ";
    std::getline(std::cin, filename);
    if (filename.empty())
        filename = "./host.mtx";

    std::cout << "B matrix file path (default - " << "./host_b.mtx): ";
    std::getline(std::cin, bmatrix);
    if (bmatrix.empty())
        bmatrix = "./host_b.mtx";

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b;
    int n, nrhs;
    if (read_sparse(filename, &n, Ap, Ai, Ax))
        return 1;

    read_bmatrix(bmatrix, b, &nrhs);

    klu_common Common;
    KLU_numeric Numeric;
    klu_symbolic Symbolic;
    klu_defaults(&Common);
    Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);

    std::cout << "nblocks=" << Symbolic.nblocks << ",nzoff=" << Symbolic.nzoff << ",maxblock=" << Symbolic.maxblock << ",nnz=" << Symbolic.nz << std::endl;

    Numeric.Xwork = (double *)malloc(n * nrhs * sizeof(double));

    std::chrono::steady_clock::time_point begin[3], end[3];
    long total[3] = {0, 0, 0};

    for (int i = 0; i < runtime; i++)
    {
        read_bmatrix(bmatrix, b, &nrhs);

        begin[1] = std::chrono::steady_clock::now();
        klu_factor(Ap.data(), Ai.data(), Ax.data(), &Symbolic, &Numeric, &Common);
        end[1] = std::chrono::steady_clock::now();
        total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

        begin[2] = std::chrono::steady_clock::now();
        klu_solve(&Symbolic, &Numeric, n, nrhs, b.data(), &Common);
        end[2] = std::chrono::steady_clock::now();
        total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs - 1; j++)
            std::cout << "x[" << i << "," << j << "] = " << b[i + n * j] << "\t";
        std::cout << "x[" << i << "," << nrhs - 1 << "] = " << b[i + n * (nrhs - 1)] << std::endl;
        if (i > 10)
            break;
    }

    std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
    std::cout << "nnz: " << Numeric.lusize_sum << std::endl;
    return 0;
}