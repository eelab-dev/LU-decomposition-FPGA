/* ========================================================================== */
/* === KLU_kernel =========================================================== */
/* ========================================================================== */

#include <iostream>
#include "klu_factor.h"
#include "klu_solve.h"
#include "mmio.h"
#include <chrono>
#include <numeric>
#include <fstream>

int main(void)
{
    std::vector<std::string> filename = {"rajat11.mtx", "rajat14.mtx", "rajat05.mtx", "oscil_dcop_01.mtx", "fpga_dcop_01.mtx"};
    std::string prefix = "../Matrix_Sample/Bench/";

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b;
    int n, nrhs = 10;
    const int runtime = 10;
    klu_common Common;
    std::ofstream data("Bench_Data.csv");

    data << "Matrix,";
    data << "Symbolic Average,Numeric Average,Solving Average (nrhs=" << nrhs << ")" << std::endl;

    for (int k = 0; k < filename.size(); k++)
    {
        if (read_sparse(prefix + filename[k], &n, Ap, Ai, Ax))
            return 1;

        data << filename[k] << ",";
        klu_defaults(&Common);

        klu_numeric Numeric;
        klu_symbolic Symbolic;
        Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);

        std::cout << "nblocks=" << Symbolic.nblocks << ",nzoff=" << Symbolic.nzoff << ",maxblock=" << Symbolic.maxblock << ",nnz=" << Symbolic.nz << std::endl;

        Numeric.Xwork = (double *)malloc(n * nrhs * sizeof(double));

        std::chrono::steady_clock::time_point begin[3], end[3];
        long total[3] = {0, 0, 0};

        b.resize(n * nrhs);

        for (int i = 0; i < runtime; i++)
        {
            std::iota(b.begin(), b.end(), 0);

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
            {
                if (j > 10)
                    break;
                printf("x[%d,%d] = %g\t", i, j, b[i + n * j]);
            }
            printf("x[%d,%d] = %g\n", i, nrhs - 1, b[i + n * (nrhs - 1)]);
            if (i > 10)
                break;
        }

        std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
        std::cout << "lusize_sum=" << Numeric.lusize_sum << std::endl;

        data << total[0] / (float)runtime << "," << total[1] / (float)runtime << "," << total[2] / (float)runtime << std::endl;
    }

    data.close();
    return 0;
}