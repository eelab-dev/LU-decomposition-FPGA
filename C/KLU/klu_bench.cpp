/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */

#include "klu.h"
#include <iostream>
#include "../myKLU/include/mmio.h"
#include <chrono>
#include <numeric>
#include <fstream>

int main(void)
{
    std::vector<std::string> filename = {"rajat11.mtx", "rajat14.mtx"};
    // std::vector<std::string> bmatrix = {"../../Matrix_Sample/host_b.mtx"};
    std::string prefix = "../../Matrix_Sample/Bench/";

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b;
    int n, nrhs = 1;
    const int runtime = 1000;
    klu_common Common;
    std::ofstream data("Bench_KLU_Data.csv");

    data << "Matrix,";
    data << "Symbolic Average,Numeric Average,Solving Average (nrhs=" << nrhs << ")" << std::endl;

    for (int k = 0; k < filename.size(); k++)
    {
        if (read_sparse(prefix + filename[k], &n, Ap, Ai, Ax))
            return 1;

        data << filename[k] << ",";
        klu_defaults(&Common);

        std::chrono::steady_clock::time_point begin[3], end[3];
        long total[3] = {0};
        klu_symbolic *Symbolic;
        klu_numeric *Numeric;
        b.resize(n * nrhs);

        for (int i = 0; i < runtime; i++)
        {
            // read_bmatrix(bmatrix[k], b, &nrhs);
            std::iota(b.begin(), b.end(), 0);

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
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nrhs - 1; j++)
                printf("x [%d,%d] = %g\t", i, j, b[i + n * j]);
            printf("x [%d,%d] = %g\n", i, nrhs - 1, b[i + n * (nrhs - 1)]);
            if (i > 10)
                break;
        }

        std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;

        data << total[0] / (float)runtime << "," << total[1] / (float)runtime << "," << total[2] / (float)runtime << std::endl;
        klu_free_symbolic(&Symbolic, &Common);
        klu_free_numeric(&Numeric, &Common);
    }

    data.close();
    return (0);
}
