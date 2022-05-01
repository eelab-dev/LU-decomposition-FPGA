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
    std::vector<std::string> filename = {"../../Matrix_Sample/host.mtx", "../../Matrix_Sample/rajat14.mtx"};

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b;
    int n, nrhs = 1000;
    const int runtime = 1000;
    klu_common Common;
    std::ofstream data("Bench_Data.csv");

    data << "Matrix,";
    for (int i = 0; i < runtime; i++)
        data << i << "," << i << "," << i << ",";
    data << "Average,Average,Average,std,std,std" << std::endl;

    for (int k = 0; k < filename.size(); k++)
    {
        if (read_sparse(filename[k], &n, Ap, Ai, Ax))
            return 1;

        data << filename[k] << ",";
        klu_defaults(&Common);

        klu_numeric Numeric;
        klu_symbolic Symbolic;
        Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);

        printf("nblocks=%d,nzoff=%d,maxblock=%d,nnz=%d\n", Symbolic.nblocks, Symbolic.nzoff, Symbolic.maxblock, Symbolic.nz);

        int nzoff1 = Symbolic.nzoff + 1, n1 = n + 1;
        int lusize = n * n;

        Numeric.n = Symbolic.n;
        Numeric.nblocks = Symbolic.nblocks;
        Numeric.nzoff = Symbolic.nzoff;
        Numeric.Pnum = (int *)malloc(n * sizeof(int));
        Numeric.Offp = (int *)malloc(n1 * sizeof(int));
        Numeric.Offi = (int *)malloc(nzoff1 * sizeof(int));
        Numeric.Offx = (double *)malloc(nzoff1 * sizeof(double));
        Numeric.Lip = (int *)calloc(n, sizeof(int));
        Numeric.Uip = (int *)malloc(n * sizeof(int));
        Numeric.Llen = (int *)malloc(n * sizeof(int));
        Numeric.Ulen = (int *)malloc(n * sizeof(int));
        Numeric.LUsize = (int *)calloc(Symbolic.nblocks, sizeof(int));
        Numeric.LUbx = (double *)calloc(lusize * 2, sizeof(double));
        Numeric.Udiag = (double *)malloc(n * sizeof(double));
        Numeric.Rs = (double *)malloc(n * sizeof(double));
        Numeric.Pinv = (int *)malloc(n * sizeof(int));
        Numeric.worksize = n * sizeof(double) + MAX(n * 3 * sizeof(double), Symbolic.maxblock * 6 * sizeof(int));
        Numeric.Xwork = (double *)calloc(n * nrhs, sizeof(double));

        std::chrono::steady_clock::time_point begin[3], end[3];
        long total[3] = {0};

        b.resize(n * nrhs);

        for (int i = 0; i < runtime; i++)
        {
            std::iota(b.begin(), b.end(), 0);
            // read_bmatrix(bmatrix, b, &nrhs);

            begin[1] = std::chrono::steady_clock::now();
            klu_factor(Ap.data(), Ai.data(), Ax.data(), &Symbolic, &Numeric, &Common);
            end[1] = std::chrono::steady_clock::now();
            total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

            begin[2] = std::chrono::steady_clock::now();
            klu_solve(&Symbolic, &Numeric, n, nrhs, b.data(), &Common);
            end[2] = std::chrono::steady_clock::now();
            total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();

            data << (end[0] - begin[0]).count() << "," << (end[1] - begin[1]).count() << "," << (end[2] - begin[2]).count() << ",";
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nrhs - 1; j++)
                printf("x[%d,%d] = %g\t", i, j, b[i + n * j]);
            printf("x[%d,%d] = %g\n", i, nrhs - 1, b[i + n * (nrhs - 1)]);
            if (i > 10)
                break;
        }

        std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
        printf("lusize=%d\n", Numeric.lusize_sum);

        data << total[0] / (float)runtime << "," << total[1] / (float)runtime << "," << total[2] / (float)runtime << "," << std::endl;
    }

    data.close();
    return 0;
}