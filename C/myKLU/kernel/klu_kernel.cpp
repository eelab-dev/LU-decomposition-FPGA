/* ========================================================================== */
/* === KLU_kernel =========================================================== */
/* ========================================================================== */

#include <iostream>
#include "klu.h"
#include "klu_factor.h"
#include "klu_solve.h"
#include "../mmio.h"
#include <chrono>
#include <numeric>

int main(void)
{
    char filename[] = "../../../Matrix_Sample/host.mtx";
    char bmatrix[] = "../../../Matrix_Sample/host_b.mtx";

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b;
    int n;
    if (read_sparse(filename, &n, Ap, Ai, Ax))
        return 1;

    // b.resize(2 * n);

    klu_common Common;
    KLU_numeric Numeric;
    klu_symbolic Symbolic;
    klu_defaults(&Common);
    // const int n = 10;
    // int Ap[] = {0, 2, 4, 7, 9, 13, 14, 18, 21, 23, 24};
    // int Ai[] = {0, 7, 1, 5, 0, 2, 7, 3, 8, 3, 4, 5, 8, 5, 0, 6, 7, 8, 4, 5, 7, 2, 8, 9};
    // double Ax[] = {8, 8, 2, 4, 10, 3, 10, 1, 5, 2, 1, 4, 10, 2, 3, 5, 3, 15, 4, 16, 3, 7, 2, 9};
    // double b2[] = {172, 18, 38, 19, 18, 118, 20, 181, 159, 9};
    Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);

    // for (int i = 0; i < n; i++)
    //     printf("P[%d]=%d,Q[%d]=%d,R[%d]=%d,Lnz[%d]=%lf\n", i, Symbolic.P[i], i, Symbolic.Q[i], i, Symbolic.R[i], i, Symbolic.Lnz[i]);
    printf("nblocks=%d,nzoff=%d,maxblock=%d\n", Symbolic.nblocks, Symbolic.nzoff, Symbolic.maxblock);

    int nzoff1 = Symbolic.nzoff + 1, n1 = n + 1;
    double lusize = Common.memgrow * (Symbolic.lnz + Symbolic.unz) + 4 * n + 1;

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
    Numeric.Xwork = (double *)calloc(n, sizeof(double));
    // Numeric.Iwork = (int *)malloc(sizeof(int) * 6 * Symbolic.maxblock);

    const int runtime = 10;
    std::chrono::steady_clock::time_point begin[3], end[3];
    long total[3] = {0};

    for (int i = 0; i < runtime; i++)
    {
        // std::iota(b.begin(), b.end(), 0);
        read_bmatrix(bmatrix, b);

        begin[1] = std::chrono::steady_clock::now();
        factor2(Ap.data(), Ai.data(), Ax.data(), Symbolic.n, Symbolic.nblocks, Symbolic.nzoff, Symbolic.P, Symbolic.Q, Symbolic.R, Symbolic.Lnz, Numeric.Pnum, Numeric.Offp, Numeric.Offi, Numeric.Lip, Numeric.Uip, Numeric.Llen, Numeric.Ulen, Numeric.Offx, Numeric.Udiag, Numeric.Rs, &Numeric.lusize_sum, Numeric.LUbx, Numeric.Pinv, Numeric.Xwork, Numeric.LUsize, &Numeric.lnz, &Numeric.unz, &Numeric.max_lnz_block, &Numeric.max_unz_block, &Common);
        end[1] = std::chrono::steady_clock::now();
        total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

        begin[2] = std::chrono::steady_clock::now();
        klu_solve2(Symbolic.Q, Symbolic.R, Numeric.Pnum, Numeric.Offp, Numeric.Offi, Numeric.Lip, Numeric.Uip, Numeric.Llen, Numeric.Ulen, Numeric.Offx, Numeric.Xwork, Numeric.Udiag, Numeric.Rs, Numeric.LUbx, Numeric.LUsize, Symbolic.nblocks, n, Numeric.lusize_sum, b.data(), &Common);
        end[2] = std::chrono::steady_clock::now();
        total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();
    }

    for (int i = 0; i < n; i++)
        printf("x [%d] = %g\n", i, b[i]);

    std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
    printf("lusize=%d\n", Numeric.lusize_sum);
    return 0;
}