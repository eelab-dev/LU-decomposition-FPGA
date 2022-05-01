#include <iostream>
#include <chrono>
#include <numeric>
#include <fstream>
#include <vector>
#include "umfpack.h"
#include "../myKLU/include/mmio.h"

int main(void)
{
    double *null = (double *)NULL;

    std::vector<std::string> filename = {"../../Matrix_Sample/host.mtx", "../../Matrix_Sample/rajat14.mtx"};

    std::vector<int> Ap, Ai;
    std::vector<double> Ax, b, x;
    int n;
    const int runtime = 1000;

    for (int k = 0; k < filename.size(); k++)
    {
        if (read_sparse(filename[k], &n, Ap, Ai, Ax))
            return 1;

        b.resize(n);
        x.resize(n);
        std::iota(b.begin(), b.end(), 0);

        std::chrono::steady_clock::time_point begin[3], end[3];
        long total[3] = {0};

        for (int i = 0; i < runtime; i++)
        {
            void *Symbolic, *Numeric;

            begin[0] = std::chrono::steady_clock::now();
            (void)umfpack_di_symbolic(n, n, Ap.data(), Ai.data(), Ax.data(), &Symbolic, null, null);
            end[0] = std::chrono::steady_clock::now();
            total[0] += std::chrono::duration_cast<std::chrono::microseconds>(end[0] - begin[0]).count();

            begin[1] = std::chrono::steady_clock::now();
            (void)umfpack_di_numeric(Ap.data(), Ai.data(), Ax.data(), Symbolic, &Numeric, null, null);
            end[1] = std::chrono::steady_clock::now();
            total[1] += std::chrono::duration_cast<std::chrono::microseconds>(end[1] - begin[1]).count();

            umfpack_di_free_symbolic(&Symbolic);

            begin[2] = std::chrono::steady_clock::now();
            (void)umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(), x.data(), b.data(), Numeric, null, null);
            end[2] = std::chrono::steady_clock::now();
            total[2] += std::chrono::duration_cast<std::chrono::microseconds>(end[2] - begin[2]).count();

            umfpack_di_free_numeric(&Numeric);
        }
        for (int i = 0; i < n; i++)
        {
            printf("x [%d] = %g\n", i, x[i]);
            if (i > 10)
                break;
        }

        std::cout << "Analyze time: " << total[0] / (float)runtime << "µs\nFactorization time: " << total[1] / (float)runtime << "µs\nSolving time: " << total[2] / (float)runtime << "µs" << std::endl;
    }

    return (0);
}
