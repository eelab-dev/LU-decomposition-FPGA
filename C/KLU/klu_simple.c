/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */

#include "klu.h"
#include <stdio.h>
#include "../myKLU/mmio.c"

// int n = 5;
// int Ap[] = {0, 2, 5, 9, 10, 12};
// int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
// double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
// double b[] = {8., 45., -3., 3., 19.};

int main(void)
{
    Mtx mtx;
    char filename[] = "../../Matrix_Sample/bcsstk06.mtx";

    if (read_sparse(filename, &mtx))
        return 1;

    int n = mtx.n;
    int *Ap = mtx.Ap;
    int *Ai = mtx.Ai;
    double *Ax = mtx.Ax;
    // for (int i = 0; i < 100; i++)
    //     printf("Ai[%d]=%d\tAx[%d]=%.20lf\n", i, mtx.Ai[i], i, mtx.Ax[i]);
    double b[n];
    for (int i = 0; i < n; i++)
        b[i] = 1;
    klu_symbolic *Symbolic;
    klu_numeric *Numeric;
    klu_common Common;
    int i;
    klu_defaults(&Common);
    Symbolic = klu_analyze(n, Ap, Ai, &Common);
    Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);
    klu_solve(Symbolic, Numeric, n, 1, b, &Common);
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);
    for (i = 0; i < n; i++)
        printf("x [%d] = %g\n", i, b[i]);
    return (0);
}
