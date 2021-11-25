/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */

#include "klu.h"
#include <iostream>
#include <stdio.h>

int n = 5;
int Ap[] = {0, 2, 5, 9, 10, 12};
int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
double b[] = {8., 45., -3., 3., 19.};

int main(void)
{
    klu_symbolic *Symbolic;
    klu_numeric *Numeric;
    klu_common Common;
    int i;

    klu_defaults(&Common);
    Symbolic = klu_analyze(n, Ap, Ai, &Common);
    Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);

    int lnz = Numeric->lnz, unz = Numeric->unz, nzoff = Numeric->nzoff;
    int ok, Lp[n + 1], Li[lnz], Up[n + 1], Ui[unz], Fp[n + 1], Fi[nzoff], P[n], Q[n], R[n];
    double Lx[lnz], Ux[unz], Fx[nzoff], Rs[n];

    ok = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, &Common);

    for (i = 0; i < n + 1; i++)
        printf("Lp[%d]=%d\n", i, Lp[i]);
    for (i = 0; i < lnz; i++)
        printf("Li[%d]=%d\n", i, Li[i]);
    for (i = 0; i < lnz; i++)
        printf("Lx[%d]=%.4lf\n", i, Lx[i]);

    for (i = 0; i < n + 1; i++)
        printf("Up[%d]=%d\n", i, Up[i]);
    for (i = 0; i < unz; i++)
        printf("Ui[%d]=%d\n", i, Ui[i]);
    for (i = 0; i < unz; i++)
        printf("Ux[%d]=%.4lf\n", i, Ux[i]);

    for (i = 0; i < n + 1; i++)
        printf("Fp[%d]=%d\n", i, Fp[i]);
    for (i = 0; i < nzoff; i++)
        printf("Fi[%d]=%d\n", i, Fi[i]);
    for (i = 0; i < nzoff; i++)
        printf("Fx[%d]=%.4lf\n", i, Fx[i]);

    for (int i = 0; i < n; i++)
        printf("P[%d]=%d\tQ[%d]=%d\tR[%d]=%d\tRs[%d]=%.4lf\n", i, P[i], i, Q[i], i, R[i], i, Rs[i]);

    klu_solve(Symbolic, Numeric, 5, 1, b, &Common);
    klu_free_symbolic(&Symbolic, &Common);
    klu_free_numeric(&Numeric, &Common);
    for (i = 0; i < n; i++)
        printf("x [%d] = %g\n", i, b[i]);
    return (0);
}
