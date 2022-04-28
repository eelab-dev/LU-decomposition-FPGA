#include "klu.h"
#include "cholmod.h"
#include <stdio.h>

// int n = 5;
// int Ap[] = {0, 2, 5, 9, 10, 12};
// int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
// double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
// double b[] = {8., 45., -3., 3., 19.};
int n = 8;
int Ap[] = {0, 2, 4, 6, 9, 13, 15, 17, 20};
int Ai[] = {0, 2, 0, 1, 1, 2, 1, 3, 4, 0, 4, 6, 7, 4, 5, 5, 6, 0, 6, 7};
double Ax[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
double b[] = {1, 1, 1, 1, 1, 1, 1, 1};

int main(void)
{
    // cholmod_sparse *A;
    // cholmod_common ch;
    // cholmod_start(&ch);
    // A = cholmod_read_sparse(stdin, &ch);

    // n = A->nrow;

    int P[n], Q[n], R[n + 1], Work[5 * n], ncomp, nfound;
    double maxwork, work;

    // ncomp = btf_order(n, A->p, A->i, maxwork, &work, P, Q, R, &nfound, Work);
    ncomp = btf_order(n, Ap, Ai, maxwork, &work, P, Q, R, &nfound, Work);

    printf("ncomp=%d\nnfound=%d\nmaxwork=%lf\nwork=%lf\n", ncomp, nfound, maxwork, work);
    for (int i = 0; i < n; i++)
    {
        printf("P[%d]=%d\tQ[%d]=%d\tR[%d]=%d\n", i, P[i], i, Q[i], i, R[i]);
    }
    printf("R[%d]=%d\n", n, R[n]);
    for (int i = 0; i < 5 * n; i++)
        printf("Work[%d]=%d\n", i, Work[i]);

    // cholmod_free_sparse(&A, &ch);
    // cholmod_finish(&ch);
    return 0;
}