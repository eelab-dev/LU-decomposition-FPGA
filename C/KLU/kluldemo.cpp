#include "klu.h"
#include <iostream>
#include <chrono>
#include "cholmod.h"
#include <numeric>

int main(void)
{
    cholmod_sparse *A;
    cholmod_common ch;
    ch.prefer_upper = false;
    cholmod_start(&ch);
    A = cholmod_read_sparse(stdin, &ch);
    int solved;

    for (int i = 0; i < A->ncol; i++)
        std::cout << "Ap[" << i << "]=" << ((int *)(A->p))[i] << std::endl;
    for (int i = 0; i < A->nzmax; i++)
        std::cout << "Ai[" << i << "]=" << ((int *)(A->i))[i] << "\tAx[" << i << "]=" << ((double *)(A->x))[i] << std::endl;
    /*
    if (A)
    {
        // if (A->nrow != A->ncol || A->stype != 0 || (!(A->xtype == CHOLMOD_REAL || A->xtype == CHOLMOD_COMPLEX)))
        // {
        //     printf("stype=%d,xtype=%d", A->stype, A->xtype);
        //     printf("invalid matrix\n");
        //     return 0;
        // }
        std::cout << "Loaded" << std::endl;
        klu_common Common;
        klu_defaults(&Common);
        // Common.btf = 0;

        int runtime = 1;
        std::chrono::steady_clock::time_point begin, end;
        long total = 0;

        int *Ap, *Ai;
        double *Ax;
        Ap = (int *)(A->p);
        Ai = (int *)(A->i);
        Ax = (double *)(A->x);

        klu_numeric *Numeric;
        klu_symbolic *Symbolic;
        Symbolic = klu_analyze(A->nrow, Ap, Ai, &Common);

        for (int i = 0; i < runtime; i++)
        {

            begin = std::chrono::steady_clock::now();
            Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);
            end = std::chrono::steady_clock::now();

            total += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

            std::cout << "Maxblock: " << Symbolic->maxblock << std::endl;
            std::cout << "nblocks: " << Symbolic->nblocks << std::endl;

            // klu_free_symbolic(&Symbolic, &Common);
            // klu_free_numeric(&Numeric, &Common);
        }

        printf("Time:%lfns\n", total / (double)runtime);
        // double *b;
        // b = klu_malloc(A->nrow, sizeof(double), &Common);
        // for (int i = 0; i < A->nrow; i++)
        //     b[i] = 3;
        double b[A->nrow];
        std::iota(b, b + A->nrow, 0);
        solved = klu_solve(Symbolic, Numeric, A->nrow, 1, b, &Common);
        // // klu_demo (A->nrow, A->p, A->i, A->x, A->xtype == CHOLMOD_REAL) ;
        for (int i = 0; i < A->nrow; i++)
            printf("x [%d] = %g\n", i, b[i]);
        cholmod_free_sparse(&A, &ch);
    }*/
    cholmod_finish(&ch);

    return (0);
}