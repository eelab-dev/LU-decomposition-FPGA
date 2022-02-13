#include "klu.h"
#include <iostream>
#include <chrono>
#include "cholmod.h"

int main(void)
{
    cholmod_sparse *A;
    cholmod_common ch;
    cholmod_start(&ch);
    A = cholmod_read_sparse(stdin, &ch);
    int solved;
    if (A)
    {
        if (A->nrow != A->ncol || A->stype != 1 || (!(A->xtype == CHOLMOD_REAL || A->xtype == CHOLMOD_COMPLEX)))
        {
            printf("stype=%d,xtype=%d", A->stype, A->xtype);
            printf("invalid matrix\n");
            return 0;
        }
        std::cout << "Loaded" << std::endl;
        klu_common Common;
        klu_defaults(&Common);
        // Common.btf = 0;

        int runtime = 1;
        std::chrono::steady_clock::time_point begin, end;
        long total = 0;

        int *Ap, *Ai;
        double *Ax;

        klu_numeric *Numeric;
        klu_symbolic *Symbolic;
        Symbolic = klu_analyze(A->nrow, Ap, Ai, &Common);

        for (int i = 0; i < runtime; i++)
        {
            Ap = (int *)(A->p);
            Ai = (int *)(A->i);
            Ax = (double *)(A->x);

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
        double b[] = {42, 61, 27, 39, 19, 60, 21, 29, 56, 66, 39, 36, 49, 46, 45, 112, 78, 42, 88, 73, 59, 101, 105, 63, 90, 135, 70, 83, 121, 32, 37, 42, 52, 54, 57, 59, 62, 67, 49};
        solved = klu_solve(Symbolic, Numeric, A->nrow, 1, b, &Common);
        // // klu_demo (A->nrow, A->p, A->i, A->x, A->xtype == CHOLMOD_REAL) ;
        // if (solved)

        //     for (int i = 0; i < A->nrow; i++)
        //     {
        //         if (b[i] != 2)
        //             printf("x [%d] = %g\n", i, b[i]);
        //     }
        // else
        //     printf("Invalid");
        // }
        cholmod_free_sparse(&A, &ch);
    }
    cholmod_finish(&ch);

    return (0);
}