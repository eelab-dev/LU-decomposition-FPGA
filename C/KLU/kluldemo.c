#include <math.h>
#include <stdio.h>
#include "klu.h"
#include <time.h>
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
        // if (A->nrow != A->ncol || A->stype != 0 || (!(A->xtype == CHOLMOD_REAL || A->xtype == CHOLMOD_COMPLEX)))
        // {
        //     printf("stype=%d,xtype=%d", A->stype, A->xtype);
        //     printf("invalid matrix\n");
        // }
        // else
        // {

        klu_common Common;
        klu_defaults(&Common);

        int runtime = 10;
        clock_t begin, end, total = 0;

        for (int i = 0; i < runtime; i++)
        {
            klu_numeric *Numeric;
            klu_symbolic *Symbolic;
            Symbolic = klu_analyze(A->nrow, A->p, A->i, &Common);

            begin = clock();
            Numeric = klu_factor(A->p, A->i, A->x, Symbolic, &Common);
            end = clock();

            total += end - begin;
            printf("%ld\n", total);
            klu_free_symbolic(&Symbolic, &Common);
            klu_free_numeric(&Numeric, &Common);
        }

        printf("Time:%lfus\n", total / (double)runtime);
        // double *b;
        // b = klu_malloc(A->nrow, sizeof(double), &Common);
        // for (int i = 0; i < A->nrow; i++)
        //     b[i] = 3;

        // solved = klu_solve(Symbolic, Numeric, A->nrow, 1, b, &Common);
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