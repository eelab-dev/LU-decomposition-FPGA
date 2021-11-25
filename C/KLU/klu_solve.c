/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */

#include "cholmod.h"
#include "klu.h"
#include <stdio.h>

int main(void)
{
    cholmod_sparse *A;
    cholmod_common ch;
    cholmod_start(&ch);
    A = cholmod_read_sparse(stdin, &ch);
    if (A)
    {
        if (A->nrow != A->ncol || A->stype != 0 || (!(A->xtype == CHOLMOD_REAL || A->xtype == CHOLMOD_COMPLEX)))
        {
            printf("invalid matrix\n");
        }
        else
        {
            if (A->xtype == CHOLMOD_REAL)
            {
                klu_symbolic *Symbolic;
                klu_numeric *Numeric;
                klu_common Common;

                int n = A->nrow;
                double b[n];
                for (int i = 0; i < n; i++)
                    b[i] = 1 + ((double)i + 1) / ((double)n);

                klu_defaults(&Common);
                Symbolic = klu_analyze(n, A->p, A->i, &Common);
                Numeric = klu_factor(A->p, A->i, A->x, Symbolic, &Common);
                klu_solve(Symbolic, Numeric, n, 1, b, &Common);
                klu_free_symbolic(&Symbolic, &Common);
                klu_free_numeric(&Numeric, &Common);
                for (int i = 0; i < n; i++)
                    printf("x [%d] = %g\n", i, b[i]);
            }
            else
                printf("Complex");
        }
        cholmod_free_sparse(&A, &ch);
    }
    cholmod_finish(&ch);

    return (0);
}
