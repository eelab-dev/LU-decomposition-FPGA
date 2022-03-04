#include "klu_kernel.h"

void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lip[],
    int Llen[],
    double LU[],
    // int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X[])
{
    double x, *Lx;
    int *Li;
    int k, p, len;

    for (k = 0; k < n; k++)
    {
        // printf("k=%d\n", k);
        x = X[k];
        GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
        /* unit diagonal of L is not stored*/
        for (p = 0; p < len; p++)
        {
            /* X [Li [p]] -= Lx [p] * x [0] ; */
            MULT_SUB(X[Li[p]], Lx[p], x);
        }
    }
}

/* ========================================================================== */
/* === KLU_usolve =========================================================== */
/* ========================================================================== */

/* Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
 * range 1 to 4. */

void KLU_usolve(
    /* inputs, not modified: */
    int n,
    int Uip[],
    int Ulen[],
    double LU[],
    double Udiag[],
    // int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    double X[])
{
    double x, *Ux;
    int *Ui, k, p, len;

    for (k = n - 1; k >= 0; k--)
    {
        GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
        /* x [0] = X [k] / Udiag [k] ; */
        DIV(x, X[k], Udiag[k]);
        X[k] = x;
        for (p = 0; p < len; p++)
        {
            /* X [Ui [p]] -= Ux [p] * x [0] ; */
            MULT_SUB(X[Ui[p]], Ux[p], x);
        }
    }
}

int KLU_solve(
    /* inputs, not modified */
    KLU_symbolic *Symbolic,
    KLU_numeric *Numeric,
    int d, /* leading dimension of B */
    // int nrhs, /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B[], /* size n*nrhs, in column-oriented form, with
                 * leading dimension d. */
    /* --------------- */
    KLU_common *Common)
{
    double x, s;
    int k1, k2, nk, pend, p, lusize_sum = Numeric->lusize_sum;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE);
    }
    // if (Numeric == NULL || Symbolic == NULL || d < Symbolic->n || nrhs < 0 ||
    //     B == NULL)
    // {
    //     Common->status = KLU_INVALID;
    //     return (FALSE);
    // }
    Common->status = KLU_OK;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    ASSERT(KLU_valid(Symbolic->n, Numeric->Offp, Numeric->Offi, Numeric->Offx));

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

    /* ------------------------------------------------------------------ */
    /* get the size of the current chunk */
    /* ------------------------------------------------------------------ */

    /* ------------------------------------------------------------------ */
    /* scale and permute the right hand side, X = P*(R\B) */
    /* ------------------------------------------------------------------ */

    if (Numeric->Rs == NULL)
    {

        /* no scaling */

        for (int k = 0; k < Symbolic->n; k++)
        {
            Numeric->Xwork[k] = B[Numeric->Pnum[k]];
        }
    }
    else
    {

        for (int k = 0; k < Symbolic->n; k++)
        {
            SCALE_DIV_ASSIGN(Numeric->Xwork[k], B[Numeric->Pnum[k]], Numeric->Rs[k]);
        }
    }

    /* ------------------------------------------------------------------ */
    /* solve X = (L*U + Off)\X */
    /* ------------------------------------------------------------------ */

    for (int block = Symbolic->nblocks - 1; block >= 0; block--)
    {

        /* -------------------------------------------------------------- */
        /* the block of size nk is from rows/columns k1 to k2-1 */
        /* -------------------------------------------------------------- */

        k1 = Symbolic->R[block];
        k2 = Symbolic->R[block + 1];
        nk = k2 - k1;
        PRINTF(("solve %d, k1 %d k2-1 %d nk %d\n", block, k1, k2 - 1, nk));

        // printf("nr=%d\n", nr);
        /* solve the block system */
        if (nk == 1)
        {
            s = Numeric->Udiag[k1];
            // printf("X[%d]=%lf,s=%lf\n", k1, X[k1], s);

            DIV(Numeric->Xwork[k1], Numeric->Xwork[k1], s);
        }
        else
        {
            lusize_sum -= Numeric->LUsize[block];
            KLU_lsolve(nk, Numeric->Lip + k1, Numeric->Llen + k1, &Numeric->LUbx[lusize_sum], Numeric->Xwork + k1);
            KLU_usolve(nk, Numeric->Uip + k1, Numeric->Ulen + k1, &Numeric->LUbx[lusize_sum], Numeric->Udiag + k1, Numeric->Xwork + k1);
        }

        /* -------------------------------------------------------------- */
        /* block back-substitution for the off-diagonal-block entries */
        /* -------------------------------------------------------------- */

        if (block > 0)
        {

            for (int k = k1; k < k2; k++)
            {
                // printf("Offp[%d]=%d,Offp[%d]=%d\n", k, Offp[k], k + 1, Offp[k + 1]);
                pend = Numeric->Offp[k + 1];
                x = Numeric->Xwork[k];
                for (p = Numeric->Offp[k]; p < pend; p++)
                {
                    // printf("X[%d]=%lf,Offx[%d]=%lf,x[0]=%lf\n", Numeric->Offi[p], Numeric->Xwork[Numeric->Offi[p]], p, Numeric->Offx[p], x);
                    MULT_SUB(Numeric->Xwork[Numeric->Offi[p]], Numeric->Offx[p], x);
                }
            }
        }
    }

    /* ------------------------------------------------------------------ */
    /* permute the result, Bz  = Q*X */
    /* ------------------------------------------------------------------ */

    for (int k = 0; k < Symbolic->n; k++)
    {
        B[Symbolic->Q[k]] = Numeric->Xwork[k];
    }

    /* ------------------------------------------------------------------ */
    /* go to the next chunk of B */
    /* ------------------------------------------------------------------ */

    return (TRUE);
}