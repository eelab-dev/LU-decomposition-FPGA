#include "klu_kernel.h"

static void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lip[],
    int Llen[],
    double LU[],
    int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X[])
{
klu_lsolve_loop:
    for (int k = 0; k < n; k++)
    {
        double *Lx;
        int len, *Li;
        GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
        /* unit diagonal of L is not stored*/
        for (int p = 0; p < len; p++)
        {
            int r = Li[p] * nrhs, s = k * nrhs;
            double lik = Lx[p];
            /* X [Li [p]] -= Lx [p] * x [0] ; */
            for (int j = 0; j < nrhs; j++)
                MULT_SUB(X[r + j], lik, X[s + j]);
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

static void KLU_usolve(
    /* inputs, not modified: */
    int n,
    int Uip[],
    int Ulen[],
    double LU[],
    double Udiag[],
    int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    double X[])
{
    double x[nrhs];
klu_usolve_loop:
    for (int k = n - 1; k >= 0; k--)
    {
        double *Ux;
        int *Ui, len;
        int r = k * nrhs;
        GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
        /* x [0] = X [k] / Udiag [k] ; */
        for (int j = 0; j < nrhs; j++)
        {
            DIV(x[j], X[r + j], Udiag[k]);
            X[r + j] = x[j];
        }

        for (int p = 0; p < len; p++)
        {
            int s = Ui[p] * nrhs;
            double uik = Ux[p];
            /* X [Ui [p]] -= Ux [p] * x [0] ; */
            for (int j = 0; j < nrhs; j++)
                MULT_SUB(X[s + j], uik, x[j]);
        }
    }
}

int KLU_solve(
    /* inputs, not modified */
    KLU_symbolic *Symbolic,
    KLU_numeric *Numeric,
    int n,    /* leading dimension of B */
    int nrhs, /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B[], /* size n*nrhs, in column-oriented form, with
                 * leading dimension d. */
    /* --------------- */
    KLU_common *Common)
{
    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE);
    }
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
            for (int j = 0; j < nrhs; j++)
                Numeric->Xwork[k * nrhs + j] = B[Numeric->Pnum[k] + n * j];
        }
    }
    else
    {
        // for (int k = 0; k < Symbolic->n; k++)
        // {
        //     for (int j = 0; j < nrhs; j++)
        //     {
        //         std::swap(B[k * nrhs + j], B[Numeric->Pnum[k] + n * j]);
        //     }
        // }

        for (int k = 0; k < Symbolic->n; k++)
        {
            for (int j = 0; j < nrhs; j++)
            {
                SCALE_DIV_ASSIGN(Numeric->Xwork[k * nrhs + j], B[Numeric->Pnum[k] + n * j], Numeric->Rs[k]);
                // printf("Xwork[%d] = %f\n", k * nrhs + j, Numeric->Xwork[k * nrhs + j]);
            }
        }
    }

    /* ------------------------------------------------------------------ */
    /* solve X = (L*U + Off)\X */
    /* ------------------------------------------------------------------ */

klu_solve_loop:
    for (int block = Symbolic->nblocks - 1; block >= 0; block--)
    {

        /* -------------------------------------------------------------- */
        /* the block of size nk is from rows/columns k1 to k2-1 */
        /* -------------------------------------------------------------- */

        int k1 = Symbolic->R[block];
        int k2 = Symbolic->R[block + 1];
        int nk = k2 - k1;
        PRINTF(("solve %d, k1 %d k2-1 %d nk %d\n", block, k1, k2 - 1, nk));

        /* solve the block system */
        if (nk == 1)
        {
            int r = k1 * nrhs;
            double s = Numeric->Udiag[k1];
            for (int j = 0; j < nrhs; j++)
                DIV(Numeric->Xwork[r + j], Numeric->Xwork[r + j], s);
        }
        else
        {
            KLU_lsolve(nk, Numeric->Lip + k1, Numeric->Llen + k1, &Numeric->LUbx[Numeric->LUsize[block]], nrhs, Numeric->Xwork + nrhs * k1);
            KLU_usolve(nk, Numeric->Uip + k1, Numeric->Ulen + k1, &Numeric->LUbx[Numeric->LUsize[block]], Numeric->Udiag + k1, nrhs, Numeric->Xwork + nrhs * k1);
        }

        /* -------------------------------------------------------------- */
        /* block back-substitution for the off-diagonal-block entries */
        /* -------------------------------------------------------------- */

        if (block > 0)
        {
            for (int k = k1; k < k2; k++)
            {
                int pend = Numeric->Offp[k + 1];
                // double x = Numeric->Xwork[k];
                for (int p = Numeric->Offp[k]; p < pend; p++)
                {
                    for (int j = 0; j < nrhs; j++)
                        MULT_SUB(Numeric->Xwork[Numeric->Offi[p] * nrhs + j], Numeric->Offx[p], Numeric->Xwork[k * nrhs + j]);
                }
            }
        }
    }

    /* ------------------------------------------------------------------ */
    /* permute the result, Bz  = Q*X */
    /* ------------------------------------------------------------------ */

    for (int k = 0; k < Symbolic->n; k++)
    {
        for (int j = 0; j < nrhs; j++)
            B[Symbolic->Q[k] + n * j] = Numeric->Xwork[k * nrhs + j];
    }

    /* ------------------------------------------------------------------ */
    /* go to the next chunk of B */
    /* ------------------------------------------------------------------ */

    return (TRUE);
}