#include "klu_kernel.h"

static void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lip[],
    int Llen[],
    double LU[],
    // int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X[])
{
klu_lsolve_loop:
    for (int k = 0; k < n; k++)
    {
#pragma HLS pipeline
        double x = X[k], *Lx;
        int len, *Li;
        GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
        /* unit diagonal of L is not stored*/
        for (int p = 0; p < len; p++)
        {
#pragma HLS pipeline
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

static void KLU_usolve(
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
klu_usolve_loop:
    for (int k = n - 1; k >= 0; k--)
    {
#pragma HLS pipeline
        double x, *Ux;
        int *Ui, len;
        GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
        /* x [0] = X [k] / Udiag [k] ; */
        DIV(x, X[k], Udiag[k]);
        X[k] = x;
        for (int p = 0; p < len; p++)
        {
#pragma HLS pipeline
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
            double s = Numeric->Udiag[k1];
            DIV(Numeric->Xwork[k1], Numeric->Xwork[k1], s);
        }
        else
        {
            KLU_lsolve(nk, Numeric->Lip + k1, Numeric->Llen + k1, &Numeric->LUbx[Numeric->LUsize[block]], Numeric->Xwork + k1);
            KLU_usolve(nk, Numeric->Uip + k1, Numeric->Ulen + k1, &Numeric->LUbx[Numeric->LUsize[block]], Numeric->Udiag + k1, Numeric->Xwork + k1);
        }

        /* -------------------------------------------------------------- */
        /* block back-substitution for the off-diagonal-block entries */
        /* -------------------------------------------------------------- */

        if (block > 0)
        {
            for (int k = k1; k < k2; k++)
            {
                int pend = Numeric->Offp[k + 1];
                double x = Numeric->Xwork[k];
                for (int p = Numeric->Offp[k]; p < pend; p++)
                {
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