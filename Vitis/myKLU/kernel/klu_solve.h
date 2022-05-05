#include "klu_kernel.h"

static void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lip[],
    int Llen[],
    double LU[],
    int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X[],
    double x[])
{
    int r, s, len, *Li;
    double lik, *Lx;
klu_lsolve_loop:
    for (int k = 0; k < n; k++)
    {
        GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
    /* unit diagonal of L is not stored*/
    klu_lsolve_loop_2:
        for (int p = 0; p < len; p++)
        {
            /* X [Li [p]] -= Lx [p] * x [0] ; */
            r = Li[p] * nrhs;
            s = k * nrhs;
            lik = Lx[p];
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS unroll factor = 16
                x[j] = X[s + j];
            }
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS DEPENDENCE variable = X type = intra dependent = false
#pragma HLS DEPENDENCE variable = X type = inter dependent = false
#pragma HLS unroll factor = 16
                X[r + j] -= lik * x[j];
            }
        }
    }
}

/* ========================================================================== */
/* === KLU_usolve =========================================================== */
/* ========================================================================== */

/* Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs. */
static void KLU_usolve(
    /* inputs, not modified: */
    int n,
    int Uip[],
    int Ulen[],
    double LU[],
    double Udiag[],
    int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    double X[],
    double x[])
{
    double *Ux;
    int *Ui, len;
klu_usolve_loop:
    for (int k = n - 1; k >= 0; k--)
    {
        GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
        /* x [0] = X [k] / Udiag [k] ; */
        int r = k * nrhs;
        for (int j = 0; j < nrhs; j++)
        {
#pragma HLS unroll factor = 16
            x[j] = X[r + j] / Udiag[k];
            X[r + j] = x[j];
        }
    klu_usolve_loop_2:
        for (int p = 0; p < len; p++)
        {
            //#pragma HLS pipeline off
            int s = Ui[p] * nrhs;
            double uik = Ux[p];
            /* X [Ui [p]] -= Ux [p] * x [0] ; */
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS DEPENDENCE variable = X type = intra dependent = false
#pragma HLS DEPENDENCE variable = X type = inter dependent = false
#pragma HLS unroll factor = 16
                X[s + j] -= uik * x[j];
            }
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
            {
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = intra dependent = false
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = inter dependent = false
#pragma HLS unroll factor = 16
                Numeric->Xwork[k * nrhs + j] = B[Numeric->Pnum[k] + n * j];
            }
        }
    }
    else
    {
        for (int k = 0; k < Symbolic->n; k++)
        {
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = intra dependent = false
#pragma HLS unroll factor = 16
                Numeric->Xwork[k * nrhs + j] = B[Numeric->Pnum[k] + n * j] / Numeric->Rs[k];
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
            double s = Numeric->Udiag[k1];
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = intra dependent = false
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = inter dependent = false
#pragma HLS unroll factor = 16
                Numeric->Xwork[k1 * nrhs + j] /= s;
            }
        }
        else
        {
            KLU_lsolve(nk, Numeric->Lip + k1, Numeric->Llen + k1, &Numeric->LUbx[Numeric->LUsize[block]], nrhs, Numeric->Xwork + nrhs * k1, Numeric->xusolve);
            KLU_usolve(nk, Numeric->Uip + k1, Numeric->Ulen + k1, &Numeric->LUbx[Numeric->LUsize[block]], Numeric->Udiag + k1, nrhs, Numeric->Xwork + nrhs * k1, Numeric->xusolve);
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
                    int r = Numeric->Offi[p] * nrhs, s = k * nrhs;
                    for (int j = 0; j < nrhs; j++)
                    {
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = intra dependent = false
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = inter dependent = false
#pragma HLS unroll factor = 16
                        Numeric->Xwork[r + j] -= Numeric->Offx[p] * Numeric->Xwork[s + j];
                    }
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
        {
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = intra dependent = false
#pragma HLS DEPENDENCE variable = Numeric->Xwork type = inter dependent = false
#pragma HLS unroll factor = 16
            B[Symbolic->Q[k] + n * j] = Numeric->Xwork[k * nrhs + j];
        }
    }

    /* ------------------------------------------------------------------ */
    /* go to the next chunk of B */
    /* ------------------------------------------------------------------ */

    return (TRUE);
}
