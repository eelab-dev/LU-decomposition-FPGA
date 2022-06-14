#include "klu_kernel.h"

static void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lip[],
    int Llen[],
    double LU[],
    int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X[][MAX_RHS])
{
    int len, *Li;
    double *Lx;

klu_lsolve_loop:
    for (int k = 0; k < n; k++)
    {
        GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
    /* unit diagonal of L is not stored*/
    klu_lsolve_loop_2:
        for (int p = 0; p < len; p++)
        {
            /* X [Li [p]] -= Lx [p] * x [0] ; */
            int r = Li[p];
            double lik = Lx[p];
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 64
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
#pragma HLS UNROLL factor = 64
                X[r][j] -= lik * X[k][j];
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
    double X[][MAX_RHS],
    int k1)
{
    double *Ux;
    int *Ui, len;

klu_usolve_loop:
    for (int k = n - 1; k >= 0; k--)
    {
        GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
        /* x [0] = X [k] / Udiag [k] ; */

        double r = Udiag[k];
        for (int j = 0; j < nrhs; j++)
        {
#pragma HLS LOOP_TRIPCOUNT max = 64
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
#pragma HLS UNROLL factor = 64
            X[k][j] /= X[MAX_SIZE - 1 - k1][j] * r;
        }

    klu_usolve_loop_2:
        for (int p = 0; p < len; p++)
        {
            int s = Ui[p];
            double uik = Ux[p];
            /* X [Ui [p]] -= Ux [p] * x [0] ; */
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 64
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
#pragma HLS UNROLL factor = 64
                X[s][j] -= uik * X[k][j];
            }
        }
    }
}

int KLU_solve(
    int *R,
    int *Q,
    int *Pnum,
    double *Rs,
    int *Lip,
    int *Llen,
    double *LUbx,
    int *LUsize,
    int *Uip,
    int *Ulen,
    double *Udiag,
    int *Offp,
    int *Offi,
    double *Offx,
    int n,    /* leading dimension of B */
    int nrhs, /* number of right-hand-sides */
    int nblocks,
    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B[] /* size n*nrhs, in column-oriented form, with
                * leading dimension d. */
    /* --------------- */
)
{
    double Xwork[MAX_SIZE][MAX_RHS];
#pragma HLS ARRAY_PARTITION variable = Xwork type = block factor = 64 dim = 2

    for (int i = 0; i < MAX_RHS; i++)
    {
#pragma HLS UNROLL
        Xwork[MAX_SIZE - 1][i] = 1;
    }

    for (int chunk = 0; chunk < nrhs; chunk += MAX_RHS)
    {
        int nr = MIN(nrhs - chunk, MAX_RHS);

        /* ------------------------------------------------------------------ */
        /* scale and permute the right hand side, X = P*(R\B) */
        /* ------------------------------------------------------------------ */

        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < nr; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 64
                Xwork[k][j] = B[Pnum[k] + n * j] / Rs[k];
            }
        }

        /* ------------------------------------------------------------------ */
        /* solve X = (L*U + Off)\X */
        /* ------------------------------------------------------------------ */

    klu_solve_loop:
        for (int block = nblocks - 1; block >= 0; block--)
        {
            /* -------------------------------------------------------------- */
            /* the block of size nk is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            int k1 = R[block];
            int k2 = R[block + 1];
            int nk = k2 - k1;
            PRINTF(("solve %d, k1 %d k2-1 %d nk %d\n", block, k1, k2 - 1, nk));

            /* solve the block system */
            if (nk == 1)
            {
                double s = Udiag[k1];
                for (int j = 0; j < nr; j++)
                {
#pragma HLS LOOP_TRIPCOUNT max = 64
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
#pragma HLS UNROLL factor = 64
                    Xwork[k1][j] /= Xwork[MAX_SIZE - 1][j] * s;
                }
            }
            else
            {
                KLU_lsolve(nk, Lip + k1, Llen + k1, &LUbx[LUsize[block]], nr, Xwork + k1);
                KLU_usolve(nk, Uip + k1, Ulen + k1, &LUbx[LUsize[block]], Udiag + k1, nr, Xwork + k1, k1);
            }

            /* -------------------------------------------------------------- */
            /* block back-substitution for the off-diagonal-block entries */
            /* -------------------------------------------------------------- */

            if (block > 0)
            {
                for (int k = k1; k < k2; k++)
                {
                    int pend = Offp[k + 1];
                    for (int p = Offp[k]; p < pend; p++)
                    {
                        int r = Offi[p];
                        for (int j = 0; j < nr; j++)
                        {
#pragma HLS LOOP_TRIPCOUNT max = 64
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
#pragma HLS UNROLL factor = 64
                            Xwork[r][j] -= Offx[p] * Xwork[k][j];
                        }
                    }
                }
            }
        }

        /* ------------------------------------------------------------------ */
        /* permute the result, Bz  = Q*X */
        /* ------------------------------------------------------------------ */

        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < nr; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 64
                B[Q[k] + n * j] = Xwork[k][j];
            }
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */
        B += n * MAX_RHS;
    }

    return TRUE;
}
