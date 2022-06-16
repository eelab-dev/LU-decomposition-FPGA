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
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
#pragma HLS unroll factor = 16
                x[j] = X[s + j];
            }
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
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
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
            // #pragma HLS unroll factor = 16
            x[j] = X[r + j] / Udiag[k];
            X[r + j] = x[j];
        }
    klu_usolve_loop_2:
        for (int p = 0; p < len; p++)
        {
            int s = Ui[p] * nrhs;
            double uik = Ux[p];
            /* X [Ui [p]] -= Ux [p] * x [0] ; */
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
                // #pragma HLS unroll factor = 16
                X[s + j] -= uik * x[j];
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
    double Xwork[MAX_NNZ], xusolve[MAX_SIZE];

#pragma HLS ARRAY_PARTITION variable=Xwork type=cyclic factor=16
#pragma HLS ARRAY_PARTITION variable=xusolve type=cyclic factor=16


    for (int chunk = 0; chunk < nrhs; chunk += 16)
    {
        int nr = MIN(nrhs - chunk, 16);

        /* ------------------------------------------------------------------ */
        /* scale and permute the right hand side, X = P*(R\B) */
        /* ------------------------------------------------------------------ */

        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < nr; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
#pragma HLS DEPENDENCE variable = B type = intra false
#pragma HLS DEPENDENCE variable = B type = inter false
                // #pragma HLS unroll factor = 16
                Xwork[k * nr + j] = B[Pnum[k] + n * j] / Rs[k];
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
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
                    // #pragma HLS unroll factor = 16
                    Xwork[k1 * nr + j] /= s;
                }
            }
            else
            {
                KLU_lsolve(nk, Lip + k1, Llen + k1, &LUbx[LUsize[block]], nr, Xwork + nr * k1, xusolve);
                KLU_usolve(nk, Uip + k1, Ulen + k1, &LUbx[LUsize[block]], Udiag + k1, nr, Xwork + nr * k1, xusolve);
            }

            /* -------------------------------------------------------------- */
            /* block back-substitution for the off-diagonal-block entries */
            /* -------------------------------------------------------------- */

            if (block > 0)
            {
                for (int k = k1; k < k2; k++)
                {
                    int pend = Offp[k + 1];
                    // double x = Xwork[k];
                    for (int p = Offp[k]; p < pend; p++)
                    {
                        int r = Offi[p] * nr, s = k * nr;
                        for (int j = 0; j < nr; j++)
                        {
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
                             #pragma HLS unroll factor = 16
                            Xwork[r + j] -= Offx[p] * Xwork[s + j];
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
#pragma HLS LOOP_TRIPCOUNT max = 16
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
#pragma HLS DEPENDENCE variable = B type = intra false
#pragma HLS DEPENDENCE variable = B type = inter false
                // #pragma HLS unroll factor = 16
                B[Q[k] + n * j] = Xwork[k * nr + j];
            }
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */
        B += n * 16;
    }

    return (TRUE);
}