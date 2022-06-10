#include "klu_kernel.h"
#include <ap_int.h>

static void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lip[],
    int Llen[],
    double LU[],
    int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X[][16],
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

            if (r >= 0 && nrhs == 16)
            {
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
                x[0] = X[s][0];
                x[1] = X[s][1];
                x[2] = X[s][2];
                x[3] = X[s][3];
                x[4] = X[s][4];
                x[5] = X[s][5];
                x[6] = X[s][6];
                x[7] = X[s][7];
                x[8] = X[s][8];
                x[9] = X[s][9];
                x[10] = X[s][10];
                x[11] = X[s][11];
                x[12] = X[s][12];
                x[13] = X[s][13];
                x[14] = X[s][14];
                x[15] = X[s][15];
                X[r][0] -= lik * x[0];
                X[r][1] -= lik * x[1];
                X[r][2] -= lik * x[2];
                X[r][3] -= lik * x[3];
                X[r][4] -= lik * x[4];
                X[r][5] -= lik * x[5];
                X[r][6] -= lik * x[6];
                X[r][7] -= lik * x[7];
                X[r][8] -= lik * x[8];
                X[r][9] -= lik * x[9];
                X[r][10] -= lik * x[10];
                X[r][11] -= lik * x[11];
                X[r][12] -= lik * x[12];
                X[r][13] -= lik * x[13];
                X[r][14] -= lik * x[14];
                X[r][15] -= lik * x[15];
            }
            else
            {
                for (int j = 0; j < nrhs; j++)
                {
#pragma HLS LOOP_TRIPCOUNT max = 16
                    x[j] = X[s][j];
                }
                for (int j = 0; j < nrhs; j++)
                {
#pragma HLS LOOP_TRIPCOUNT max = 16
                    X[r][j] -= lik * x[j];
                }
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
    double X[][16],
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

        if (r >= 0 && nrhs == 16)
        {
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
            x[0] = X[r][0] / Udiag[k];
            x[1] = X[r][1] / Udiag[k];
            x[2] = X[r][2] / Udiag[k];
            x[3] = X[r][3] / Udiag[k];
            x[4] = X[r][4] / Udiag[k];
            x[5] = X[r][5] / Udiag[k];
            x[6] = X[r][6] / Udiag[k];
            x[7] = X[r][7] / Udiag[k];
            x[8] = X[r][8] / Udiag[k];
            x[9] = X[r][9] / Udiag[k];
            x[10] = X[r][10] / Udiag[k];
            x[11] = X[r][11] / Udiag[k];
            x[12] = X[r][12] / Udiag[k];
            x[13] = X[r][13] / Udiag[k];
            x[14] = X[r][14] / Udiag[k];
            x[15] = X[r][15] / Udiag[k];
            X[r][0] = x[0];
            X[r][1] = x[1];
            X[r][2] = x[2];
            X[r][3] = x[3];
            X[r][4] = x[4];
            X[r][5] = x[5];
            X[r][6] = x[6];
            X[r][7] = x[7];
            X[r][8] = x[8];
            X[r][9] = x[9];
            X[r][10] = x[10];
            X[r][11] = x[11];
            X[r][12] = x[12];
            X[r][13] = x[13];
            X[r][14] = x[14];
            X[r][15] = x[15];
        }
        else
        {
            for (int j = 0; j < nrhs; j++)
            {
#pragma HLS LOOP_TRIPCOUNT max = 16
                x[j] = X[r][j] / Udiag[k];
                X[r][j] = x[j];
            }
        }

    klu_usolve_loop_2:
        for (int p = 0; p < len; p++)
        {
            int s = Ui[p] * nrhs;
            double uik = Ux[p];
            /* X [Ui [p]] -= Ux [p] * x [0] ; */

            if (s >= 0 && nrhs == 16)
            {
#pragma HLS DEPENDENCE variable = X type = intra false
#pragma HLS DEPENDENCE variable = X type = inter false
                X[s][0] -= uik * x[0];
                X[s][1] -= uik * x[1];
                X[s][2] -= uik * x[2];
                X[s][3] -= uik * x[3];
                X[s][4] -= uik * x[4];
                X[s][5] -= uik * x[5];
                X[s][6] -= uik * x[6];
                X[s][7] -= uik * x[7];
                X[s][8] -= uik * x[8];
                X[s][9] -= uik * x[9];
                X[s][10] -= uik * x[10];
                X[s][11] -= uik * x[11];
                X[s][12] -= uik * x[12];
                X[s][13] -= uik * x[13];
                X[s][14] -= uik * x[14];
                X[s][15] -= uik * x[15];
            }
            else
            {
                for (int j = 0; j < nrhs; j++)
                {
#pragma HLS LOOP_TRIPCOUNT max = 16
                    X[s][j] -= uik * x[j];
                }
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
    double Xwork[MAX_NNZ][16], xusolve[16];
#pragma HLS ARRAY_PARTITION variable = Xwork type = block factor = 16 dim = 2
#pragma HLS ARRAY_PARTITION variable = xusolve

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
                // #pragma HLS DEPENDENCE variable = Xwork type = intra dependent = false
                // #pragma HLS unroll factor = 16
                Xwork[k * nr][j] = B[Pnum[k] + n * j] / Rs[k];
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
                int r = k1 * nr;
                double s = Udiag[k1];

                if (nrhs == 16)
                {
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
                    Xwork[r][0] /= s;
                    Xwork[r][1] /= s;
                    Xwork[r][2] /= s;
                    Xwork[r][3] /= s;
                    Xwork[r][4] /= s;
                    Xwork[r][5] /= s;
                    Xwork[r][6] /= s;
                    Xwork[r][7] /= s;
                    Xwork[r][8] /= s;
                    Xwork[r][9] /= s;
                    Xwork[r][10] /= s;
                    Xwork[r][11] /= s;
                    Xwork[r][12] /= s;
                    Xwork[r][13] /= s;
                    Xwork[r][14] /= s;
                    Xwork[r][15] /= s;
                }
                else
                {
                    for (int j = 0; j < nr; j++)
                    {
#pragma HLS LOOP_TRIPCOUNT max = 16
                        Xwork[r][j] /= s;
                    }
                }
            }
            else
            {
                KLU_lsolve(nk, Lip + k1, Llen + k1, &LUbx[LUsize[block]], nr, Xwork + nr, xusolve);
                KLU_usolve(nk, Uip + k1, Ulen + k1, &LUbx[LUsize[block]], Udiag + k1, nr, Xwork + nr, xusolve);
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
                        int r = Offi[p] * nr, s = k * nr;

                        if (nrhs == 16)
                        {
#pragma HLS DEPENDENCE variable = Xwork type = intra false
#pragma HLS DEPENDENCE variable = Xwork type = inter false
                            Xwork[r][0] -= Offx[p] * Xwork[s][0];
                            Xwork[r][1] -= Offx[p] * Xwork[s][1];
                            Xwork[r][2] -= Offx[p] * Xwork[s][2];
                            Xwork[r][3] -= Offx[p] * Xwork[s][3];
                            Xwork[r][4] -= Offx[p] * Xwork[s][4];
                            Xwork[r][5] -= Offx[p] * Xwork[s][5];
                            Xwork[r][6] -= Offx[p] * Xwork[s][6];
                            Xwork[r][7] -= Offx[p] * Xwork[s][7];
                            Xwork[r][8] -= Offx[p] * Xwork[s][8];
                            Xwork[r][9] -= Offx[p] * Xwork[s][9];
                            Xwork[r][10] -= Offx[p] * Xwork[s][10];
                            Xwork[r][11] -= Offx[p] * Xwork[s][11];
                            Xwork[r][12] -= Offx[p] * Xwork[s][12];
                            Xwork[r][13] -= Offx[p] * Xwork[s][13];
                            Xwork[r][14] -= Offx[p] * Xwork[s][14];
                            Xwork[r][15] -= Offx[p] * Xwork[s][15];
                        }
                        else
                        {
                            for (int j = 0; j < nr; j++)
                            {
#pragma HLS LOOP_TRIPCOUNT max = 16
                                Xwork[r][j] -= Offx[p] * Xwork[s][j];
                            }
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
                // #pragma HLS DEPENDENCE variable = Xwork type = intra dependent = false
                // #pragma HLS DEPENDENCE variable = Xwork type = inter dependent = false
                // #pragma HLS unroll factor = 16
                B[Q[k] + n * j] = Xwork[k * nr][j];
            }
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */
        B += n * 16;
    }

    return (TRUE);
}
