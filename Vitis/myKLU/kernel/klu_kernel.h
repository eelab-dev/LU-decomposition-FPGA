/* ========================================================================== */
/* === KLU/Include/klu_internal.h =========================================== */
/* ========================================================================== */

/* For internal use in KLU routines only, not for user programs */

#ifndef _KLU_KERNEL_H
#define _KLU_KERNEL_H

#include <limits.h>
#include <stdlib.h>
#include <math.h>

#ifndef NPRINT
#define NPRINT
#endif

#define MAX_SIZE 2048
#define MAX_RHS 32
#define MAX_NNZ 16384

#ifndef NPRINT
#define PRINTF(s) printf s
#else
#define PRINTF(s)
#endif

#define KLU_solve klu_solve
#define KLU_lsolve klu_lsolve
#define KLU_usolve klu_usolve

#define TRUE 1
#define FALSE 0

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define BYTES(type, n) (sizeof(type) * (n))
#define CEILING(b, u) (((b) + (u)-1) / (u))
#define UNITS(type, n) (CEILING(BYTES(type, n), sizeof(double)))

#define GET_POINTER(LU, Xip, Xlen, Xi, Xx, k, xlen) \
    {                                               \
        double *xp = LU + Xip[k];                   \
        xlen = Xlen[k];                             \
        Xi = (int *)xp;                             \
        Xx = (double *)(xp + UNITS(int, xlen));     \
    }

#endif
