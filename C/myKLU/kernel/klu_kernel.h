/* ========================================================================== */
/* === KLU/Include/klu_internal.h =========================================== */
/* ========================================================================== */

/* For internal use in KLU routines only, not for user programs */

#ifndef _KLU_KERNEL_H
#define _KLU_KERNEL_H

// #include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 65536

#ifndef NDEBUG
#define NDEBUG
#endif
#ifndef NPRINT
#define NPRINT
#endif

/* To enable debugging and assertions, uncomment this line:
 #undef NDEBUG
 */

/* To enable diagnostic printing, uncomment this line:
#undef NPRINT
*/

#undef ASSERT
#ifndef NDEBUG
#define ASSERT(a) assert(a)
#else
#define ASSERT(a)
#endif

#ifndef NPRINT
#define PRINTF(s) printf s
#else
#define PRINTF(s)
#endif

#define PRINT_SCALAR(a)

/* true if an integer (stored in double x) would overflow (or if x is NaN) */
#define INT_OVERFLOW(x) ((!((x) * (1.0 + 1e-8) <= (double)__INT_MAX__)) || SCALAR_IS_NAN(x))

#define TRUE 1
#define FALSE 0
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP(i) : (i))

#define SCALAR_IS_NAN(x) ((x) != (x))
#define SCALAR_IS_ZERO(x) ((x) == 0.)
#define SCALAR_IS_NONZERO(x) ((x) != 0.)
#define SCALAR_IS_LTZERO(x) ((x) < 0.)

#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO(x)) ? -(x) : (x))

#define SPLIT(s) (1)
#define REAL(c) (c)
#define IMAG(c) (0.)
#define ASSIGN(c, s1, s2, p, split) \
    {                               \
        (c) = (s1)[p];              \
    }
#define CLEAR(c)  \
    {             \
        (c) = 0.; \
    }
#define CLEAR_AND_INCREMENT(p) \
    {                          \
        *p++ = 0.;             \
    }
#define IS_NAN(a) SCALAR_IS_NAN(a)
#define IS_ZERO(a) SCALAR_IS_ZERO(a)
#define IS_NONZERO(a) SCALAR_IS_NONZERO(a)
#define SCALE_DIV(c, s) \
    {                   \
        (c) /= (s);     \
    }
#define SCALE_DIV_ASSIGN(a, c, s) \
    {                             \
        a = c / s;                \
    }
#define SCALE(c, s) \
    {               \
        (c) *= (s); \
    }
#define ASSEMBLE(c, a) \
    {                  \
        (c) += (a);    \
    }
#define ASSEMBLE_AND_INCREMENT(c, p) \
    {                                \
        (c) += *p++;                 \
    }
#define DECREMENT(c, a) \
    {                   \
        (c) -= (a);     \
    }
#define MULT(c, a, b)    \
    {                    \
        (c) = (a) * (b); \
    }
#define MULT_CONJ(c, a, b) \
    {                      \
        (c) = (a) * (b);   \
    }
#define MULT_SUB(c, a, b) \
    {                     \
        (c) -= (a) * (b); \
    }
#define MULT_SUB_CONJ(c, a, b) \
    {                          \
        (c) -= (a) * (b);      \
    }
#define DIV(c, a, b)     \
    {                    \
        (c) = (a) / (b); \
    }
#define RECIPROCAL(c)    \
    {                    \
        (c) = 1.0 / (c); \
    }
#define DIV_CONJ(c, a, b) \
    {                     \
        (c) = (a) / (b);  \
    }
#define APPROX_ABS(s, a)     \
    {                        \
        (s) = SCALAR_ABS(a); \
    }
#define ABS(s, a)            \
    {                        \
        (s) = SCALAR_ABS(a); \
    }

#define GET_POINTER(LU, Xip, Xlen, Xi, Xx, k, xlen) \
    {                                               \
        double *xp = LU + Xip[k];                   \
        xlen = Xlen[k];                             \
        Xi = (int *)xp;                             \
        Xx = (double *)(xp + UNITS(int, xlen));     \
    }
#define PRINT_ENTRY(a) PRINT_SCALAR(a)
#define CONJ(a, x) a = x

#define KLU_scale klu_scale
#define KLU_solve klu_solve
#define KLU_factor klu_factor
#define KLU_lsolve klu_lsolve
#define KLU_usolve klu_usolve
#define KLU_kernel klu_kernel
#define KLU_defaults klu_defaults

#define KLU_symbolic klu_symbolic
#define KLU_numeric klu_numeric2
#define KLU_common klu_common

#define BYTES(type, n) (sizeof(type) * (n))
#define CEILING(b, u) (((b) + (u)-1) / (u))
#define UNITS(type, n) (CEILING(BYTES(type, n), sizeof(double)))
#define DUNITS(type, n) (ceil(BYTES(type, (double)n) / sizeof(double)))

#define KLU_OK 0
#define KLU_SINGULAR (1) /* status > 0 is a warning, not an error */
#define KLU_OUT_OF_MEMORY (-2)
#define KLU_INVALID (-3)
#define KLU_TOO_LARGE (-4) /* integer overflow has occured */

// typedef struct klu_common_struct2
// {

//     /* ---------------------------------------------------------------------- */
//     /* parameters */
//     /* ---------------------------------------------------------------------- */

//     double tol;         /* pivot tolerance for diagonal preference */
//     double memgrow;     /* realloc memory growth size for LU factors */
//     double initmem_amd; /* init. memory size with AMD: c*nnz(L) + n */
//     double initmem;     /* init. memory size: c*nnz(A) + n */
//     double maxwork;     /* maxwork for BTF, <= 0 if no limit */

//     int btf;      /* use BTF pre-ordering, or not */
//     int ordering; /* 0: AMD, 1: COLAMD, 2: user P and Q,
//                    * 3: user function */
//     int scale;    /* row scaling: -1: none (and no error check),
//                    * 0: none, 1: sum, 2: max */

//     /* pointer to user ordering function */
//     int (*user_order)(int, int *, int *, int *, struct klu_common_struct2 *);

//     /* pointer to user data, passed unchanged as the last parameter to the
//      * user ordering function (optional, the user function need not use this
//      * information). */
//     void *user_data;

//     int halt_if_singular; /* how to handle a singular matrix:
//                            * FALSE: keep going.  Return a Numeric object with a zero U(k,k).  A
//                            *   divide-by-zero may occur when computing L(:,k).  The Numeric object
//                            *   can be passed to klu_solve (a divide-by-zero will occur).  It can
//                            *   also be safely passed to klu_refactor.
//                            * TRUE: stop quickly.  klu_factor will free the partially-constructed
//                            *   Numeric object.  klu_refactor will not free it, but will leave the
//                            *   numerical values only partially defined.  This is the default. */

//     /* ---------------------------------------------------------------------- */
//     /* statistics */
//     /* ---------------------------------------------------------------------- */

//     int status;   /* KLU_OK if OK, < 0 if error */
//     int nrealloc; /* # of reallocations of L and U */

//     int structural_rank; /* 0 to n-1 if the matrix is structurally rank
//                           * deficient (as determined by maxtrans).  -1 if not computed.  n if the
//                           * matrix has full structural rank.  This is computed by klu_analyze
//                           * if a BTF preordering is requested. */

//     int numerical_rank; /* First k for which a zero U(k,k) was found,
//                          * if the matrix was singular (in the range 0 to n-1).  n if the matrix
//                          * has full rank. This is not a true rank-estimation.  It just reports
//                          * where the first zero pivot was found.  -1 if not computed.
//                          * Computed by klu_factor and klu_refactor. */

//     int singular_col; /* n if the matrix is not singular.  If in the
//                        * range 0 to n-1, this is the column index of the original matrix A that
//                        * corresponds to the column of U that contains a zero diagonal entry.
//                        * -1 if not computed.  Computed by klu_factor and klu_refactor. */

//     int noffdiag; /* # of off-diagonal pivots, -1 if not computed */

//     double flops;   /* actual factorization flop count, from klu_flops */
//     double rcond;   /* crude reciprocal condition est., from klu_rcond */
//     double condest; /* accurate condition est., from klu_condest */
//     double rgrowth; /* reciprocal pivot rgrowth, from klu_rgrowth */
//     double work;    /* actual work done in BTF, in klu_analyze */

//     int memusage; /* current memory usage, in bytes */
//     int mempeak;  /* peak memory usage, in bytes */

// } klu_common;

// typedef struct
// {
//     /* A (P,Q) is in upper block triangular form.  The kth block goes from
//      * row/col index R [k] to R [k+1]-1.  The estimated number of nonzeros
//      * in the L factor of the kth block is Lnz [k].
//      */

//     /* only computed if the AMD ordering is chosen: */
//     double symmetry;  /* symmetry of largest block */
//     double est_flops; /* est. factorization flop count */
//     double lnz, unz;  /* estimated nz in L and U, including diagonals */
//     double *Lnz;      /* size n, but only Lnz [0..nblocks-1] is used */

//     /* computed for all orderings: */
//     int
//         n,        /* input matrix A is n-by-n */
//         nz,       /* # entries in input matrix */
//         *P,       /* size n */
//         *Q,       /* size n */
//         *R,       /* size n+1, but only R [0..nblocks] is used */
//         nzoff,    /* nz in off-diagonal blocks */
//         nblocks,  /* number of blocks */
//         maxblock, /* size of largest block */
//         ordering, /* ordering used (AMD, COLAMD, or GIVEN) */
//         do_btf;   /* whether or not BTF preordering was requested */

//     /* only computed if BTF preordering requested */
//     int structural_rank; /* 0 to n-1 if the matrix is structurally rank
//                           * deficient.  -1 if not computed.  n if the matrix has
//                           * full structural rank */

// } klu_symbolic;

typedef struct
{
    /* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    int n;             /* A is n-by-n */
    int nblocks;       /* number of diagonal blocks */
    int lnz;           /* actual nz in L, including diagonal */
    int unz;           /* actual nz in U, including diagonal */
    int max_lnz_block; /* max actual nz in L in any one block, incl. diag */
    int max_unz_block; /* max actual nz in U in any one block, incl. diag */
    int *Pnum;         /* size n. final pivot permutation */
    int *Pinv;         /* size n. inverse of final pivot permutation */

    /* LU factors of each block */
    int *Lip;      /* size n. pointers into LUbx[block] for L */
    int *Uip;      /* size n. pointers into LUbx[block] for U */
    int *Llen;     /* size n. Llen [k] = # of entries in kth column of L */
    int *Ulen;     /* size n. Ulen [k] = # of entries in kth column of U */
    double *LUbx;  /* L and U indices and entries (excl. diagonal of U) */
    int *LUsize;   /* size of each LUbx [block], in sizeof (Unit) */
    double *Udiag; /* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    double *Rs; /* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    int worksize; /* size (in bytes) of Work */
    // void *Work;      /* workspace */
    double *Xwork; /* alias into Numeric->Work */
    int *Iwork;    /* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    int *Offp;    /* size n+1, column pointers */
    int *Offi;    /* size nzoff, row indices */
    double *Offx; /* size nzoff, numerical values */
    int nzoff;

    int lusize_sum;

} klu_numeric2;

int KLU_kernel /* final size of LU on output */
    (
        /* input, not modified */
        int n,       /* A is n-by-n */
        int Ap[],    /* size n+1, column pointers for A */
        int Ai[],    /* size nz = Ap [n], row indices for A */
        double Ax[], /* size nz, values of A */
        int Q[],     /* size n, optional input permutation */
        int lusize,  /* initial size of LU */

        /* output, not defined on input */
        int Pinv[],     /* size n */
        int P[],        /* size n */
        double *p_LU,   /* size lusize on input, size Uxp[n] on output*/
        double Udiag[], /* size n, diagonal of U */
        int Llen[],     /* size n, column length of L */
        int Ulen[],     /* size n, column length of U */
        int Lip[],      /* size n+1 */
        int Uip[],      /* size n+1 */
        int *lnz,       /* size of L */
        int *unz,       /* size of U */

        /* workspace, not defined on input */
        double X[], /* size n, zero on output */

        /* workspace, not defined on input or output */
        int Stack[],   /* size n */
        int Flag[],    /* size n */
        int adj_pos[], /* size n */

        /* workspace for pruning only */
        int Lpend[], /* size n workspace */

        /* inputs, not modified on output */
        int k1,      /* the block of A is from k1 to k2-1 */
        int PSinv[], /* inverse of P from symbolic factorization */
        double Rs[], /* scale factors for A */

        /* inputs, modified on output */
        int Offp[], /* off-diagonal matrix (modified by this routine) */
        int Offi[],
        double Offx[],
        KLU_common *Common /* the control input/output structure */
    );

void KLU_lsolve(
    /* inputs, not modified: */
    int n,
    int Lp[],
    int Li[],
    double LU[],
    /* right-hand-side on input, solution to Lx=b on output */
    double X[]);

void KLU_usolve(
    /* inputs, not modified: */
    int n,
    int Up[],
    int Ui[],
    double LU[],
    double Udiag[],
    /* right-hand-side on input, solution to Ux=b on output */
    double X[]);

int KLU_defaults(
    KLU_common *Common)
{
    if (Common == NULL)
    {
        return (FALSE);
    }

    /* parameters */
    Common->tol = 0.001;             /* pivot tolerance for diagonal */
    Common->memgrow = 1.2;           /* realloc size ratio increase for LU factors */
    Common->initmem_amd = 1.2;       /* init. mem with AMD:  c*nnz(L) + n */
    Common->initmem = 10;            /* init. mem otherwise: c*nnz(A) + n */
    Common->btf = TRUE;              /* use BTF pre-ordering, or not */
    Common->maxwork = 0;             /* no limit to work done by btf_order */
    Common->ordering = 0;            /* 0: AMD, 1: COLAMD, 2: user-provided P and Q,
                                      * 3: user-provided function */
    Common->scale = 2;               /* scale: -1: none, and do not check for errors
                                      * in the input matrix in KLU_refactor.
                                      * 0: none, but check for errors,
                                      * 1: sum, 2: max */
    Common->halt_if_singular = TRUE; /* quick halt if matrix is singular */

    /* user ordering function and optional argument */
    Common->user_order = NULL;
    Common->user_data = NULL;

    /* statistics */
    Common->status = KLU_OK;
    Common->nrealloc = 0;
    Common->structural_rank = EMPTY;
    Common->numerical_rank = EMPTY;
    Common->noffdiag = EMPTY;
    Common->flops = EMPTY;
    Common->rcond = EMPTY;
    Common->condest = EMPTY;
    Common->rgrowth = EMPTY;
    Common->work = 0; /* work done by btf_order */

    Common->memusage = 0;
    Common->mempeak = 0;

    return (TRUE);
}

#endif