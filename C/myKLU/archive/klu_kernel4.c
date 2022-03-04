/* ========================================================================== */
/* === KLU_kernel =========================================================== */
/* ========================================================================== */

/* Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with a non-recursive DFS and with Eisenstat &
 * Liu's symmetric pruning.  No user-callable routines are in this file.
 */
#include "klu.h"
#include "klu_internal4.h"
#include "mmio.h"

/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* Does a depth-first-search, starting at node j. */

static int dfs(
    /* input, not modified on output: */
    int j,      /* node at which to start the DFS */
    int k,      /* mark value, for the Flag array */
    int Pinv[], /* Pinv [i] = k if row i is kth pivot row, or EMPTY if
                 * row i is not yet pivotal.  */
    int Llen[], /* size n, Llen [k] = # nonzeros in column k of L */
    int Lip[],  /* size n, Lip [k] is position in LU of column k of L */

    /* workspace, not defined on input or output */
    int Stack[], /* size n */

    /* input/output: */
    int Flag[],  /* Flag [i] == k means i is marked */
    int Lpend[], /* for symmetric pruning */
    int top,     /* top of stack on input*/
    double LU[],
    int *Lik, /* Li row index array of the kth column */
    int *plength,

    /* other, not defined on input or output */
    int Ap_pos[] /* keeps track of position in adj list during DFS */
)
{
    int pos, jnew, head, l_length;
    int *Li;

    l_length = *plength;

    head = 0;
    Stack[0] = j;
    ASSERT(Flag[j] != k);

    while (head >= 0)
    {
        j = Stack[head];
        jnew = Pinv[j];
        ASSERT(jnew >= 0 && jnew < k); /* j is pivotal */

        if (Flag[j] != k) /* a node is not yet visited */
        {
            /* first time that j has been visited */
            Flag[j] = k;
            PRINTF(("[ start dfs at %d : new %d\n", j, jnew));
            /* set Ap_pos [head] to one past the last entry in col j to scan */
            Ap_pos[head] =
                (Lpend[jnew] == EMPTY) ? Llen[jnew] : Lpend[jnew];
        }

        /* add the adjacent nodes to the recursive stack by iterating through
         * until finding another non-visited pivotal node */
        Li = (int *)(LU + Lip[jnew]);
        for (pos = --Ap_pos[head]; pos >= 0; --pos)
        {
            int i = Li[pos];
            if (Flag[i] != k)
            {
                /* node i is not yet visited */
                if (Pinv[i] >= 0)
                {
                    /* keep track of where we left off in the scan of the
                     * adjacency list of node j so we can restart j where we
                     * left off. */
                    Ap_pos[head] = pos;

                    /* node i is pivotal; push it onto the recursive stack
                     * and immediately break so we can recurse on node i. */
                    Stack[++head] = i;
                    break;
                }
                else
                {
                    /* node i is not pivotal (no outgoing edges). */
                    /* Flag as visited and store directly into L,
                     * and continue with current node j. */
                    Flag[i] = k;
                    Lik[l_length] = i;
                    l_length++;
                }
            }
        }

        if (pos == -1)
        {
            /* if all adjacent nodes of j are already visited, pop j from
             * recursive stack and push j onto output stack */
            head--;
            Stack[--top] = j;
            PRINTF(("  end   dfs at %d ] head : %d\n", j, head));
        }
    }

    *plength = l_length;
    return (top);
}

/* ========================================================================== */
/* === lsolve_symbolic ====================================================== */
/* ========================================================================== */

/* Finds the pattern of x, for the solution of Lx=b */

static int lsolve_symbolic(
    /* input, not modified on output: */
    int n, /* L is n-by-n, where n >= 0 */
    int k, /* also used as the mark value, for the Flag array */
    int Ap[],
    int Ai[],
    int Q[],
    int Pinv[], /* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
                 * is not yet pivotal.  */

    /* workspace, not defined on input or output */
    int Stack[], /* size n */

    /* workspace, defined on input and output */
    int Flag[], /* size n.  Initially, all of Flag [0..n-1] < k.  After
                 * lsolve_symbolic is done, Flag [i] == k if i is in
                 * the pattern of the output, and Flag [0..n-1] <= k. */

    /* other */
    int Lpend[],  /* for symmetric pruning */
    int Ap_pos[], /* workspace used in dfs */

    double LU[], /* LU factors (pattern and values) */
    int lup,     /* pointer to free space in LU */
    int Llen[],  /* size n, Llen [k] = # nonzeros in column k of L */
    int Lip[],   /* size n, Lip [k] is position in LU of column k of L */

    /* ---- the following are only used in the BTF case --- */

    int k1,     /* the block of A is from k1 to k2-1 */
    int PSinv[] /* inverse of P from symbolic factorization */
)
{
    int *Lik;
    int i, p, pend, oldcol, kglobal, top, l_length;

    top = n;
    l_length = 0;
    Lik = (int *)(LU + lup);

    /* ---------------------------------------------------------------------- */
    /* BTF factorization of A (k1:k2-1, k1:k2-1) */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1;    /* column k of the block is col kglobal of A */
    oldcol = Q[kglobal]; /* Q must be present for BTF case */
    pend = Ap[oldcol + 1];
    for (p = Ap[oldcol]; p < pend; p++)
    {
        i = PSinv[Ai[p]] - k1;
        if (i < 0)
            continue; /* skip entry outside the block */

        /* (i,k) is an entry in the block.  start a DFS at node i */
        PRINTF(("\n ===== DFS at node %d in b, inew: %d\n", i, Pinv[i]));
        if (Flag[i] != k)
        {
            if (Pinv[i] >= 0)
            {
                top = dfs(i, k, Pinv, Llen, Lip, Stack, Flag,
                          Lpend, top, LU, Lik, &l_length, Ap_pos);
            }
            else
            {
                /* i is not pivotal, and not flagged. Flag and put in L */
                Flag[i] = k;
                Lik[l_length] = i;
                l_length++;
            }
        }
    }

    /* If Llen [k] is zero, the matrix is structurally singular */
    Llen[k] = l_length;
    return (top);
}

/* ========================================================================== */
/* === construct_column ===================================================== */
/* ========================================================================== */

/* Construct the kth column of A, and the off-diagonal part, if requested.
 * Scatter the numerical values into the workspace X, and construct the
 * corresponding column of the off-diagonal matrix. */

static void construct_column(
    /* inputs, not modified on output */
    int k, /* the column of A (or the column of the block) to get */
    int Ap[],
    int Ai[],
    double Ax[],
    int Q[], /* column pre-ordering */

    /* zero on input, modified on output */
    double X[],

    /* ---- the following are only used in the BTF case --- */

    /* inputs, not modified on output */
    int k1,      /* the block of A is from k1 to k2-1 */
    int PSinv[], /* inverse of P from symbolic factorization */
    double Rs[], /* scale factors for A */
    int scale,   /* 0: no scaling, nonzero: scale the rows with Rs */

    /* inputs, modified on output */
    int Offp[], /* off-diagonal matrix (modified by this routine) */
    int Offi[],
    double Offx[])
{
    double aik;
    int i, p, pend, oldcol, kglobal, poff, oldrow;

    /* ---------------------------------------------------------------------- */
    /* Scale and scatter the column into X. */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1;     /* column k of the block is col kglobal of A */
    poff = Offp[kglobal]; /* start of off-diagonal column */
    oldcol = Q[kglobal];
    pend = Ap[oldcol + 1];

    if (scale <= 0)
    {
        /* no scaling */
        for (p = Ap[oldcol]; p < pend; p++)
        {
            oldrow = Ai[p];
            i = PSinv[oldrow] - k1;
            aik = Ax[p];
            if (i < 0)
            {
                /* this is an entry in the off-diagonal part */
                Offi[poff] = oldrow;
                Offx[poff] = aik;
                poff++;
            }
            else
            {
                /* (i,k) is an entry in the block.  scatter into X */
                X[i] = aik;
            }
        }
    }
    else
    {
        /* row scaling */
        for (p = Ap[oldcol]; p < pend; p++)
        {
            oldrow = Ai[p];
            i = PSinv[oldrow] - k1;
            aik = Ax[p];
            SCALE_DIV(aik, Rs[oldrow]);
            if (i < 0)
            {
                /* this is an entry in the off-diagonal part */
                Offi[poff] = oldrow;
                Offx[poff] = aik;
                poff++;
            }
            else
            {
                /* (i,k) is an entry in the block.  scatter into X */
                X[i] = aik;
            }
        }
    }

    Offp[kglobal + 1] = poff; /* start of the next col of off-diag part */
}

/* ========================================================================== */
/* === lsolve_numeric ======================================================= */
/* ========================================================================== */

/* Computes the numerical values of x, for the solution of Lx=b.  Note that x
 * may include explicit zeros if numerical cancelation occurs.  L is assumed
 * to be unit-diagonal, with possibly unsorted columns (but the first entry in
 * the column must always be the diagonal entry). */

static void lsolve_numeric(
    /* input, not modified on output: */
    int Pinv[],  /* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
                  * is not yet pivotal.  */
    double *LU,  /* LU factors (pattern and values) */
    int Stack[], /* stack for dfs */
    int Lip[],   /* size n, Lip [k] is position in LU of column k of L */
    int top,     /* top of stack on input */
    int n,       /* A is n-by-n */
    int Llen[],  /* size n, Llen [k] = # nonzeros in column k of L */

    /* output, must be zero on input: */
    double X[] /* size n, initially zero.  On output,
                * X [Ui [up1..up-1]] and X [Li [lp1..lp-1]]
                * contains the solution. */

)
{
    double xj;
    double *Lx;
    int *Li;
    int p, s, j, jnew, len;

    /* solve Lx=b */
    for (s = top; s < n; s++)
    {
        /* forward solve with column j of L */
        j = Stack[s];
        jnew = Pinv[j];
        ASSERT(jnew >= 0);
        xj = X[j];
        GET_POINTER(LU, Lip, Llen, Li, Lx, jnew, len);
        ASSERT(Lip[jnew] <= Lip[jnew + 1]);
        for (p = 0; p < len; p++)
        {
            /*X [Li [p]] -= Lx [p] * xj ; */
            MULT_SUB(X[Li[p]], Lx[p], xj);
        }
    }
}

/* ========================================================================== */
/* === lpivot =============================================================== */
/* ========================================================================== */

/* Find a pivot via partial pivoting, and scale the column of L. */

static int lpivot(
    int diagrow,
    int *p_pivrow,
    double *p_pivot,
    double *p_abs_pivot,
    double tol,
    double X[],
    double *LU, /* LU factors (pattern and values) */
    int Lip[],
    int Llen[],
    int k,
    int n,

    int Pinv[], /* Pinv [i] = k if row i is kth pivot row, or EMPTY if
                 * row i is not yet pivotal.  */

    int *p_firstrow,
    KLU_common *Common)
{
    double x, pivot, *Lx;
    double abs_pivot, xabs;
    int p, i, ppivrow, pdiag, pivrow, *Li, last_row_index, firstrow, len;

    pivrow = EMPTY;
    if (Llen[k] == 0)
    {
        /* matrix is structurally singular */
        if (Common->halt_if_singular)
        {
            return (FALSE);
        }
        for (firstrow = *p_firstrow; firstrow < n; firstrow++)
        {
            PRINTF(("check %d\n", firstrow));
            if (Pinv[firstrow] < 0)
            {
                /* found the lowest-numbered non-pivotal row.  Pick it. */
                pivrow = firstrow;
                PRINTF(("Got pivotal row: %d\n", pivrow));
                break;
            }
        }
        ASSERT(pivrow >= 0 && pivrow < n);
        CLEAR(pivot);
        *p_pivrow = pivrow;
        *p_pivot = pivot;
        *p_abs_pivot = 0;
        *p_firstrow = firstrow;
        return (FALSE);
    }

    pdiag = EMPTY;
    ppivrow = EMPTY;
    abs_pivot = EMPTY;
    i = Llen[k] - 1;
    GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
    last_row_index = Li[i];

    /* decrement the length by 1 */
    Llen[k] = i;
    GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);

    /* look in Li [0 ..Llen [k] - 1 ] for a pivot row */
    for (p = 0; p < len; p++)
    {
        /* gather the entry from X and store in L */
        i = Li[p];
        x = X[i];
        CLEAR(X[i]);

        Lx[p] = x;
        /* xabs = ABS (x) ; */
        ABS(xabs, x);

        /* find the diagonal */
        if (i == diagrow)
        {
            pdiag = p;
        }

        /* find the partial-pivoting choice */
        if (xabs > abs_pivot)
        {
            abs_pivot = xabs;
            ppivrow = p;
        }
    }

    /* xabs = ABS (X [last_row_index]) ;*/
    ABS(xabs, X[last_row_index]);
    if (xabs > abs_pivot)
    {
        abs_pivot = xabs;
        ppivrow = EMPTY;
    }

    /* compare the diagonal with the largest entry */
    if (last_row_index == diagrow)
    {
        if (xabs >= tol * abs_pivot)
        {
            abs_pivot = xabs;
            ppivrow = EMPTY;
        }
    }
    else if (pdiag != EMPTY)
    {
        /* xabs = ABS (Lx [pdiag]) ;*/
        ABS(xabs, Lx[pdiag]);
        if (xabs >= tol * abs_pivot)
        {
            /* the diagonal is large enough */
            abs_pivot = xabs;
            ppivrow = pdiag;
        }
    }

    if (ppivrow != EMPTY)
    {
        pivrow = Li[ppivrow];
        pivot = Lx[ppivrow];
        /* overwrite the ppivrow values with last index values */
        Li[ppivrow] = last_row_index;
        Lx[ppivrow] = X[last_row_index];
    }
    else
    {
        pivrow = last_row_index;
        pivot = X[last_row_index];
    }
    CLEAR(X[last_row_index]);

    *p_pivrow = pivrow;
    *p_pivot = pivot;
    *p_abs_pivot = abs_pivot;
    ASSERT(pivrow >= 0 && pivrow < n);

    if (IS_ZERO(pivot) && Common->halt_if_singular)
    {
        /* numerically singular case */
        return (FALSE);
    }

    /* divide L by the pivot value */
    for (p = 0; p < Llen[k]; p++)
    {
        /* Lx [p] /= pivot ; */
        DIV(Lx[p], Lx[p], pivot);
    }

    return (TRUE);
}

/* ========================================================================== */
/* === prune ================================================================ */
/* ========================================================================== */

/* Prune the columns of L to reduce work in subsequent depth-first searches */
static void prune(
    /* input/output: */
    int Lpend[], /* Lpend [j] marks symmetric pruning point for L(:,j) */

    /* input: */
    int Pinv[], /* Pinv [i] = k if row i is kth pivot row, or EMPTY if
                 * row i is not yet pivotal.  */
    int k,      /* prune using column k of U */
    int pivrow, /* current pivot row */

    /* input/output: */
    double *LU, /* LU factors (pattern and values) */

    /* input */
    int Uip[],  /* size n, column pointers for U */
    int Lip[],  /* size n, column pointers for L */
    int Ulen[], /* size n, column length of U */
    int Llen[]  /* size n, column length of L */
)
{
    double x;
    double *Lx, *Ux;
    int *Li, *Ui;
    int p, i, j, p2, phead, ptail, llen, ulen;

    /* check to see if any column of L can be pruned */
    /* Ux is set but not used.  This OK. */
    GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, ulen);
    for (p = 0; p < ulen; p++)
    {
        j = Ui[p];
        ASSERT(j < k);
        PRINTF(("%d is pruned: %d. Lpend[j] %d Lip[j+1] %d\n",
                j, Lpend[j] != EMPTY, Lpend[j], Lip[j + 1]));
        if (Lpend[j] == EMPTY)
        {
            /* scan column j of L for the pivot row */
            GET_POINTER(LU, Lip, Llen, Li, Lx, j, llen);
            for (p2 = 0; p2 < llen; p2++)
            {
                if (pivrow == Li[p2])
                {
                    /* found it!  This column can be pruned */

                    /* partition column j of L.  The unit diagonal of L
                     * is not stored in the column of L. */
                    phead = 0;
                    ptail = Llen[j];
                    while (phead < ptail)
                    {
                        i = Li[phead];
                        if (Pinv[i] >= 0)
                        {
                            /* leave at the head */
                            phead++;
                        }
                        else
                        {
                            /* swap with the tail */
                            ptail--;
                            Li[phead] = Li[ptail];
                            Li[ptail] = i;
                            x = Lx[phead];
                            Lx[phead] = Lx[ptail];
                            Lx[ptail] = x;
                        }
                    }

                    /* set Lpend to one past the last entry in the
                     * first part of the column of L.  Entries in
                     * Li [0 ... Lpend [j]-1] are the only part of
                     * column j of L that needs to be scanned in the DFS.
                     * Lpend [j] was EMPTY; setting it >= 0 also flags
                     * column j as pruned. */
                    Lpend[j] = ptail;
                }
            }
        }
    }
}

/* ========================================================================== */
/* === KLU_kernel =========================================================== */
/* ========================================================================== */

int KLU_kernel /* final size of LU on output */
    (
        /* input, not modified */
        int n,       /* A is n-by-n */
        int Ap[],    /* size n+1, column pointers for A */
        int Ai[],    /* size nz = Ap [n], row indices for A */
        double Ax[], /* size nz, values of A */
        int Q[],     /* size n, optional input permutation */
        int lusize,  /* initial size of LU on input */

        /* output, not defined on input */
        int Pinv[],     /* size n, inverse row permutation, where Pinv [i] = k if
                         * row i is the kth pivot row */
        int P[],        /* size n, row permutation, where P [k] = i if row i is the
                         * kth pivot row. */
        double *p_LU,   /* LU array, size lusize on input */
        double Udiag[], /* size n, diagonal of U */
        int Llen[],     /* size n, column length of L */
        int Ulen[],     /* size n, column length of U */
        int Lip[],      /* size n, column pointers for L */
        int Uip[],      /* size n, column pointers for U */
        int *lnz,       /* size of L*/
        int *unz,       /* size of U*/
        /* workspace, not defined on input */
        double X[], /* size n, undefined on input, zero on output */

        /* workspace, not defined on input or output */
        int Stack[],  /* size n */
        int Flag[],   /* size n */
        int Ap_pos[], /* size n */

        /* other workspace: */
        int Lpend[], /* size n workspace, for pruning only */

        /* inputs, not modified on output */
        int k1,      /* the block of A is from k1 to k2-1 */
        int PSinv[], /* inverse of P from symbolic factorization */
        double Rs[], /* scale factors for A */

        /* inputs, modified on output */
        int Offp[], /* off-diagonal matrix (modified by this routine) */
        int Offi[],
        double Offx[],
        /* --------------- */
        KLU_common *Common)
{
    double pivot;
    double abs_pivot, xsize, nunits, tol, memgrow;
    double *Ux;
    int *Li, *Ui;
    double *LU; /* LU factors (pattern and values) */
    int k, p, i, j, pivrow = 0, kbar, diagrow, firstrow, lup, top, scale, len;
    int newlusize;

    ASSERT(Common != NULL);
    scale = Common->scale;
    tol = Common->tol;
    memgrow = Common->memgrow;
    *lnz = 0;
    *unz = 0;
    CLEAR(pivot);

    /* ---------------------------------------------------------------------- */
    /* get initial Li, Lx, Ui, and Ux */
    /* ---------------------------------------------------------------------- */

    PRINTF(("input: lusize %d \n", lusize));
    ASSERT(lusize > 0);
    LU = p_LU;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    firstrow = 0;
    lup = 0;

    for (k = 0; k < n; k++)
    {
        /* X [k] = 0 ; */
        CLEAR(X[k]);
        Flag[k] = EMPTY;
        Lpend[k] = EMPTY; /* flag k as not pruned */
        // printf("Lpend[%d]=%d\n", k, Lpend[k]);
    }

    /* ---------------------------------------------------------------------- */
    /* mark all rows as non-pivotal and determine initial diagonal mapping */
    /* ---------------------------------------------------------------------- */

    /* PSinv does the symmetric permutation, so don't do it here */
    for (k = 0; k < n; k++)
    {
        P[k] = k;
        Pinv[k] = FLIP(k); /* mark all rows as non-pivotal */
    }
    /* initialize the construction of the off-diagonal matrix */
    Offp[0] = 0;

    /* P [k] = row means that UNFLIP (Pinv [row]) = k, and visa versa.
     * If row is pivotal, then Pinv [row] >= 0.  A row is initially "flipped"
     * (Pinv [k] < EMPTY), and then marked "unflipped" when it becomes
     * pivotal. */

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */
    // for (int i = 0; i < n; i++)
    //     printf("Pinv[%d]=%d,Stack[%d]=%d,Flag[%d]=%d,Lpend[%d]=%d,Ap_pos[%d]=%d,Llen[%d]=%d,Lip[%d]=%d\n", i, Pinv[i], i, Stack[i], i, Flag[i], i, Lpend[i], i, Ap_pos[i], i, Llen[i], i, Lip[i]);
    for (k = 0; k < n; k++)
    {

        PRINTF(("\n\n==================================== k: %d\n", k));

        /* ------------------------------------------------------------------ */
        /* determine if LU factors have grown too big */
        /* ------------------------------------------------------------------ */

        /* (n - k) entries for L and k entries for U */
        nunits = DUNITS(int, n - k) + DUNITS(int, k) +
                 DUNITS(double, n - k) + DUNITS(double, k);

        /* LU can grow by at most 'nunits' entries if the column is dense */
        PRINTF(("lup %d lusize %g lup+nunits: %g\n", lup, (double)lusize,
                lup + nunits));
        // xsize = ((double)lup) + nunits;
        // if (xsize > (double)lusize)
        // {
        //     /* check here how much to grow */
        //     xsize = (memgrow * ((double)lusize) + 4 * n + 1);
        //     if (INT_OVERFLOW(xsize))
        //     {
        //         PRINTF(("Matrix is too large (int overflow)\n"));
        //         Common->status = KLU_TOO_LARGE;
        //         return (lusize);
        //     }
        //     newlusize = memgrow * lusize + 2 * n + 1;
        //     /* Future work: retry mechanism in case of malloc failure */
        //     LU = KLU_realloc(newlusize, lusize, sizeof(double), LU, Common);
        //     Common->nrealloc++;
        //     p_LU = LU;
        //     if (Common->status == KLU_OUT_OF_MEMORY)
        //     {
        //         PRINTF(("Matrix is too large (LU)\n"));
        //         return (lusize);
        //     }
        //     lusize = newlusize;
        //     PRINTF(("inc LU to %d done\n", lusize));
        // }

        /* ------------------------------------------------------------------ */
        /* start the kth column of L and U */
        /* ------------------------------------------------------------------ */

        Lip[k] = lup;

        /* ------------------------------------------------------------------ */
        /* compute the nonzero pattern of the kth column of L and U */
        /* ------------------------------------------------------------------ */

        top = lsolve_symbolic(n, k, Ap, Ai, Q, Pinv, Stack, Flag,
                              Lpend, Ap_pos, LU, lup, Llen, Lip, k1, PSinv);
        // for (int i = 0; i < n; i++)
        //     printf("Pinv[%d]=%d,Stack[%d]=%d,Flag[%d]=%d,Lpend[%d]=%d,Ap_pos[%d]=%d,Llen[%d]=%d,Lip[%d]=%d\n", i, Pinv[i], i, Stack[i], i, Flag[i], i, Lpend[i], i, Ap_pos[i], i, Llen[i], i, Lip[i]);
        /* ------------------------------------------------------------------ */
        /* get the column of the matrix to factorize and scatter into X */
        /* ------------------------------------------------------------------ */

        construct_column(k, Ap, Ai, Ax, Q, X,
                         k1, PSinv, Rs, scale, Offp, Offi, Offx);

        /* ------------------------------------------------------------------ */
        /* compute the numerical values of the kth column (s = L \ A (:,k)) */
        /* ------------------------------------------------------------------ */

        lsolve_numeric(Pinv, LU, Stack, Lip, top, n, Llen, X);

        /* ------------------------------------------------------------------ */
        /* partial pivoting with diagonal preference */
        /* ------------------------------------------------------------------ */

        /* determine what the "diagonal" is */
        diagrow = P[k]; /* might already be pivotal */
        PRINTF(("k %d, diagrow = %d, UNFLIP (diagrow) = %d\n",
                k, diagrow, UNFLIP(diagrow)));

        /* find a pivot and scale the pivot column */
        if (!lpivot(diagrow, &pivrow, &pivot, &abs_pivot, tol, X, LU, Lip,
                    Llen, k, n, Pinv, &firstrow, Common))
        {
            /* matrix is structurally or numerically singular */
            Common->status = KLU_SINGULAR;
            if (Common->numerical_rank == EMPTY)
            {
                Common->numerical_rank = k + k1;
                Common->singular_col = Q[k + k1];
            }
            if (Common->halt_if_singular)
            {
                /* do not continue the factorization */
                return (lusize);
            }
        }

        /* we now have a valid pivot row, even if the column has NaN's or
         * has no entries on or below the diagonal at all. */
        PRINTF(("\nk %d : Pivot row %d : ", k, pivrow));
        PRINT_ENTRY(pivot);
        ASSERT(pivrow >= 0 && pivrow < n);
        ASSERT(Pinv[pivrow] < 0);

        /* set the Uip pointer */
        Uip[k] = Lip[k] + UNITS(int, Llen[k]) + UNITS(double, Llen[k]);

        /* move the lup pointer to the position where indices of U
         * should be stored */
        lup += UNITS(int, Llen[k]) + UNITS(double, Llen[k]);

        Ulen[k] = n - top;

        /* extract Stack [top..n-1] to Ui and the values to Ux and clear X */
        GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
        for (p = top, i = 0; p < n; p++, i++)
        {
            j = Stack[p];
            Ui[i] = Pinv[j];
            Ux[i] = X[j];
            CLEAR(X[j]);
        }

        /* position the lu index at the starting point for next column */
        lup += UNITS(int, Ulen[k]) + UNITS(double, Ulen[k]);

        /* U(k,k) = pivot */
        Udiag[k] = pivot;

        /* ------------------------------------------------------------------ */
        /* log the pivot permutation */
        /* ------------------------------------------------------------------ */

        ASSERT(UNFLIP(Pinv[diagrow]) < n);
        ASSERT(P[UNFLIP(Pinv[diagrow])] == diagrow);

        if (pivrow != diagrow)
        {
            /* an off-diagonal pivot has been chosen */
            Common->noffdiag++;
            PRINTF((">>>>>>>>>>>>>>>>> pivrow %d k %d off-diagonal\n",
                    pivrow, k));
            if (Pinv[diagrow] < 0)
            {
                /* the former diagonal row index, diagrow, has not yet been
                 * chosen as a pivot row.  Log this diagrow as the "diagonal"
                 * entry in the column kbar for which the chosen pivot row,
                 * pivrow, was originally logged as the "diagonal" */
                kbar = FLIP(Pinv[pivrow]);
                P[kbar] = diagrow;
                Pinv[diagrow] = FLIP(kbar);
            }
        }
        P[k] = pivrow;
        Pinv[pivrow] = k;

        /* ------------------------------------------------------------------ */
        /* symmetric pruning */
        /* ------------------------------------------------------------------ */

        prune(Lpend, Pinv, k, pivrow, LU, Uip, Lip, Ulen, Llen);

        *lnz += Llen[k] + 1; /* 1 added to lnz for diagonal */
        *unz += Ulen[k] + 1; /* 1 added to unz for diagonal */
    }

    /* ---------------------------------------------------------------------- */
    /* finalize column pointers for L and U, and put L in the pivotal order */
    /* ---------------------------------------------------------------------- */

    for (p = 0; p < n; p++)
    {
        Li = (int *)(LU + Lip[p]);
        for (i = 0; i < Llen[p]; i++)
        {
            Li[i] = Pinv[Li[i]];
        }
    }

    /* ---------------------------------------------------------------------- */
    /* shrink the LU factors to just the required size */
    /* ---------------------------------------------------------------------- */

    newlusize = lup;
    ASSERT((int)newlusize <= lusize);

    /* this cannot fail, since the block is descreasing in size */
    // LU = KLU_realloc(newlusize, lusize, sizeof(double), LU, Common);
    p_LU = LU;
    // for (int i = 0; i < newlusize; i++)
    //     printf("LU[%d]=%lf", i, LU[i]);
    return (newlusize);
}

Int klu_solve2(
    /* inputs, not modified */
    KLU_symbolic *Symbolic,
    KLU_numeric *Numeric,
    Int d,    /* leading dimension of B */
    Int nrhs, /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B[], /* size n*nrhs, in column-oriented form, with
                 * leading dimension d. */
    /* --------------- */
    KLU_common *Common)
{
    Entry x[4], offik, s;
    double rs;
    Int k1, k2, nk, k, pend, p, nr, i;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE);
    }
    if (Numeric == NULL || Symbolic == NULL || d < Symbolic->n || nrhs < 0 ||
        B == NULL)
    {
        Common->status = KLU_INVALID;
        return (FALSE);
    }
    Common->status = KLU_OK;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    int lusize_sum = 0;
    for (int i = 0, j = 0; i < Symbolic->nblocks; i++)
    {
        // printf("Numeric->LUsize[%d]=%d\n", i, Numeric->LUsize[i]);
        if (Numeric->LUsize[i])
        {
            if (j)
                lusize_sum += Numeric->LUsize[i];
            else
                j++;
        }
    }
    // printf("lusize_sum=%d\n", lusize_sum);

    // for (int i = 0; i < Numeric->LUsize[1]; i++)
    //     printf("LU[%d]=%lf\n", i, LUbx[i]);

    // for (int i = 0; i < n; i++)
    //     printf("Llen[%d]=%d,Lip[%d]=%d,Ulen[%d]=%d,Uip[%d]=%d\n", i, Llen[i], i, Lip[i], i, Ulen[i], i, Uip[i]);

    ASSERT(KLU_valid(Symbolic->n, Numeric->Offp, Numeric->Offi, Numeric->Offx));

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

    for (int chunk = 0; chunk < nrhs; chunk += 4)
    {

        /* ------------------------------------------------------------------ */
        /* get the size of the current chunk */
        /* ------------------------------------------------------------------ */

        nr = MIN(nrhs - chunk, 4);

        /* ------------------------------------------------------------------ */
        /* scale and permute the right hand side, X = P*(R\B) */
        /* ------------------------------------------------------------------ */

        if (Numeric->Rs == NULL)
        {

            /* no scaling */
            switch (nr)
            {

            case 1:

                for (k = 0; k < Symbolic->n; k++)
                {
                    Numeric->Xwork[k] = B[Numeric->Pnum[k]];
                }
                break;

            case 2:

                for (k = 0; k < Symbolic->n; k++)
                {
                    i = Numeric->Pnum[k];
                    Numeric->Xwork[2 * k] = B[i];
                    Numeric->Xwork[2 * k + 1] = B[i + d];
                }
                break;

            case 3:

                for (k = 0; k < Symbolic->n; k++)
                {
                    i = Numeric->Pnum[k];
                    Numeric->Xwork[3 * k] = B[i];
                    Numeric->Xwork[3 * k + 1] = B[i + d];
                    Numeric->Xwork[3 * k + 2] = B[i + d * 2];
                }
                break;

            case 4:

                for (k = 0; k < Symbolic->n; k++)
                {
                    i = Numeric->Pnum[k];
                    Numeric->Xwork[4 * k] = B[i];
                    Numeric->Xwork[4 * k + 1] = B[i + d];
                    Numeric->Xwork[4 * k + 2] = B[i + d * 2];
                    Numeric->Xwork[4 * k + 3] = B[i + d * 3];
                }
                break;
            }
        }
        else
        {
            switch (nr)
            {

            case 1:

                for (k = 0; k < Symbolic->n; k++)
                {
                    SCALE_DIV_ASSIGN(Numeric->Xwork[k], B[Numeric->Pnum[k]], Numeric->Rs[k]);
                }
                break;

            case 2:

                for (k = 0; k < Symbolic->n; k++)
                {
                    i = Numeric->Pnum[k];
                    rs = Numeric->Rs[k];
                    SCALE_DIV_ASSIGN(Numeric->Xwork[2 * k], B[i], rs);
                    SCALE_DIV_ASSIGN(Numeric->Xwork[2 * k + 1], B[i + d], rs);
                }
                break;

            case 3:

                for (k = 0; k < Symbolic->n; k++)
                {
                    i = Numeric->Pnum[k];
                    rs = Numeric->Rs[k];
                    SCALE_DIV_ASSIGN(Numeric->Xwork[3 * k], B[i], rs);
                    SCALE_DIV_ASSIGN(Numeric->Xwork[3 * k + 1], B[i + d], rs);
                    SCALE_DIV_ASSIGN(Numeric->Xwork[3 * k + 2], B[i + d * 2], rs);
                }
                break;

            case 4:

                for (k = 0; k < Symbolic->n; k++)
                {
                    i = Numeric->Pnum[k];
                    rs = Numeric->Rs[k];
                    SCALE_DIV_ASSIGN(Numeric->Xwork[4 * k], B[i], rs);
                    SCALE_DIV_ASSIGN(Numeric->Xwork[4 * k + 1], B[i + d], rs);
                    SCALE_DIV_ASSIGN(Numeric->Xwork[4 * k + 2], B[i + d * 2], rs);
                    SCALE_DIV_ASSIGN(Numeric->Xwork[4 * k + 3], B[i + d * 3], rs);
                }
                break;
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
                switch (nr)
                {

                case 1:
                    DIV(Numeric->Xwork[k1], Numeric->Xwork[k1], s);
                    break;

                case 2:
                    DIV(Numeric->Xwork[2 * k1], Numeric->Xwork[2 * k1], s);
                    DIV(Numeric->Xwork[2 * k1 + 1], Numeric->Xwork[2 * k1 + 1], s);
                    break;

                case 3:
                    DIV(Numeric->Xwork[3 * k1], Numeric->Xwork[3 * k1], s);
                    DIV(Numeric->Xwork[3 * k1 + 1], Numeric->Xwork[3 * k1 + 1], s);
                    DIV(Numeric->Xwork[3 * k1 + 2], Numeric->Xwork[3 * k1 + 2], s);
                    break;

                case 4:
                    DIV(Numeric->Xwork[4 * k1], Numeric->Xwork[4 * k1], s);
                    DIV(Numeric->Xwork[4 * k1 + 1], Numeric->Xwork[4 * k1 + 1], s);
                    DIV(Numeric->Xwork[4 * k1 + 2], Numeric->Xwork[4 * k1 + 2], s);
                    DIV(Numeric->Xwork[4 * k1 + 3], Numeric->Xwork[4 * k1 + 3], s);
                    break;
                }
            }
            else
            {
                KLU_lsolve(nk, Numeric->Lip + k1, Numeric->Llen + k1, &Numeric->LUbx[lusize_sum], nr,
                           Numeric->Xwork + nr * k1);
                KLU_usolve(nk, Numeric->Uip + k1, Numeric->Ulen + k1, &Numeric->LUbx[lusize_sum],
                           Numeric->Udiag + k1, nr, Numeric->Xwork + nr * k1);
                lusize_sum -= Numeric->LUsize[block];
            }

            /* -------------------------------------------------------------- */
            /* block back-substitution for the off-diagonal-block entries */
            /* -------------------------------------------------------------- */

            if (block > 0)
            {
                switch (nr)
                {

                case 1:

                    for (k = k1; k < k2; k++)
                    {
                        // printf("Offp[%d]=%d,Offp[%d]=%d\n", k, Offp[k], k + 1, Offp[k + 1]);
                        pend = Numeric->Offp[k + 1];
                        x[0] = Numeric->Xwork[k];
                        for (p = Numeric->Offp[k]; p < pend; p++)
                        {
                            // printf("X[%d]=%lf,Offx[%d]=%lf,x[0]=%lf\n", Offi[p], X[Offi[p]], p, Offx[p], x[0]);
                            MULT_SUB(Numeric->Xwork[Numeric->Offi[p]], Numeric->Offx[p], x[0]);
                        }
                    }
                    break;

                case 2:

                    for (k = k1; k < k2; k++)
                    {
                        pend = Numeric->Offp[k + 1];
                        x[0] = Numeric->Xwork[2 * k];
                        x[1] = Numeric->Xwork[2 * k + 1];
                        for (p = Numeric->Offp[k]; p < pend; p++)
                        {
                            i = Numeric->Offi[p];
                            offik = Numeric->Offx[p];
                            MULT_SUB(Numeric->Xwork[2 * i], offik, x[0]);
                            MULT_SUB(Numeric->Xwork[2 * i + 1], offik, x[1]);
                        }
                    }
                    break;

                case 3:

                    for (k = k1; k < k2; k++)
                    {
                        pend = Numeric->Offp[k + 1];
                        x[0] = Numeric->Xwork[3 * k];
                        x[1] = Numeric->Xwork[3 * k + 1];
                        x[2] = Numeric->Xwork[3 * k + 2];
                        for (p = Numeric->Offp[k]; p < pend; p++)
                        {
                            i = Numeric->Offi[p];
                            offik = Numeric->Offx[p];
                            MULT_SUB(Numeric->Xwork[3 * i], offik, x[0]);
                            MULT_SUB(Numeric->Xwork[3 * i + 1], offik, x[1]);
                            MULT_SUB(Numeric->Xwork[3 * i + 2], offik, x[2]);
                        }
                    }
                    break;

                case 4:

                    for (k = k1; k < k2; k++)
                    {
                        pend = Numeric->Offp[k + 1];
                        x[0] = Numeric->Xwork[4 * k];
                        x[1] = Numeric->Xwork[4 * k + 1];
                        x[2] = Numeric->Xwork[4 * k + 2];
                        x[3] = Numeric->Xwork[4 * k + 3];
                        for (p = Numeric->Offp[k]; p < pend; p++)
                        {
                            i = Numeric->Offi[p];
                            offik = Numeric->Offx[p];
                            MULT_SUB(Numeric->Xwork[4 * i], offik, x[0]);
                            MULT_SUB(Numeric->Xwork[4 * i + 1], offik, x[1]);
                            MULT_SUB(Numeric->Xwork[4 * i + 2], offik, x[2]);
                            MULT_SUB(Numeric->Xwork[4 * i + 3], offik, x[3]);
                        }
                    }
                    break;
                }
            }
        }

        /* ------------------------------------------------------------------ */
        /* permute the result, Bz  = Q*X */
        /* ------------------------------------------------------------------ */

        switch (nr)
        {

        case 1:

            for (k = 0; k < Symbolic->n; k++)
            {
                B[Symbolic->Q[k]] = Numeric->Xwork[k];
            }
            break;

        case 2:

            for (k = 0; k < Symbolic->n; k++)
            {
                i = Symbolic->Q[k];
                B[i] = Numeric->Xwork[2 * k];
                B[i + d] = Numeric->Xwork[2 * k + 1];
            }
            break;

        case 3:

            for (k = 0; k < Symbolic->n; k++)
            {
                i = Symbolic->Q[k];
                B[i] = Numeric->Xwork[3 * k];
                B[i + d] = Numeric->Xwork[3 * k + 1];
                B[i + d * 2] = Numeric->Xwork[3 * k + 2];
            }
            break;

        case 4:

            for (k = 0; k < Symbolic->n; k++)
            {
                i = Symbolic->Q[k];
                B[i] = Numeric->Xwork[4 * k];
                B[i + d] = Numeric->Xwork[4 * k + 1];
                B[i + d * 2] = Numeric->Xwork[4 * k + 2];
                B[i + d * 3] = Numeric->Xwork[4 * k + 3];
            }
            break;
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */

        B += d * 4;
    }
    return (TRUE);
}

void KLU_lsolve(
    /* inputs, not modified: */
    Int n,
    Int Lip[],
    Int Llen[],
    Unit LU[],
    Int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    Entry X[])
{
    Entry x[4], lik;
    Int *Li;
    Entry *Lx;
    Int k, p, len, i;

    switch (nrhs)
    {

    case 1:
        for (k = 0; k < n; k++)
        {
            x[0] = X[k];
            GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
            /* unit diagonal of L is not stored*/
            for (p = 0; p < len; p++)
            {
                /* X [Li [p]] -= Lx [p] * x [0] ; */
                MULT_SUB(X[Li[p]], Lx[p], x[0]);
                // printf("X[%d]=%lf,Lx[%d]=%lf,x[0]=%lf\n", Li[p], X[Li[p]], p, Lx[p], x[0]);
            }
        }
        break;

    case 2:

        for (k = 0; k < n; k++)
        {
            x[0] = X[2 * k];
            x[1] = X[2 * k + 1];
            GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
            for (p = 0; p < len; p++)
            {
                i = Li[p];
                lik = Lx[p];
                MULT_SUB(X[2 * i], lik, x[0]);
                MULT_SUB(X[2 * i + 1], lik, x[1]);
            }
        }
        break;

    case 3:

        for (k = 0; k < n; k++)
        {
            x[0] = X[3 * k];
            x[1] = X[3 * k + 1];
            x[2] = X[3 * k + 2];
            GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
            for (p = 0; p < len; p++)
            {
                i = Li[p];
                lik = Lx[p];
                MULT_SUB(X[3 * i], lik, x[0]);
                MULT_SUB(X[3 * i + 1], lik, x[1]);
                MULT_SUB(X[3 * i + 2], lik, x[2]);
            }
        }
        break;

    case 4:

        for (k = 0; k < n; k++)
        {
            x[0] = X[4 * k];
            x[1] = X[4 * k + 1];
            x[2] = X[4 * k + 2];
            x[3] = X[4 * k + 3];
            GET_POINTER(LU, Lip, Llen, Li, Lx, k, len);
            for (p = 0; p < len; p++)
            {
                i = Li[p];
                lik = Lx[p];
                MULT_SUB(X[4 * i], lik, x[0]);
                MULT_SUB(X[4 * i + 1], lik, x[1]);
                MULT_SUB(X[4 * i + 2], lik, x[2]);
                MULT_SUB(X[4 * i + 3], lik, x[3]);
            }
        }
        break;
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
    Int n,
    Int Uip[],
    Int Ulen[],
    Unit LU[],
    Entry Udiag[],
    Int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    Entry X[])
{
    Entry x[4], uik, ukk;
    Int *Ui;
    Entry *Ux;
    Int k, p, len, i;

    switch (nrhs)
    {

    case 1:

        for (k = n - 1; k >= 0; k--)
        {
            GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
            /* x [0] = X [k] / Udiag [k] ; */
            DIV(x[0], X[k], Udiag[k]);
            X[k] = x[0];
            for (p = 0; p < len; p++)
            {
                /* X [Ui [p]] -= Ux [p] * x [0] ; */
                MULT_SUB(X[Ui[p]], Ux[p], x[0]);
            }
        }

        break;

    case 2:

        for (k = n - 1; k >= 0; k--)
        {
            GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
            ukk = Udiag[k];
            /* x [0] = X [2*k    ] / ukk ;
            x [1] = X [2*k + 1] / ukk ; */
            DIV(x[0], X[2 * k], ukk);
            DIV(x[1], X[2 * k + 1], ukk);

            X[2 * k] = x[0];
            X[2 * k + 1] = x[1];
            for (p = 0; p < len; p++)
            {
                i = Ui[p];
                uik = Ux[p];
                /* X [2*i    ] -= uik * x [0] ;
                X [2*i + 1] -= uik * x [1] ; */
                MULT_SUB(X[2 * i], uik, x[0]);
                MULT_SUB(X[2 * i + 1], uik, x[1]);
            }
        }

        break;

    case 3:

        for (k = n - 1; k >= 0; k--)
        {
            GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
            ukk = Udiag[k];

            DIV(x[0], X[3 * k], ukk);
            DIV(x[1], X[3 * k + 1], ukk);
            DIV(x[2], X[3 * k + 2], ukk);

            X[3 * k] = x[0];
            X[3 * k + 1] = x[1];
            X[3 * k + 2] = x[2];
            for (p = 0; p < len; p++)
            {
                i = Ui[p];
                uik = Ux[p];
                MULT_SUB(X[3 * i], uik, x[0]);
                MULT_SUB(X[3 * i + 1], uik, x[1]);
                MULT_SUB(X[3 * i + 2], uik, x[2]);
            }
        }

        break;

    case 4:

        for (k = n - 1; k >= 0; k--)
        {
            GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, len);
            ukk = Udiag[k];

            DIV(x[0], X[4 * k], ukk);
            DIV(x[1], X[4 * k + 1], ukk);
            DIV(x[2], X[4 * k + 2], ukk);
            DIV(x[3], X[4 * k + 3], ukk);

            X[4 * k] = x[0];
            X[4 * k + 1] = x[1];
            X[4 * k + 2] = x[2];
            X[4 * k + 3] = x[3];
            for (p = 0; p < len; p++)
            {
                i = Ui[p];
                uik = Ux[p];

                MULT_SUB(X[4 * i], uik, x[0]);
                MULT_SUB(X[4 * i + 1], uik, x[1]);
                MULT_SUB(X[4 * i + 2], uik, x[2]);
                MULT_SUB(X[4 * i + 3], uik, x[3]);
            }
        }

        break;
    }
}

Int KLU_scale /* return TRUE if successful, FALSE otherwise */
    (
        /* inputs, not modified */
        Int scale, /* 0: none, 1: sum, 2: max */
        Int n,
        Int Ap[], /* size n+1, column pointers */
        Int Ai[], /* size nz, row indices */
        double Ax[],
        /* outputs, not defined on input */
        double Rs[], /* size n, can be NULL if scale <= 0 */
        /* workspace, not defined on input or output */
        Int W[], /* size n, can be NULL */
        /* --------------- */
        KLU_common *Common)
{
    double a;
    Entry *Az;
    Int row, col, p, pend, check_duplicates;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE);
    }
    Common->status = KLU_OK;

    if (scale < 0)
    {
        /* return without checking anything and without computing the
         * scale factors */
        return (TRUE);
    }

    Az = (Entry *)Ax;

    if (n <= 0 || Ap == NULL || Ai == NULL || Az == NULL ||
        (scale > 0 && Rs == NULL))
    {
        /* Ap, Ai, Ax and Rs must be present, and n must be > 0 */
        Common->status = KLU_INVALID;
        return (FALSE);
    }
    if (Ap[0] != 0 || Ap[n] < 0)
    {
        /* nz = Ap [n] must be >= 0 and Ap [0] must equal zero */
        Common->status = KLU_INVALID;
        return (FALSE);
    }
    for (col = 0; col < n; col++)
    {
        if (Ap[col] > Ap[col + 1])
        {
            /* column pointers must be non-decreasing */
            Common->status = KLU_INVALID;
            return (FALSE);
        }
    }

    /* ---------------------------------------------------------------------- */
    /* scale */
    /* ---------------------------------------------------------------------- */

    if (scale > 0)
    {
        /* initialize row sum or row max */
        for (row = 0; row < n; row++)
        {
            Rs[row] = 0;
        }
    }

    /* check for duplicates only if W is present */
    check_duplicates = (W != (Int *)NULL);
    if (check_duplicates)
    {
        for (row = 0; row < n; row++)
        {
            W[row] = EMPTY;
        }
    }

    for (col = 0; col < n; col++)
    {
        pend = Ap[col + 1];
        for (p = Ap[col]; p < pend; p++)
        {
            row = Ai[p];
            if (row < 0 || row >= n)
            {
                /* row index out of range, or duplicate entry */
                Common->status = KLU_INVALID;
                return (FALSE);
            }
            if (check_duplicates)
            {
                if (W[row] == col)
                {
                    /* duplicate entry */
                    Common->status = KLU_INVALID;
                    return (FALSE);
                }
                /* flag row i as appearing in column col */
                W[row] = col;
            }
            /* a = ABS (Az [p]) ;*/
            ABS(a, Az[p]);
            if (scale == 1)
            {
                /* accumulate the abs. row sum */
                Rs[row] += a;
            }
            else if (scale > 1)
            {
                /* find the max abs. value in the row */
                Rs[row] = MAX(Rs[row], a);
            }
        }
    }

    if (scale > 0)
    {
        /* do not scale empty rows */
        for (row = 0; row < n; row++)
        {
            /* matrix is singular */
            PRINTF(("Rs [%d] = %g\n", row, Rs[row]));

            if (Rs[row] == 0.0)
            {
                PRINTF(("Row %d of A is all zero\n", row));
                Rs[row] = 1.0;
            }
        }
    }

    return (TRUE);
}

static void factor2(
    /* inputs, not modified */
    Int Ap[], /* size n+1, column pointers */
    Int Ai[], /* size nz, row indices */
    Entry Ax[],
    KLU_symbolic *Symbolic,

    /* inputs, modified on output: */
    KLU_numeric *Numeric,
    KLU_common *Common)
{
    Int *Pblock, k1, k2, nk, oldcol, pend, oldrow, lnz = 0, unz = 0, newrow, poff, max_lnz_block = 1, max_unz_block = 1, lusize_sum = 0;
    double s;
    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    // Iwork = Numeric->Iwork; /* 5*maxblock for KLU_factor */
    //                         /* 1*maxblock for Pblock */
    Pblock = Numeric->Iwork + 5 * ((int)Symbolic->maxblock);

    /* compute the inverse of P from symbolic analysis.  Will be updated to
     * become the inverse of the numerical factorization when the factorization
     * is done, for use in KLU_refactor */

    for (int k = 0; k < Symbolic->n; k++)
    {
        ASSERT(Symbolic->P[k] >= 0 && Symbolic->P[k] < Symbolic->n);
        Numeric->Pinv[Symbolic->P[k]] = k;
    }

    Common->noffdiag = 0;
    Numeric->Offp[0] = 0;

    /* ---------------------------------------------------------------------- */
    /* optionally check input matrix and compute scale factors */
    /* ---------------------------------------------------------------------- */

    if (Common->scale >= 0)
    {
        /* use Pnum as workspace. NOTE: scale factors are not yet permuted
         * according to the final pivot row ordering, so Rs [oldrow] is the
         * scale factor for A (oldrow,:), for the user's matrix A.  Pnum is
         * used as workspace in KLU_scale.  When the factorization is done,
         * the scale factors are permuted according to the final pivot row
         * permutation, so that Rs [k] is the scale factor for the kth row of
         * A(p,q) where p and q are the final row and column permutations. */
        KLU_scale(Common->scale, Symbolic->n, Ap, Ai, (double *)Ax, Numeric->Rs, Numeric->Pnum, Common);
        if (Common->status < KLU_OK)
        {
            /* matrix is invalid */
            return;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* factor each block using klu */
    /* ---------------------------------------------------------------------- */

    for (int block = 0; block < Symbolic->nblocks; block++)
    {

        /* ------------------------------------------------------------------ */
        /* the block is from rows/columns k1 to k2-1 */
        /* ------------------------------------------------------------------ */

        k1 = Symbolic->R[block];
        k2 = Symbolic->R[block + 1];
        nk = k2 - k1;
        PRINTF(("FACTOR BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2 - 1, nk));

        if (nk == 1)
        {

            /* -------------------------------------------------------------- */
            /* singleton case */
            /* -------------------------------------------------------------- */

            poff = Numeric->Offp[k1];
            oldcol = Symbolic->Q[k1];
            pend = Ap[oldcol + 1];
            CLEAR(s);

            if (Common->scale <= 0)
            {
                /* no scaling */
                for (int p = Ap[oldcol]; p < pend; p++)
                {
                    oldrow = Ai[p];
                    newrow = Numeric->Pinv[oldrow];
                    if (newrow < k1)
                    {
                        Numeric->Offi[poff] = oldrow;
                        Numeric->Offx[poff] = Ax[p];
                        poff++;
                    }
                    else
                    {
                        ASSERT(newrow == k1);
                        PRINTF(("singleton block %d", block));
                        PRINT_ENTRY(Ax[p]);
                        s = Ax[p];
                    }
                }
            }
            else
            {
                /* row scaling.  NOTE: scale factors are not yet permuted
                 * according to the pivot row permutation, so Rs [oldrow] is
                 * used below.  When the factorization is done, the scale
                 * factors are permuted, so that Rs [newrow] will be used in
                 * klu_solve, klu_tsolve, and klu_rgrowth */
                for (int p = Ap[oldcol]; p < pend; p++)
                {
                    oldrow = Ai[p];
                    newrow = Numeric->Pinv[oldrow];
                    if (newrow < k1)
                    {
                        Numeric->Offi[poff] = oldrow;
                        /* Offx [poff] = Ax [p] / Rs [oldrow] ; */
                        SCALE_DIV_ASSIGN(Numeric->Offx[poff], Ax[p], Numeric->Rs[oldrow]);
                        poff++;
                    }
                    else
                    {
                        ASSERT(newrow == k1);
                        PRINTF(("singleton block %d ", block));
                        PRINT_ENTRY(Ax[p]);
                        SCALE_DIV_ASSIGN(s, Ax[p], Numeric->Rs[oldrow]);
                    }
                }
            }

            Numeric->Udiag[k1] = s;

            if (IS_ZERO(s))
            {
                /* singular singleton */
                Common->status = KLU_SINGULAR;
                Common->numerical_rank = k1;
                Common->singular_col = oldcol;
                if (Common->halt_if_singular)
                {
                    return;
                }
            }

            Numeric->Offp[k1 + 1] = poff;
            Numeric->Pnum[k1] = Symbolic->P[k1];
            lnz++;
            unz++;
        }
        else
        {
            double lsize;
            /* -------------------------------------------------------------- */
            /* construct and factorize the kth block */
            /* -------------------------------------------------------------- */

            if (Symbolic->Lnz[block] < 0)
            {
                /* COLAMD was used - no estimate of fill-in */
                /* use 10 times the nnz in A, plus n */
                lsize = -(Common->initmem);
            }
            else
            {
                lsize = Common->initmem_amd * Symbolic->Lnz[block] + nk;
            }

            nk = MAX(1, nk);

            int anz = Ap[nk + k1] - Ap[k1];

            int Lsize, maxlnz;

            if (lsize <= 0)
            {
                lsize = -lsize;
                lsize = MAX(lsize, 1.0);
                Lsize = lsize * anz + nk;
            }
            else
            {
                Lsize = lsize;
            }

            int Usize = Lsize;

            Lsize = MAX(nk + 1, Lsize);
            Usize = MAX(nk + 1, Usize);

            maxlnz = (((double)nk) * ((double)nk) + ((double)nk)) / 2.;
            Lsize = MIN(maxlnz, Lsize);
            Usize = MIN(maxlnz, Usize);

            int lusize = DUNITS(Int, Lsize) + DUNITS(Entry, Lsize) +
                         DUNITS(Int, Usize) + DUNITS(Entry, Usize);

            int *pinv, *Stack, *Flag, *Ap_pos, *Lpend, *W, lnz_block, unz_block;
            W = Numeric->Iwork;
            pinv = (Int *)W;
            W += nk;
            Stack = (Int *)W;
            W += nk;
            Flag = (Int *)W;
            W += nk;
            Lpend = (Int *)W;
            W += nk;
            Ap_pos = (Int *)W;
            W += nk;

            Numeric->LUsize[block] = KLU_kernel(nk, Ap, Ai, Ax, Symbolic->Q, lusize,
                                                pinv, Pblock, &Numeric->LUbx[lusize_sum], Numeric->Udiag + k1, Numeric->Llen + k1, Numeric->Ulen + k1, Numeric->Lip + k1, Numeric->Uip + k1, &lnz_block, &unz_block,
                                                Numeric->Xwork, Stack, Flag, Ap_pos, Lpend,
                                                k1, Numeric->Pinv, Numeric->Rs, Numeric->Offp, Numeric->Offi, Numeric->Offx, Common);
            if (Common->status < KLU_OK || (Common->status == KLU_SINGULAR && Common->halt_if_singular))
            {
                /* out of memory, invalid inputs, or singular */
                return;
            }

            lusize_sum += Numeric->LUsize[block];

            PRINTF(("\n----------------------- L %d:\n", block));
            ASSERT(KLU_valid_LU(nk, TRUE, Numeric->Lip + k1, Numeric->Llen + k1, Numeric->LUbx[block]));
            PRINTF(("\n----------------------- U %d:\n", block));
            ASSERT(KLU_valid_LU(nk, FALSE, Numeric->Uip + k1, Numeric->Ulen + k1, Numeric->LUbx[block]));

            /* -------------------------------------------------------------- */
            /* get statistics */
            /* -------------------------------------------------------------- */

            lnz += lnz_block;
            unz += unz_block;
            max_lnz_block = MAX(max_lnz_block, lnz_block);
            max_unz_block = MAX(max_unz_block, unz_block);

            if (Symbolic->Lnz[block] == EMPTY)
            {
                /* revise estimate for subsequent factorization */
                Symbolic->Lnz[block] = MAX(lnz_block, unz_block);
            }

            /* -------------------------------------------------------------- */
            /* combine the klu row ordering with the symbolic pre-ordering */
            /* -------------------------------------------------------------- */

            PRINTF(("Pnum, 1-based:\n"));
            for (int k = 0; k < nk; k++)
            {
                ASSERT(k + k1 < Symbolic->n);
                ASSERT(Pblock[k] + k1 < Symbolic->n);
                // printf("P[%d]=%d\n", k, P[k]);
                Numeric->Pnum[k + k1] = Symbolic->P[Pblock[k] + k1];
                PRINTF(("Pnum (%d + %d + 1 = %d) = %d + 1 = %d\n",
                        k, k1, k + k1 + 1, Numeric->Pnum[k + k1], Numeric->Pnum[k + k1] + 1));
            }

            /* the local pivot row permutation Pblock is no longer needed */
        }
    }
    ASSERT(Symbolic->nzoff == Numeric->Offp[Symbolic->n]);
    PRINTF(("\n------------------- Off diagonal entries:\n"));
    ASSERT(KLU_valid(Symbolic->n, Numeric->Offp, Numeric->Offi, Numeric->Offx));

    Numeric->lnz = lnz;
    Numeric->unz = unz;
    Numeric->max_lnz_block = max_lnz_block;
    Numeric->max_unz_block = max_unz_block;

    /* compute the inverse of Pnum */
    for (int k = 0; k < Symbolic->n; k++)
    {
        ASSERT(Numeric->Pnum[k] >= 0 && Numeric->Pnum[k] < Symbolic->n);
        Numeric->Pinv[Numeric->Pnum[k]] = k;
    }

    /* permute scale factors Rs according to pivotal row order */
    if (Common->scale > 0)
    {
        for (int k = 0; k < Symbolic->n; k++)
        {
            REAL(Numeric->Xwork[k]) = Numeric->Rs[Numeric->Pnum[k]];
        }
        for (int k = 0; k < Symbolic->n; k++)
        {
            Numeric->Rs[k] = REAL(Numeric->Xwork[k]);
        }
    }

    PRINTF(("\n------------------- Off diagonal entries, old:\n"));
    ASSERT(KLU_valid(Symbolic->n, Numeric->Offp, Numeric->Offi, Numeric->Offx));

    /* apply the pivot row permutations to the off-diagonal entries */
    for (int p = 0; p < Symbolic->nzoff; p++)
    {
        ASSERT(Numeric->Offi[p] >= 0 && Numeric->Offi[p] < Symbolic->n);
        Numeric->Offi[p] = Numeric->Pinv[Numeric->Offi[p]];
    }

    PRINTF(("\n------------------- Off diagonal entries, new:\n"));
    ASSERT(KLU_valid(Symbolic->n, Numeric->Offp, Numeric->Offi, Numeric->Offx));
}

int main(void)
{
    Mtx mtx;
    char filename[] = "../../Matrix_Sample/arc130.mtx";

    if (read_sparse(filename, &mtx))
        return 1;

    klu_common Common;
    KLU_numeric Numeric;
    klu_symbolic Symbolic;
    klu_defaults(&Common);
    // const int n = 10;
    // int Ap[] = {0, 2, 4, 7, 9, 13, 14, 18, 21, 23, 24};
    // int Ai[] = {0, 7, 1, 5, 0, 2, 7, 3, 8, 3, 4, 5, 8, 5, 0, 6, 7, 8, 4, 5, 7, 2, 8, 9};
    // double Ax[] = {8, 8, 2, 4, 10, 3, 10, 1, 5, 2, 1, 4, 10, 2, 3, 5, 3, 15, 4, 16, 3, 7, 2, 9};
    // double b[] = {172, 18, 38, 19, 18, 118, 20, 181, 159, 9};

    int n = mtx.n;
    int *Ap = mtx.Ap;
    int *Ai = mtx.Ai;
    double *Ax = mtx.Ax;
    double b[] = {279.584743202215, -324.558579717574, 41.6857698124613, 71.7796552191998, -21.8080536390688, -12.1206651085711, 6.96409248369851, 8.06655960258013, 9.07550956120881, 10.0242759948783, 10.9356800700396, 11.2575737075623, 13.0755095612088, 14.6143625004934, 14.6622950007692, 16, 17.0665596025801, 18.0011830061659, 19.0153260491416, -27475.3902422406, -80429157.8945312, -79328370.8867188, -40946052.6914062, -74347681.515625, -72174052.6171875, 45.2518649101257, 25.4358819499556, 29.3744153231358, 33.4031240940092, 32.1007794141767, 68.682388305664, 27.5902938740802, 34.2920764182203, 41.8601366281505, 38.4772992134089, 85.2251358032226, 29.9290522684125, 38.2619187734309, 48.8282384276385, 44.3459486961358, 91.8335390090942, 33.3837329002542, 41.6764824338056, 53.9081966876978, 49.5371791720383, 89.9676032066345, 38.418592582108, 45.5973847664851, 57.5117719173426, 54.1664958000176, 83.7884101867675, 44.8544080049303, 50.2856411694522, 60.4747216701503, 58.5157126188272, 77.5720725059509, 52.0548983215017, 55.7535365455683, 63.5658667087552, 62.8661473095413, 73.8493753075599, 59.3249367251971, 61.7303034401082, 67.1980535984037, 67.3882255144415, 73.2623927593231, 66.1970339380203, 67.8179687568906, 71.427210163325, 72.1252968907355, 75.2432642094791, 72.5106880366802, 73.7411950565874, 76.1072703003883, 77.0394964143633, 78.8693265914917, 78.3303363509476, 79.4039598405361, 81.0516444817185, 82.06717222929, 83.3696878105402, 83.8095736652612, 84.8357935622334, 86.1200445592403, 87.1544153615832, 88.2661053612828, 89.093873038888, 90.1123360395431, 91.2343204244971, 92.2676514834165, 93.3175967484713, 94.2826037406922, 95.3031234964728, 96.3608182370663, 97.390259783715, 98.4220633506775, 99.4306319057942, 100.453658305109, 101.48865038529, 102.515559270978, 103.54242381081, 104.563432127237, 105.587807916105, 106.615641176701, 107.641398645937, 108.666984058917, 109.691199589521, 110.716114014387, 111.741927120835, 112.767261900008, 113.792525097728, 114.817488133907, 115.842579871416, 116.867881961167, 117.893085516989, 118.918268278241, 119.943389151245, 120.968530938029, 121.993717238307, 123.01888525486, 124.044048041105, 125.069199554622, 126.094353720546, 127.119516149163, 128.144675865769, 129.169834211469, 130.194990679622, 131.220147609711, 132.245305974036, 133.270463384688};

    // printf("n=%d\n", n);
    // for (int i = 0; i <= n; i++)
    //     printf("Ap[%d]=%d\n", i, mtx.Ap[i]);
    // for (int i = 0; i < mtx.nz; i++)
    //     printf("Ai[%d]=%d\tAx[%d]=%lf\n", i, mtx.Ai[i], i, mtx.Ax[i]);

    Symbolic = *klu_analyze(n, Ap, Ai, &Common);

    // for (int i = 0; i < n; i++)
    //     printf("P[%d]=%d,Q[%d]=%d,R[%d]=%d,Lnz[%d]=%lf\n", i, Symbolic.P[i], i, Symbolic.Q[i], i, Symbolic.R[i], i, Symbolic.Lnz[i]);
    // printf("nblocks=%d,nzoff=%d\n", Symbolic.nblocks, Symbolic.nzoff);

    int nzoff1 = Symbolic.nzoff + 1,
        n1 = n + 1;
    double lusize = Common.memgrow * (Symbolic.lnz + Symbolic.unz) + 4 * n + 1;

    printf("lusize=%.0lf\n", lusize);

    Numeric.n = Symbolic.n;
    Numeric.nblocks = Symbolic.nblocks;
    Numeric.nzoff = Symbolic.nzoff;
    Numeric.Pnum = malloc(n * sizeof(Int));
    Numeric.Offp = malloc(n1 * sizeof(Int));
    Numeric.Offi = malloc(nzoff1 * sizeof(Int));
    Numeric.Offx = malloc(nzoff1 * sizeof(Entry));
    Numeric.Lip = malloc(n * sizeof(Int));
    Numeric.Uip = malloc(n * sizeof(Int));
    Numeric.Llen = malloc(n * sizeof(Int));
    Numeric.Ulen = malloc(n * sizeof(Int));
    Numeric.LUsize = calloc(Symbolic.nblocks, sizeof(int));
    Numeric.LUbx = malloc(lusize * sizeof(double));
    Numeric.Udiag = malloc(n * sizeof(Entry));
    Numeric.Rs = malloc(n * sizeof(double));
    Numeric.Pinv = malloc(n * sizeof(Int));

    // int ok = 1;
    // int s = KLU_mult_size_t(n, sizeof(Entry), &ok);
    // int n3 = KLU_mult_size_t(n, 3 * sizeof(Entry), &ok);
    // int b6 = KLU_mult_size_t(maxblock, 6 * sizeof(Int), &ok);
    // Numeric.worksize = KLU_add_size_t(s, MAX(n3, b6), &ok);

    // int s = n * sizeof(double);
    // int n3 = n * 3 * sizeof(double);
    // int b6 = Symbolic.maxblock * 6 * sizeof(int);
    Numeric.worksize = n * sizeof(double) + MAX(n * 3 * sizeof(double), Symbolic.maxblock * 6 * sizeof(int));

    // Numeric.Work = malloc(Numeric.worksize);
    Numeric.Xwork = malloc(Numeric.worksize);
    Numeric.Iwork = (int *)(Numeric.Xwork + n);

    factor2(Ap, Ai, Ax, &Symbolic, &Numeric, &Common);

    klu_solve2(&Symbolic, &Numeric, n, 1, b, &Common);

    for (int i = 0; i < n; i++)
        printf("x [%d] = %g\n", i, b[i]);
    return 0;
}