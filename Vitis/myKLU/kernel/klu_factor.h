#include "klu_kernel.h"

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
    int pos, l_length = *plength;
    int *Li;

    Stack[0] = j;
    ASSERT(Flag[j] != k);

dfs_loop:
    for (int head = 0; head >= 0;)
    {
        j = Stack[head];
        int jnew = Pinv[j];
        ASSERT(jnew >= 0 && jnew < k); /* j is pivotal */

        if (Flag[j] != k) /* a node is not yet visited */
        {
            /* first time that j has been visited */
            Flag[j] = k;
            PRINTF(("[ start dfs at %d : new %d\n", j, jnew));
            /* set Ap_pos [head] to one past the last entry in col j to scan */
            Ap_pos[head] = (Lpend[jnew] == EMPTY) ? Llen[jnew] : Lpend[jnew];
        }

        /* add the adjacent nodes to the recursive stack by iterating through
         * until finding another non-visited pivotal node */
        Li = (int *)(LU + Lip[jnew]);

    dfs_loop_2:
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
    int pend, oldcol, kglobal, top = n, l_length = 0;

    Lik = (int *)(LU + lup);

    /* ---------------------------------------------------------------------- */
    /* BTF factorization of A (k1:k2-1, k1:k2-1) */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1;    /* column k of the block is col kglobal of A */
    oldcol = Q[kglobal]; /* Q must be present for BTF case */
    pend = Ap[oldcol + 1];

lsolve_symbolic_loop:
    for (int p = Ap[oldcol]; p < pend; p++)
    {
        int i = PSinv[Ai[p]] - k1;
        if (i < 0)
            continue; /* skip entry outside the block */

        /* (i,k) is an entry in the block.  start a DFS at node i */
        PRINTF(("\n ===== DFS at node %d in b, inew: %d\n", i, Pinv[i]));
        if (Flag[i] != k)
        {
            if (Pinv[i] >= 0)
            {
                top = dfs(i, k, Pinv, Llen, Lip, Stack, Flag, Lpend, top, LU, Lik, &l_length, Ap_pos);
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
    int pend, oldcol, kglobal, poff, oldrow;

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
        for (int p = Ap[oldcol]; p < pend; p++)
        {
            oldrow = Ai[p];
            int i = PSinv[oldrow] - k1;
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
        for (int p = Ap[oldcol]; p < pend; p++)
        {
            oldrow = Ai[p];
            int i = PSinv[oldrow] - k1;
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
    double *Lx;
    int *Li, len;

    /* solve Lx=b */
lsolve_numeric_loop:
    for (int s = top; s < n; s++)
    {
        /* forward solve with column j of L */
        int j = Stack[s];
        int jnew = Pinv[j];
        ASSERT(jnew >= 0);
        double xj = X[j];
        GET_POINTER(LU, Lip, Llen, Li, Lx, jnew, len);
        ASSERT(Lip[jnew] <= Lip[jnew + 1]);
        for (int p = 0; p < len; p++)
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
    int i, ppivrow, pdiag, pivrow, *Li, last_row_index, firstrow, len;

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
lpivot_loop:
    for (int p = 0; p < len; p++)
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
    for (int p = 0; p < Llen[k]; p++)
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
    double *Lx, *Ux;
    int *Li, *Ui;
    int llen, ulen;

    /* check to see if any column of L can be pruned */
    /* Ux is set but not used.  This OK. */
    GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, ulen);

prune_loop:
    for (int p = 0; p < ulen; p++)
    {
        int j = Ui[p];
        ASSERT(j < k);
        PRINTF(("%d is pruned: %d. Lpend[j] %d Lip[j+1] %d\n",
                j, Lpend[j] != EMPTY, Lpend[j], Lip[j + 1]));
        if (Lpend[j] == EMPTY)
        {
            /* scan column j of L for the pivot row */
            GET_POINTER(LU, Lip, Llen, Li, Lx, j, llen);
            for (int p2 = 0; p2 < llen; p2++)
            {
                if (pivrow == Li[p2])
                {
                    /* found it!  This column can be pruned */

                    /* partition column j of L.  The unit diagonal of L
                     * is not stored in the column of L. */
                    int phead = 0;
                    int ptail = Llen[j];
                    while (phead < ptail)
                    {
                        int i = Li[phead];
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
                            double x = Lx[phead];
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

/* Scale a matrix and check to see if it is valid.  Can be called by the user.
 * This is called by KLU_factor and KLU_refactor.  Returns TRUE if the input
 * matrix is valid, FALSE otherwise.  If the W input argument is non-NULL,
 * then the input matrix is checked for duplicate entries.
 *
 * scaling methods:
 *      <0: no scaling, do not compute Rs, and do not check input matrix.
 *      0: no scaling
 *      1: the scale factor for row i is sum (abs (A (i,:)))
 *      2 or more: the scale factor for row i is max (abs (A (i,:)))
 */
static int KLU_scale /* return TRUE if successful, FALSE otherwise */
    (
        /* inputs, not modified */
        int scale, /* 0: none, 1: sum, 2: max */
        int n,
        int Ap[], /* size n+1, column pointers */
        int Ai[], /* size nz, row indices */
        double Ax[],
        /* outputs, not defined on input */
        double Rs[], /* size n, can be NULL if scale <= 0 */
        /* workspace, not defined on input or output */
        int W[], /* size n, can be NULL */
        /* --------------- */
        KLU_common *Common)
{
    if (scale < 0)
    {
        /* return without checking anything and without computing the
         * scale factors */
        return (TRUE);
    }

    /* ---------------------------------------------------------------------- */
    /* scale */
    /* ---------------------------------------------------------------------- */

    if (scale > 0)
    {
        /* initialize row sum or row max */
        for (int row = 0; row < n; row++)
        {
            Rs[row] = 0;
        }
    }

    /* check for duplicates only if W is present */
    int check_duplicates = (W != (int *)NULL);
    if (check_duplicates)
    {
        for (int row = 0; row < n; row++)
        {
            W[row] = EMPTY;
        }
    }

    for (int col = 0; col < n; col++)
    {
        int pend = Ap[col + 1];
        for (int p = Ap[col]; p < pend; p++)
        {
            int row = Ai[p];
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
            double a;
            ABS(a, Ax[p]);
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
        for (int row = 0; row < n; row++)
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

/* ========================================================================== */
/* === KLU_kernel =========================================================== */
/* ========================================================================== */

static int KLU_kernel /* final size of LU on output */
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
    double pivot = 0, abs_pivot, xsize, nunits;
    double *Ux;
    int *Li, *Ui;
    double *LU; /* LU factors (pattern and values) */
    int pivrow = 0, diagrow, firstrow = 0, lup = 0, len, newlusize;

    ASSERT(Common != NULL);
    *lnz = 0;
    *unz = 0;

    /* ---------------------------------------------------------------------- */
    /* get initial Li, Lx, Ui, and Ux */
    /* ---------------------------------------------------------------------- */

    PRINTF(("input: lusize %d \n", lusize));
    ASSERT(lusize > 0);
    LU = p_LU;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    for (int k = 0; k < n; k++)
    {
        /* X [k] = 0 ; */
        CLEAR(X[k]);
        Flag[k] = EMPTY;
        Lpend[k] = EMPTY; /* flag k as not pruned */
    }

    /* ---------------------------------------------------------------------- */
    /* mark all rows as non-pivotal and determine initial diagonal mapping */
    /* ---------------------------------------------------------------------- */

    /* PSinv does the symmetric permutation, so don't do it here */
    for (int k = 0; k < n; k++)
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
klu_kernel_factor_loop:
    for (int k = 0; k < n; k++)
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

        /* ------------------------------------------------------------------ */
        /* start the kth column of L and U */
        /* ------------------------------------------------------------------ */

        Lip[k] = lup;

        /* ------------------------------------------------------------------ */
        /* compute the nonzero pattern of the kth column of L and U */
        /* ------------------------------------------------------------------ */

        int top = lsolve_symbolic(n, k, Ap, Ai, Q, Pinv, Stack, Flag, Lpend, Ap_pos, LU, lup, Llen, Lip, k1, PSinv);
        /* ------------------------------------------------------------------ */
        /* get the column of the matrix to factorize and scatter into X */
        /* ------------------------------------------------------------------ */

        construct_column(k, Ap, Ai, Ax, Q, X, k1, PSinv, Rs, Common->scale, Offp, Offi, Offx);

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
        if (!lpivot(diagrow, &pivrow, &pivot, &abs_pivot, Common->tol, X, LU, Lip, Llen, k, n, Pinv, &firstrow, Common))
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
        for (int p = top, i = 0; p < n; p++, i++)
        {
            int j = Stack[p];
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
                int kbar = FLIP(Pinv[pivrow]);
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

    for (int p = 0; p < n; p++)
    {
        Li = (int *)(LU + Lip[p]);
        for (int i = 0; i < Llen[p]; i++)
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

    return (newlusize);
}

void KLU_factor(
    /* inputs, not modified */
    int Ap[], /* size n+1, column pointers */
    int Ai[], /* size nz, row indices */
    double Ax[],
    KLU_symbolic *Symbolic,

    /* inputs, modified on output: */
    KLU_numeric *Numeric,
    KLU_common *Common)
{
    int oldcol, pend, oldrow, lnz = 0, unz = 0, newrow, poff, max_lnz_block = 1, max_unz_block = 1;

    Numeric->lusize_sum = 0;
    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    // Iwork = Numeric->Iwork; /* 5*maxblock for KLU_factor */
    //                         /* 1*maxblock for Pblock */

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
        KLU_scale(Common->scale, Symbolic->n, Ap, Ai, Ax, Numeric->Rs, Numeric->Pnum, Common);
        if (Common->status < KLU_OK)
        {
            /* matrix is invalid */
            return;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* factor each block using klu */
    /* ---------------------------------------------------------------------- */
klu_factor_loop:
    for (int block = 0; block < Symbolic->nblocks; block++)
    {

        /* ------------------------------------------------------------------ */
        /* the block is from rows/columns k1 to k2-1 */
        /* ------------------------------------------------------------------ */

        int k1 = Symbolic->R[block];
        int k2 = Symbolic->R[block + 1];
        int nk = k2 - k1;
        PRINTF(("FACTOR BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2 - 1, nk));

        if (nk == 1)
        {

            /* -------------------------------------------------------------- */
            /* singleton case */
            /* -------------------------------------------------------------- */

            poff = Numeric->Offp[k1];
            oldcol = Symbolic->Q[k1];
            pend = Ap[oldcol + 1];
            double s = 0;

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

            int lusize = DUNITS(int, Lsize) + DUNITS(double, Lsize) +
                         DUNITS(int, Usize) + DUNITS(double, Usize);

            int lnz_block, unz_block;
            // W = Numeric->Iwork;
            // pinv = (int *)W;
            // W += nk;
            // Stack = (int *)W;
            // W += nk;
            // Flag = (int *)W;
            // W += nk;
            // Lpend = (int *)W;
            // W += nk;
            // Ap_pos = (int *)W;
            // W += nk;

            int pinv[MAX_SIZE], Stack[MAX_SIZE], Flag[MAX_SIZE], Ap_pos[MAX_SIZE], Lpend[MAX_SIZE], Pblock[MAX_SIZE];

            Numeric->LUsize[block] = KLU_kernel(nk, Ap, Ai, Ax, Symbolic->Q, lusize, pinv, Pblock, &Numeric->LUbx[Numeric->lusize_sum], Numeric->Udiag + k1, Numeric->Llen + k1, Numeric->Ulen + k1, Numeric->Lip + k1, Numeric->Uip + k1, &lnz_block, &unz_block, Numeric->Xwork, Stack, Flag, Ap_pos, Lpend, k1, Numeric->Pinv, Numeric->Rs, Numeric->Offp, Numeric->Offi, Numeric->Offx, Common);
            if (Common->status < KLU_OK || (Common->status == KLU_SINGULAR && Common->halt_if_singular))
            {
                /* out of memory, invalid inputs, or singular */
                return;
            }

            Numeric->lusize_sum += Numeric->LUsize[block];

            Numeric->LUsize[block] = Numeric->lusize_sum - Numeric->LUsize[block];

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
