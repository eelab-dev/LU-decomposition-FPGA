#include <iostream>
#include <assert.h>
#include <vector>

#define Int int
#define BTF_FLIP(j) (-(j)-2)
#define BTF_ISFLIPPED(j) ((j) < -1)
#define BTF_UNFLIP(j) ((BTF_ISFLIPPED(j)) ? BTF_FLIP(j) : (j))
#define UNVISITED (-2)  /* Flag [j] = UNVISITED if node j not visited yet */
#define UNASSIGNED (-1) /* Flag [j] = UNASSIGNED if node j has been visited, \
                         * but not yet assigned to a strongly-connected      \
                         * component (aka block).  Flag [j] = k (k in the    \
                         * range 0 to nblocks-1) if node j has been visited  \
                         * (and completed, with its postwork done) and       \
                         * assigned to component k. */
#define ASSERT(a) assert(a)
#define TRUE 1
#define FALSE 0
#define EMPTY (-1)
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

static void dfs(
    /* inputs, not modified on output: */
    Int j,    /* start the DFS at node j */
    Int Ap[], /* size n+1, column pointers for the matrix A */
    Int Ai[], /* row indices, size nz = Ap [n] */
    Int Q[],  /* input column permutation */

    /* inputs, modified on output (each array is of size n): */
    Int Time[],       /* Time [j] = "time" that node j was first visited */
    Int Flag[],       /* Flag [j]: see above */
    Int Low[],        /* Low [j]: see definition below */
    Int *p_nblocks,   /* number of blocks (aka strongly-connected-comp.)*/
    Int *p_timestamp, /* current "time" */

    /* workspace, not defined on input or output: */
    Int Cstack[], /* size n, output stack to hold nodes of components */
    Int Jstack[], /* size n, stack for the variable j */
    Int Pstack[]  /* size n, stack for the variable p */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables, and initializations */
    /* ---------------------------------------------------------------------- */

    /* local variables, but "global" to all DFS levels: */
    Int chead; /* top of Cstack */
    Int jhead; /* top of Jstack and Pstack */

    /* variables that are purely local to any one DFS level: */
    Int i;      /* edge (j,i) considered; i can be next node to traverse */
    Int parent; /* parent of node j in the DFS tree */
    Int pend;   /* one past the end of the adjacency list for node j */
    Int jj;     /* column j of A*Q is column jj of the input matrix A */

    /* variables that need to be pushed then popped from the stack: */
    Int p; /* current index into the adj. list for node j */
    /* the variables j and p are stacked in Jstack and Pstack */

    /* local copies of variables in the calling routine */
    Int nblocks = *p_nblocks;
    Int timestamp = *p_timestamp;

    /* ---------------------------------------------------------------------- */
    /* start a DFS at node j (same as the recursive call dfs (EMPTY, j)) */
    /* ---------------------------------------------------------------------- */

    chead = 0;     /* component stack is empty */
    jhead = 0;     /* Jstack and Pstack are empty */
    Jstack[0] = j; /* put the first node j on the Jstack */
    ASSERT(Flag[j] == UNVISITED);

    while (jhead >= 0)
    {
        j = Jstack[jhead]; /* grab the node j from the top of Jstack */

        /* determine which column jj of the A is column j of A*Q */
        jj = (Q == (Int *)NULL) ? (j) : (BTF_UNFLIP(Q[j]));
        pend = Ap[jj + 1]; /* j's row index list ends at Ai [pend-1] */

        if (Flag[j] == UNVISITED)
        {

            /* -------------------------------------------------------------- */
            /* prework at node j */
            /* -------------------------------------------------------------- */

            /* node j is being visited for the first time */
            Cstack[++chead] = j; /* push j onto the stack */
            timestamp++;         /* get a timestamp */
            Time[j] = timestamp; /* give the timestamp to node j */
            Low[j] = timestamp;
            Flag[j] = UNASSIGNED; /* flag node j as visited */

            /* -------------------------------------------------------------- */
            /* set Pstack [jhead] to the first entry in column j to scan */
            /* -------------------------------------------------------------- */

            Pstack[jhead] = Ap[jj];
        }

        /* ------------------------------------------------------------------ */
        /* DFS rooted at node j (start it, or continue where left off) */
        /* ------------------------------------------------------------------ */

        for (p = Pstack[jhead]; p < pend; p++)
        {
            i = Ai[p]; /* examine the edge from node j to node i */
            if (Flag[i] == UNVISITED)
            {
                /* Node i has not been visited - start a DFS at node i.
                 * Keep track of where we left off in the scan of adjacency list
                 * of node j so we can restart j where we left off. */
                Pstack[jhead] = p + 1;
                /* Push i onto the stack and immediately break
                 * so we can recurse on node i. */
                Jstack[++jhead] = i;
                ASSERT(Time[i] == EMPTY);
                ASSERT(Low[i] == EMPTY);
                /* break here to do what the recursive call dfs (j,i) does */
                break;
            }
            else if (Flag[i] == UNASSIGNED)
            {
                /* Node i has been visited, but still unassigned to a block
                 * this is a back or cross edge if Time [i] < Time [j].
                 * Note that i might equal j, in which case this code does
                 * nothing. */
                ASSERT(Time[i] > 0);
                ASSERT(Low[i] > 0);
                Low[j] = MIN(Low[j], Time[i]);
            }
        }

        if (p == pend)
        {
            /* If all adjacent nodes of j are already visited, pop j from
             * Jstack and do the post work for node j.  This also pops p
             * from the Pstack. */
            jhead--;

            /* -------------------------------------------------------------- */
            /* postwork at node j */
            /* -------------------------------------------------------------- */

            /* determine if node j is the head of a component */
            if (Low[j] == Time[j])
            {
                /* pop all nodes in this SCC from Cstack */
                while (TRUE)
                {
                    ASSERT(chead >= 0);  /* stack not empty (j in it) */
                    i = Cstack[chead--]; /* pop a node from the Cstack */
                    ASSERT(i >= 0);
                    ASSERT(Flag[i] == UNASSIGNED);
                    Flag[i] = nblocks; /* assign i to current block */
                    if (i == j)
                        break; /* current block ends at j */
                }
                nblocks++; /* one more block has been found */
            }
            /* update Low [parent], if the parent exists */
            if (jhead >= 0)
            {
                parent = Jstack[jhead];
                Low[parent] = MIN(Low[parent], Low[j]);
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* cleanup: update timestamp and nblocks */
    /* ---------------------------------------------------------------------- */

    *p_timestamp = timestamp;
    *p_nblocks = nblocks;
}

Int augment(
    Int k,         /* which stage of the main loop we're in */
    Int Ap[],      /* column pointers, size n+1 */
    Int Ai[],      /* row indices, size nz = Ap [n] */
    Int Match[],   /* size n,  Match [i] = j if col j matched to i */
    Int Cheap[],   /* rows Ai [Ap [j] .. Cheap [j]-1] alread matched */
    Int Flag[],    /* Flag [j] = k if j already visited this stage */
    Int Istack[],  /* size n.  Row index stack. */
    Int Jstack[],  /* size n.  Column index stack. */
    Int Pstack[],  /* size n.  Keeps track of position in adjacency list */
    double *work,  /* work performed by the depth-first-search */
    double maxwork /* maximum work allowed */
)
{
    /* local variables, but "global" to all DFS levels: */
    Int found; /* true if match found.  */
    Int head;  /* top of stack */

    /* variables that are purely local to any one DFS level: */
    Int j2;   /* the next DFS goes to node j2 */
    Int pend; /* one past the end of the adjacency list for node j */
    Int pstart;
    Int quick;

    /* variables that need to be pushed then popped from the stack: */
    Int i; /* the row tentatively matched to i if DFS successful */
    Int j; /* the DFS is at the current node j */
    Int p; /* current index into the adj. list for node j */
    /* the variables i, j, and p are stacked in Istack, Jstack, and Pstack */

    quick = (maxwork > 0);

    /* start a DFS to find a match for column k */
    found = FALSE;
    i = EMPTY;
    head = 0;
    Jstack[0] = k;
    ASSERT(Flag[k] != k);

    while (head >= 0)
    {
        j = Jstack[head];
        pend = Ap[j + 1];

        if (Flag[j] != k) /* a node is not yet visited */
        {

            /* -------------------------------------------------------------- */
            /* prework for node j */
            /* -------------------------------------------------------------- */

            /* first time that j has been visited */
            Flag[j] = k;
            /* cheap assignment: find the next unmatched row in col j.  This
             * loop takes at most O(nnz(A)) time for the sum total of all
             * calls to augment. */
            for (p = Cheap[j]; p < pend && !found; p++)
            {
                i = Ai[p];
                found = (Match[i] == EMPTY);
            }
            Cheap[j] = p;

            /* -------------------------------------------------------------- */

            /* prepare for DFS */
            if (found)
            {
                /* end of augmenting path, column j matched with row i */
                Istack[head] = i;
                break;
            }
            /* set Pstack [head] to the first entry in column j to scan */
            Pstack[head] = Ap[j];
        }

        /* ------------------------------------------------------------------ */
        /* quick return if too much work done */
        /* ------------------------------------------------------------------ */

        if (quick && *work > maxwork)
        {
            /* too much work has been performed; abort the search */
            return (EMPTY);
        }

        /* ------------------------------------------------------------------ */
        /* DFS for nodes adjacent to j */
        /* ------------------------------------------------------------------ */

        /* If cheap assignment not made, continue the depth-first search.  All
         * rows in column j are already matched.  Add the adjacent nodes to the
         * stack by iterating through until finding another non-visited node.
         *
         * It is the following loop that can force maxtrans to take
         * O(n*nnz(A)) time. */

        pstart = Pstack[head];
        for (p = pstart; p < pend; p++)
        {
            i = Ai[p];
            j2 = Match[i];
            ASSERT(j2 != EMPTY);
            if (Flag[j2] != k)
            {
                /* Node j2 is not yet visited, start a depth-first search on
                 * node j2.  Keep track of where we left off in the scan of adj
                 * list of node j so we can restart j where we left off. */
                Pstack[head] = p + 1;
                /* Push j2 onto the stack and immediately break so we can
                 * recurse on node j2.  Also keep track of row i which (if this
                 * search for an augmenting path works) will be matched with the
                 * current node j. */
                Istack[head] = i;
                Jstack[++head] = j2;
                break;
            }
        }

        /* ------------------------------------------------------------------ */
        /* determine how much work was just performed */
        /* ------------------------------------------------------------------ */

        *work += (p - pstart + 1);

        /* ------------------------------------------------------------------ */
        /* node j is done, but the postwork is postponed - see below */
        /* ------------------------------------------------------------------ */

        if (p == pend)
        {
            /* If all adjacent nodes of j are already visited, pop j from
             * stack and continue.  We failed to find a match. */
            head--;
        }
    }

    /* postwork for all nodes j in the stack */
    /* unwind the path and make the corresponding matches */
    if (found)
    {
        for (p = head; p >= 0; p--)
        {
            j = Jstack[p];
            i = Istack[p];

            /* -------------------------------------------------------------- */
            /* postwork for node j */
            /* -------------------------------------------------------------- */
            /* if found, match row i with column j */
            Match[i] = j;
        }
    }
    return (found);
}

int main(void)
{
    int n = 10;
    int Ap[] = {0, 2, 4, 7, 9, 13, 14, 18, 21, 23, 24};
    int Ai[] = {0, 7, 1, 5, 0, 2, 7, 3, 8, 3, 4, 5, 8, 5, 0, 6, 7, 8, 4, 5, 7, 2, 8, 9};
    double Ax[] = {8, 8, 2, 4, 10, 3, 10, 1, 5, 2, 1, 4, 10, 2, 3, 5, 3, 15, 4, 16, 3, 7, 2, 9};

    int j;
    int Q[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<int> Flag(n), Low(n), Time(n), Cstack(n), Jstack(n), Pstack(n);
    int timestamp = 0; /* each node given a timestamp when it is visited */
    int nblocks = 0;   /* number of blocks found so far */

    for (int j = 0; j < n; j++)
    {
        Flag[j] = UNVISITED;
        Low[j] = EMPTY;
        Time[j] = EMPTY;
    }

    for (j = 0; j < n; j++)
    {
        /* node j is unvisited or assigned to a block. Cstack is empty. */
        ASSERT(Flag[j] == UNVISITED || (Flag[j] >= 0 && Flag[j] < nblocks));
        if (Flag[j] == UNVISITED)
        {
            /* non-recursive dfs (default) */
            dfs(j, Ap, Ai, Q, Time.data(), Flag.data(), Low.data(), &nblocks, &timestamp, Cstack.data(), Jstack.data(), Pstack.data());

            std::cout << j << ":" << std::endl;
            for (int i = 0; i < n; i++)
            {
                std::cout << "Time[" << i << "]= " << Time[i] << "\tFlag[" << i << "]= " << Flag[i] << "\tLow[" << i << "]= " << Low[i] << "\tCstack[" << i << "]= " << Cstack[i] << "\tJstack[" << i << "]= " << Jstack[i] << "\tPstack[" << i << "]= " << Pstack[i] << std::endl;
            }
        }
    }
    // assert(1 == 2);

    return 0;
}