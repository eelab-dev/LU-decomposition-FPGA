#include <iostream>
#include <chrono>
#include <vector>
#include <assert.h>

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

    std::vector<int> Flag(n), Cheap(n), Istack(n), Jstack(n), Pstack(n), Match(n);
    int maxwork = 0, result, nmatch, work_limit_reached;
    double *work;
    for (int j = 0; j < n; j++)
    {
        Cheap[j] = Ap[j];
        Flag[j] = EMPTY;
    }

    /* all rows and columns are currently unmatched */
    for (int i = 0; i < n; i++)
    {
        Match[i] = EMPTY;
    }

    if (maxwork > 0)
    {
        maxwork *= Ap[n];
    }
    work = 0;

    for (int k = 0; k < n; k++)
    {
        /* find an augmenting path to match some row i to column k */
        result = augment(k, Ap, Ai, Match.data(), Cheap.data(), Flag.data(), Istack.data(), Jstack.data(), Pstack.data(), work, maxwork);
        if (result == TRUE)
        {
            /* we found it.  Match [i] = k for some row i has been done. */
            nmatch++;
        }
        else if (result == EMPTY)
        {
            /* augment gave up because of too much work, and no match found */
            work_limit_reached = TRUE;
        }

        std::cout << k << ":" << std::endl;
        for (int i = 0; i < n; i++)
        {
            std::cout << "Match[" << i << "]= " << Match[i] << "\tCheap[" << i << "]= " << Cheap[i] << "\tFlag[" << i << "]= " << Flag[i] << "\tIstack[" << i << "]= " << Istack[i] << "\tJstack[" << i << "]= " << Jstack[i] << "\tPstack[" << i << "]= " << Pstack[i] << std::endl;
        }
    }

    std::cout << nmatch << std::endl;

    return 0;
}