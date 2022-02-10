#include <stdio.h>
#include <malloc.h>
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
    int *Lip;       /* size n. pointers into LUbx[block] for L */
    int *Uip;       /* size n. pointers into LUbx[block] for U */
    int *Llen;      /* size n. Llen [k] = # of entries in kth column of L */
    int *Ulen;      /* size n. Ulen [k] = # of entries in kth column of U */
    void **LUbx;    /* L and U indices and entries (excl. diagonal of U) */
    size_t *LUsize; /* size of each LUbx [block], in sizeof (Unit) */
    void *Udiag;    /* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    double *Rs; /* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    size_t worksize; /* size (in bytes) of Work */
    void *Work;      /* workspace */
    void *Xwork;     /* alias into Numeric->Work */
    int *Iwork;      /* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    int *Offp;  /* size n+1, column pointers */
    int *Offi;  /* size nzoff, row indices */
    void *Offx; /* size nzoff, numerical values */
    int nzoff;

} klu_numeric;

void get_pointer(int *x, int l);

void get_pointer2(int **x);

int main(void)
{
    klu_numeric *x;
    // int **lu = malloc(5 * sizeof(int *));
    // int a[100];
    // for (int i = 0; i < 100; i++)
    //     a[i] = 2 * i;
    // for (int i = 0; i < 5; i++)
    //     *(lu + i) = &a[5 * i];

    // int *xp;
    // xp = *(lu + 2) + 2;
    // printf("a[0]=%d\n", *xp);

    // printf("a[20]=%d\n", *(&a[19] + 1));
    int *lu = malloc(100 * sizeof(int));
    printf("0\n");
    // x->Llen = malloc(1 * sizeof(int));
    printf("1\n");
    int a[5];
    for (int i = 0; i < 100; i++)
        lu[i] = 2 * i;
    for (int i = 0; i < 5; i++)
        a[i] = 20 * i;
    printf("2\n");
    get_pointer(&lu[a[2]], 5);

    printf("sizeof(klu_numeric)=%lu\n", sizeof(x));

    printf("sizeof(int *)=%lu\n", sizeof(int *));

    return 0;
}

void get_pointer(int *x, int l)
{
    for (int i = 0; i < l; i++)
        printf("x[%d]=%d\n", i, x[i]);
}

void get_pointer2(int **x)
{
    printf("Entered\n");
    int *p;
    p = *x;
    for (int i = 0; i < 5; i++)
        printf("p[%d]=%d\n", i, p[i]);
}