/*
 *   Matrix Market I/O example program
 *
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.  This porgram does nothing useful, but
 *   illustrates common usage of the Matrix Matrix I/O routines.
 *   (See http://math.nist.gov/MatrixMarket for details.)
 *
 *   Usage:  a.out [filename] > output
 *
 *
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.
 */

#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"

typedef struct
{
    int n;
    int *Ap;
    int *Ai;
    double *Ax;
} Mtx;

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J;
    double *val;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        exit(1);

    if (M != N)
    {
        printf("Sorry, the matrix is not square\n");
        exit(1);
    }

    /* reseve memory for matrices */

    I = (int *)malloc(nz * sizeof(int));
    J = (int *)malloc(nz * sizeof(int));
    val = (double *)malloc(nz * sizeof(double));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    char num[150];
    for (int i = 0, temp; i < nz; i++)
    {
        if (fgets(num, 150, f) == NULL)
        {
            printf("Inconsistent line number\n");
            exit(1);
        }

        temp = sscanf(num, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        if (temp == 2)
            val[i] = 0;
        else if (temp < 2)
        {
            printf("No enough value\n");
            exit(1);
        }

        I[i]--; /* adjust from 1-based to 0-based */
        J[i]--;

        if (i && (J[i] - J[i - 1] < 0))
        {
            printf("Invalid column number\n");
            exit(1);
        }
    }

    if (f != stdin)
        fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i = 0; i < nz; i++)
        fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);

    Mtx mtx;
    mtx.Ap = malloc((M + 1) * sizeof(int));
    mtx.Ai = malloc(nz * sizeof(int));
    mtx.Ax = malloc(nz * sizeof(double));

    int *temp;
    double nz2 = 0;
    temp = calloc(nz, sizeof(int));

    for (int k = 0; k < nz; k++)
        temp[J[k]]++;

    int temp2 = 0;
    for (int i = 0; i < M; i++)
    {
        mtx.Ap[i] = temp2;
        temp2 += temp[i];
        nz2 += temp[i];      /* also in double to avoid csi overflow */
        temp[i] = mtx.Ap[i]; /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    mtx.Ap[M] = temp2;

    int p;
    for (int k = 0; k < nz; k++)
    {
        mtx.Ai[p = temp[J[k]]++] = I[k]; /* A(i,j) is the pth entry in C */
        if (mtx.Ax)
            mtx.Ax[p] = val[k];
    }

    for (int i = 0; i <= M; i++)
        printf("Ap[%d]=%d\n", i, mtx.Ap[i]);
    for (int i = 0; i < nz; i++)
        printf("Ai[%d]=%d\tAx[%d]=%lf\n", i, mtx.Ai[i], i, mtx.Ax[i]);

    return 0;
}
