#include <iostream>
#include "../myKLU/include/mmio.h"

int main(void)
{
    char filename[] = "./host.mtx";
    char bname[] = "./host_b.mtx";

    FILE *f;
    MM_typecode matcode;
    int M, N, nz;

    if ((f = fopen(filename, "r")) == NULL)
        return -1;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", filename);
        return -1;
    }

    if (!((mm_is_real(matcode) || mm_is_pattern(matcode)) && mm_is_matrix(matcode) &&
          mm_is_sparse(matcode)))
    {
        fprintf(stderr, "Sorry, this application does not support ");
        fprintf(stderr, "Market Market type: [%s]\n",
                mm_typecode_to_str(matcode));
        return -1;
    }

    /* find out size of sparse matrix: M, N, nz .... */

    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0)
    {
        fprintf(stderr, "read_unsymmetric_sparse(): could not parse matrix size.\n");
        return -1;
    }

    fclose(f);

    int nrhs;
    std::cout << "B matrix size: ";
    std::cin >> nrhs;
    std::cout << "Size: " << nrhs << std::endl;

    if ((f = fopen(bname, "w")) == NULL)
        return -1;

    mm_write_banner(f, matcode);

    fprintf(f, "%d %d\n", M, nrhs);
    for (int j = 0; j < nrhs; j++)
        for (int i = 0; i < M; i++)
            fprintf(f, "%d\n", i + j);

    fclose(f);
    return 0;
}