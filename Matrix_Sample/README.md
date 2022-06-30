# Matrix Samples
---

This directory contains some matrices to test the LU decomposition. All the circuit matrices are from the University of Florida Sparse Matrix Collection.

Timothy A. Davis and Yifan Hu. 2011. The University of Florida Sparse Matrix Collection. ACM Transactions on Mathematical Software 38, 1, Article 1 (December 2011), 25 pages. DOI: https://doi.org/10.1145/2049662.2049663

- `host.mtx`: Matrix to be solved.
- `host_b.mtx`: B matrix (right-hand side vectors) file.
- `gen_b.cpp`: Code to generate the B matrix file.
```bash
g++ gen_b.cpp -o gen_b && ./gen_b
```