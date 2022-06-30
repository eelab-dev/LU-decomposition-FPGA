# Software implementation of LU Decomposition
---

- AMD: a set of routines for permuting sparse matrices prior to factorization.
  Requires SuiteSparse_config, in the ./SuiteSparse_config directory relative to this directory.

- BTF: a software package for permuting a matrix into block upper triangular form.
  Requires SuiteSparse_config, in the ./SuiteSparse_config directory relative to this directory.

- KLU: a set of routines for solving sparse linear systems of equations. It is particularly well-suited to matrices arising in SPICE-like circuit simulation applications.
  Here, only the symbolic analysis routines are implemented in this directory.

- myKLU: Modified KLU. Numeric factorisation and solving routines are implemented here. With SIMD for parallelising multiple right-hand side solving.

- SuiteSparse_config: This directory contains a default SuiteSparse_config.mk file.  It tries to detect your system (Linux, SunOS, or Mac), which compiler to use (icc or cc), which BLAS and LAPACK library to use (Intel MKL is strongly preferred), and whether or not to compile with CUDA.

## Building Instructions
1. Go to directory ./myKLU
2. Check if clang is installed. If not, modify `Makefile` to use appropriate compiler. Typically `GCC` is acceptable.
3. To make static library only, run `make library`
4. To run software version of LU decomposition, run `make klu_kernel`. The default matrix file should be put in the same directory as the executable program.
5. To benchmark it with matrices in directory ../Matrix_Sample, run `make klu_bench`