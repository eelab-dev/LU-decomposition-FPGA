# Acceleration of LU Decomposition on FPGAs

## Background
Simulation Program with Integrated Circuit Emphasis or SPICE has now been widely used in the IC design and verification. Solving of sparse matrices often takes up most of the SPICE simulation time. Lower–upper (LU) decomposition is the most commonly used method to solve matrices. It factorizes a matrix into two factors – a lower triangular matrix L and an upper triangular matrix U. In this way, we only need to solve triangular systems to get results. However, the sparse-matrix computation is hard to parallelize on regular processors due to the irregular structure of the matrices. Modern FPGAs, however, have the potential to compute these hard-to-parallelise problems more efficiently due to its flexible reconfigurability.

## Implement Dense LU Decomposition in C++
### Code
  [here](C/vector2D2.cpp)

### Result

## Implement LU Decomposition with Partial Pivoting in C++
### Code
  [here](C/lupivot_cmd.cpp)

## Implementation in Vitis HLS


## Sparse LU Decomposition