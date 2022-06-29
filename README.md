# Acceleration of LU Decomposition on FPGAs

Team Number: xohw22-006
***

## Background
Simulation Program with Integrated Circuit Emphasis or SPICE has now been widely used in the IC design and verification. Solving of sparse matrices often takes up most of the SPICE simulation time. Lower–upper (LU) decomposition is the most commonly used method to solve matrices. It factorizes a matrix into two factors – a lower triangular matrix L and an upper triangular matrix U. In this way, we only need to solve triangular systems to get results. However, the sparse-matrix computation is hard to parallelize on regular processors due to the irregular structure of the matrices. Modern FPGAs, however, have the potential to compute these hard-to-parallelise problems more efficiently due to its flexible reconfigurability.

## Team and Project Information
- University name: the University of Edinburgh
- Supervisor: Dr Danial Chitnis ([chitnis@ed.ac.uk]((mailto:chitnis@ed.ac.uk)))
- Student: Yichen Zhang ([s2130520@ed.ac.uk](mailto:s2130520@ed.ac.uk))
- Board: Alveo U280 Data Center Accelerator Card
- Software version: v2021.2
- Video link: [text](https://)
- Report: [report_xohw22-006.pdf](../Documents/xohw/out/report_xohw22-006.pdf)

## File Organisation
- Matrix_Sample/
  Some matrices used for test.
- LU/
  Software version for LU decomposition

## Experimental Results
