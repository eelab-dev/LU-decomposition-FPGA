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
- Video link: [https://youtu.be/Pth-gqt7MiY](https://youtu.be/Pth-gqt7MiY)
- Report: [report_xohw22-006.pdf](../Documents/xohw/out/report_xohw22-006.pdf)

## File Organisation
- Matrix_Sample/
  Some matrices used for test.
- myKLU/
  Software version for KLU decomposition
- Vitis/myKLU/host
  Host code for FPGA implementation
- Vitis/myKLU/host
  Kernel code for FPGA implementation

## Experimental Results
Matrices used for test
|      Matrix     | Order |  NNZ | Sparsity | Pattern Symmetry | Numeric Symmetry |
|:---------------:|:-----:|:----:|:--------:|:----------------:|:----------------:|
|     rajat11     |  135  |  665 |  3.65\%  |      89.10\%     |       63\%       |
|     rajat14     |  180  | 1475 |  4.55\%  |       100\%      |      2.50\%      |
|     rajat05     |  301  | 1250 |  1.38\%  |       77\%       |      70.60\%     |
| oscil\_dcop\_01 |  430  | 1544 |  0.84\%  |      97.60\%     |      69.80\%     |
|  fpga\_dcop\_01 |  1220 | 5892 |  0.40\%  |      81.80\%     |      27.30\%     |

![cpu](Images/klu_cpu.svg)
CPU

![fpga](Images/klu_fpga.svg)
FPGA

For smaller matrices, FPGA tends to take longer time to solve per right hand side vectors than CPU. However, when the matrices becomes larger, FPGA tends to be faster than CPU, with a speedup of about 1.2.