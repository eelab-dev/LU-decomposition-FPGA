Team number: xohw22-006

Project name: Acceleration of LU Decomposition on FPGAs

Link to YouTube Video(s): https://youtu.be/Pth-gqt7MiY

Link to project repository: https://github.com/danchitnis/LU-decomposition-FPGA

University name: the University of Edinburgh

Participant(s): Yichen Zhang

Email: s2130520@ed.ac.uk

Supervisor name: Dr Danial Chitnis

Supervisor e-mail: d.chitnis@ed.ac.uk

Board used: Alveo U280 Data Center Accelerator Card

Software Version: v2021.2

Brief description of project: This project aims to accelerate the solving of linear systems by parallelizing the matrix decomposition on FPGAs. Solving large linear systems is part of everyday engineering design, including the integrated circuit simulation, where the majority of the matrices are sparse. Lower–upper (LU) decomposition is the most commonly used method to solve the sparse linear system. However, the sparse-matrix decomposition is hard to parallelize on regular processors due to the irregular structure of the input matrices. Modern FPGAs have the potential to compute these hard-to-parallelize problems more efficiently due to their reconfigurable structure, such as flexible memory access, loop flattening and unrolling. As a result, these FPGA implementations may lead to higher compute throughput and efficiency.

Description of archive (explain directory structure, documents and source files):
- Matrix_Sample/
  Some matrices used for test.
- myKLU/
  Software version for KLU decomposition
- Vitis/myKLU/host
  Host code for FPGA implementation
- Vitis/myKLU/host
  Kernel code for FPGA implementation

Instructions to build and test project

For CPU version:
Step 1: Go to directory ./myKLU
Step 2: Check if clang is installed. If not, modify `Makefile` to use appropriate compiler. Typically `GCC` is acceptable.
Step 3: To make static library only, run `make library`
Step 4: To run software version of LU decomposition, run `make klu_kernel`. The default matrix file should be put in the same directory as the executable program.
Step 5: To benchmark it with matrices in directory ../Matrix_Sample, run `make klu_bench`

For FPGA version:
Step 1: Import the project to Vitis 2021.2.
Step 2: Compile the CPU version above first, FPGA version requires the results from CPU version to verify that the results are correct.
Step 3: Check the host C/C++ build settings. Ensure the include and lib settings are correct.
Step 4: Build the project.