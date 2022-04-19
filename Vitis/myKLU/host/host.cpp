/**
 * Copyright (C) 2019-2021 Xilinx, Inc
 *
 * Licensed under the Apache License, Version 2.0 (the "License"). You may
 * not use this file except in compliance with the License. A copy of the
 * License is located at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

#include "xcl2.hpp"
#include <algorithm>
#include <vector>
#include "klu.h"
#include "mmio.h"

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string binaryFile = argv[1];

    cl_int err;
    cl::Context context;
    cl::Kernel krnl_vector_add;
    cl::CommandQueue q;
    // Allocate Memory in Host Memory
    // When creating a buffer with user pointer (CL_MEM_USE_HOST_PTR), under the
    // hood user ptr
    // is used if it is properly aligned. when not aligned, runtime had no choice
    // but to create
    // its own host side buffer. So it is recommended to use this allocator if
    // user wish to
    // create buffer using CL_MEM_USE_HOST_PTR to align user buffer to page
    // boundary. It will
    // ensure that user buffer is used when user create Buffer/Mem object with
    // CL_MEM_USE_HOST_PTR
    std::vector<int, aligned_allocator<int>> Ap, Ai;
    std::vector<double, aligned_allocator<double>> Ax, b;

    std::string homeDir = getenv("HOME");

    std::string filename, bmatrix;
    std::cout << "Left matrix file path (default - " << homeDir << "/beng-project/Matrix_Sample/host.mtx): ";
    std::getline(std::cin, filename);
    if (filename.empty())
        filename = homeDir + "/beng-project/Matrix_Sample/host.mtx";

    std::cout << "B matrix file path (default - " << homeDir << "/beng-project/Matrix_Sample/host_b.mtx): ";
    std::getline(std::cin, bmatrix);
    if (bmatrix.empty())
        bmatrix = homeDir + "/beng-project/Matrix_Sample/host_b.mtx";

    int n, nrhs;

    if (read_sparse(filename, &n, Ap, Ai, Ax))
    {
        printf("Error\n");
        return 1;
    }

    if (read_bmatrix(bmatrix, b, &nrhs))
    {
        printf("Error bmatrix\n");
        return 1;
    }

    printf("nrhs: %d, bnumber:%d\n", nrhs, b.size());

    std::vector<double> b2(b.begin(), b.end());

    klu_symbolic Symbolic;
    klu_numeric *Numeric;
    klu_common Common;

    klu_defaults(&Common);
    Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);
    Numeric = klu_factor(Ap.data(), Ai.data(), Ax.data(), &Symbolic, &Common);
    klu_solve(&Symbolic, Numeric, n, nrhs, b2.data(), &Common);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs - 1; j++)
            printf("x [%d,%d] = %g\t", i, j, b2[i + n * j]);
        printf("x [%d,%d] = %g\n", i, nrhs - 1, b2[i + n * (nrhs - 1)]);
    }

    std::vector<int, aligned_allocator<int>> P(Symbolic.P, Symbolic.P + n), Q(Symbolic.Q, Symbolic.Q + n), R(Symbolic.R, Symbolic.R + n + 1);
    std::vector<double, aligned_allocator<double>> Lnz(Symbolic.Lnz, Symbolic.Lnz + n);

    // OPENCL HOST CODE AREA START
    // get_xil_devices() is a utility API which will find the xilinx
    // platforms and will return list of devices connected to Xilinx platform
    auto devices = xcl::get_xil_devices();
    // read_binary_file() is a utility API which will load the binaryFile
    // and will return the pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++)
    {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));
        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS)
        {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        }
        else
        {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, krnl_vector_add = cl::Kernel(program, "lu", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device)
    {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    int p_size_bytes = sizeof(int) * Ap.size();
    int i_size_bytes = sizeof(int) * Ai.size();
    int x_size_bytes = sizeof(double) * Ax.size();
    int pq_size_bytes = sizeof(int) * n;
    int r_size_bytes = sizeof(int) * (n + 1);
    int lnz_size_bytes = sizeof(double) * n;
    int b_size_bytes = sizeof(double) * b.size();

    // Allocate Buffer in Global Memory
    // Buffers are allocated using CL_MEM_USE_HOST_PTR for efficient memory and
    // Device-to-host communication
    OCL_CHECK(err, cl::Buffer buffer_in0(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, p_size_bytes, Ap.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in1(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, i_size_bytes, Ai.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in2(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, x_size_bytes, Ax.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in3(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, pq_size_bytes, P.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in4(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, pq_size_bytes, Q.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in5(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, r_size_bytes, R.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in6(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, lnz_size_bytes, Lnz.data(), &err));

    OCL_CHECK(err, cl::Buffer buffer_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, b_size_bytes, b.data(), &err));

    OCL_CHECK(err, err = krnl_vector_add.setArg(0, buffer_in0));
    OCL_CHECK(err, err = krnl_vector_add.setArg(1, buffer_in1));
    OCL_CHECK(err, err = krnl_vector_add.setArg(2, buffer_in2));
    OCL_CHECK(err, err = krnl_vector_add.setArg(3, buffer_in3));
    OCL_CHECK(err, err = krnl_vector_add.setArg(4, buffer_in4));
    OCL_CHECK(err, err = krnl_vector_add.setArg(5, buffer_in5));
    OCL_CHECK(err, err = krnl_vector_add.setArg(6, buffer_in6));

    OCL_CHECK(err, err = krnl_vector_add.setArg(7, n));
    OCL_CHECK(err, err = krnl_vector_add.setArg(8, Symbolic.nblocks));
    OCL_CHECK(err, err = krnl_vector_add.setArg(9, Symbolic.maxblock));
    OCL_CHECK(err, err = krnl_vector_add.setArg(10, Symbolic.nzoff));
    OCL_CHECK(err, err = krnl_vector_add.setArg(11, Symbolic.nz));
    OCL_CHECK(err, err = krnl_vector_add.setArg(12, nrhs));

    OCL_CHECK(err, err = krnl_vector_add.setArg(13, buffer_output));

    // Copy input data to device global memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in0, buffer_in1, buffer_in2, buffer_in3, buffer_in4, buffer_in5, buffer_in6, buffer_output}, 0 /* 0 means from host*/));

    // Launch the Kernel
    // For HLS kernels global and local size is always (1,1,1). So, it is
    // recommended
    // to always use enqueueTask() for invoking HLS kernel
    OCL_CHECK(err, err = q.enqueueTask(krnl_vector_add));

    // Copy Result from Device Global Memory to Host Local Memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_output}, CL_MIGRATE_MEM_OBJECT_HOST));
    q.finish();
    // OPENCL HOST CODE AREA END

    // Compare the results of the Device to the simulation
    bool match = true;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs; j++)
        {
            if (std::abs(b2[i + n * j] - b[i + n * j]) > 1e-5)
            {
                std::cout << "Mismatched result x[" << i << "][" << j << "]: CPU x[" << i << "][" << j << "]=" << b2[i + n * j] << ", FPGA x[" << i << "][" << j << "]=" << b[i + n * j] << std::endl;
                match = false;
                //            break;
            }
            else
            {
                if (j < nrhs - 1)
                    printf("x[%d,%d] = %g\t", i, j, b[i + n * j]);
                else
                    printf("x[%d,%d] = %g\n", i, nrhs - 1, b[i + n * (nrhs - 1)]);
            }
        }
    }

    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << std::endl;
    return (match ? EXIT_SUCCESS : EXIT_FAILURE);
}
