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
#include <iostream>
#include "xcl2.hpp"
#include <algorithm>
#include <vector>
#include "klu_factor.h"
#include "klu_solve.h"
#include "mmio.h"

auto constexpr num_cu = 2;

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
    std::vector<cl::Kernel> krnl_lu(num_cu);
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
        std::cout << "Error reading matrix file" << std::endl;
        return 1;
    }

    if (read_bmatrix(bmatrix, b, &nrhs))
    {
        std::cout << "Error reading bmatrix file" << std::endl;
        return 1;
    }

    std::cout << "nrhs: " << nrhs << ", bnumber: " << b.size() << std::endl;

    std::vector<double> b2(b.begin(), b.end());

    klu_symbolic Symbolic;
    klu_numeric Numeric;
    klu_common Common;
    klu_defaults(&Common);
    Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);

    Numeric.Xwork = (double *)malloc(n * nrhs * sizeof(double));

    klu_factor(Ap.data(), Ai.data(), Ax.data(), &Symbolic, &Numeric, &Common);
    klu_solve(&Symbolic, &Numeric, n, nrhs, b2.data(), &Common);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs - 1; j++)
            std::cout << "x [" << i << "," << j << "] = " << b2[i + n * j] << " ";

        std::cout << "x [" << i << "," << nrhs - 1 << "] = " << b2[i + n * (nrhs - 1)] << std::endl;
        if (i > 10)
            break;
    }

    // I/O data vectors
    std::vector<int, aligned_allocator<int>> R(Symbolic.R, Symbolic.R + n + 1),
        Q(Symbolic.Q, Symbolic.Q + n),
        Pnum(Numeric.Pnum, Numeric.Pnum + n),
        Lip(Numeric.Lip, Numeric.Lip + n),
        Llen(Numeric.Llen, Numeric.Llen + n),
        LUsize(Numeric.LUsize, Numeric.LUsize + n),
        Uip(Numeric.Uip, Numeric.Uip + n),
        Ulen(Numeric.Ulen, Numeric.Ulen + n),
        Offp(Numeric.Offp, Numeric.Offp + n + 1),
        Offi(Numeric.Offi, Numeric.Offi + Symbolic.nzoff);
    std::vector<double, aligned_allocator<double>> Rs(Numeric.Rs, Numeric.Rs + n),
        LUbx(Numeric.LUbx, Numeric.LUbx + Numeric.lusize_sum),
        Udiag(Numeric.Udiag, Numeric.Udiag + n),
        Offx(Numeric.Offx, Numeric.Offx + Symbolic.nzoff);

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
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err));

        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS)
        {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        }
        else
        {
            std::cout << "Device[" << i << "]: program successful!\n";
            // Creating Kernel objects
            for (int i = 0; i < num_cu; i++)
            {
                OCL_CHECK(err, krnl_lu[i] = cl::Kernel(program, "lu", &err));
            }
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device)
    {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    std::vector<cl::Buffer> buffer_in0(num_cu), buffer_in1(num_cu), buffer_in2(num_cu), buffer_in3(num_cu), buffer_in4(num_cu), buffer_in5(num_cu), buffer_in6(num_cu), buffer_in7(num_cu), buffer_in8(num_cu), buffer_in9(num_cu), buffer_in10(num_cu), buffer_in11(num_cu), buffer_in12(num_cu), buffer_in13(num_cu);
    std::vector<cl::Buffer> buffer_inout(num_cu);

    // Allocate Buffer in Global Memory
    // Buffers are allocated using CL_MEM_USE_HOST_PTR for efficient memory and
    // Device-to-host communication

    int chunk_size = 1;
    if (nrhs > 1)
        chunk_size = nrhs / num_cu;

    std::cout << "chunk_size: " << chunk_size << std::endl;

    for (int i = 0; i < num_cu; i++)
    {
        OCL_CHECK(err, buffer_in0[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * R.size(), R.data(), &err));
        OCL_CHECK(err, buffer_in1[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Q.size(), Q.data(), &err));
        OCL_CHECK(err, buffer_in2[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Pnum.size(), Pnum.data(), &err));
        OCL_CHECK(err, buffer_in3[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * Rs.size(), Rs.data(), &err));
        OCL_CHECK(err, buffer_in4[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Lip.size(), Lip.data(), &err));
        OCL_CHECK(err, buffer_in5[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Llen.size(), Llen.data(), &err));
        OCL_CHECK(err, buffer_in6[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * LUbx.size(), LUbx.data(), &err));
        OCL_CHECK(err, buffer_in7[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * LUsize.size(), LUsize.data(), &err));
        OCL_CHECK(err, buffer_in8[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Uip.size(), Uip.data(), &err));
        OCL_CHECK(err, buffer_in9[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Ulen.size(), Ulen.data(), &err));
        OCL_CHECK(err, buffer_in10[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * Udiag.size(), Udiag.data(), &err));
        OCL_CHECK(err, buffer_in11[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Offp.size(), Offp.data(), &err));
        OCL_CHECK(err, buffer_in12[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * Offi.size(), Offi.data(), &err));
        OCL_CHECK(err, buffer_in13[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * Offx.size(), Offx.data(), &err));

        if (nrhs > 1)
        {
            if (i == num_cu - 1)
            {
                OCL_CHECK(err, buffer_inout[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * n * (nrhs - i * chunk_size), b.data() + i * n * chunk_size, &err));
                OCL_CHECK(err, err = krnl_lu[i].setArg(17, nrhs - i * chunk_size));
            }
            else
            {
                OCL_CHECK(err, buffer_inout[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * n * chunk_size, b.data() + i * n * chunk_size, &err));
                OCL_CHECK(err, err = krnl_lu[i].setArg(17, chunk_size));
            }
        }
        else
        {
            if (i)
            {
                OCL_CHECK(err, buffer_inout[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * b.size(), b.data(), &err));
                OCL_CHECK(err, err = krnl_lu[i].setArg(17, 0));
            }
            else
            {
                OCL_CHECK(err, buffer_inout[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * b.size(), b.data(), &err));
                OCL_CHECK(err, err = krnl_lu[i].setArg(17, 1));
            }
        }
    }

    for (int i = 0; i < num_cu; i++)
    {
        OCL_CHECK(err, err = krnl_lu[i].setArg(0, buffer_in0[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(1, buffer_in1[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(2, buffer_in2[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(3, buffer_in3[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(4, buffer_in4[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(5, buffer_in5[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(6, buffer_in6[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(7, buffer_in7[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(8, buffer_in8[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(9, buffer_in9[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(10, buffer_in10[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(11, buffer_in11[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(12, buffer_in12[i]));
        OCL_CHECK(err, err = krnl_lu[i].setArg(13, buffer_in13[i]));

        OCL_CHECK(err, err = krnl_lu[i].setArg(14, n));
        OCL_CHECK(err, err = krnl_lu[i].setArg(15, Numeric.lusize_sum));
        OCL_CHECK(err, err = krnl_lu[i].setArg(16, Symbolic.nblocks));

        OCL_CHECK(err, err = krnl_lu[i].setArg(18, buffer_inout[i]));

        // Copy input data to device global memory
        OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in0[i], buffer_in1[i], buffer_in2[i], buffer_in3[i], buffer_in4[i], buffer_in5[i], buffer_in6[i], buffer_in7[i], buffer_in8[i], buffer_in9[i], buffer_in10[i], buffer_in11[i], buffer_in12[i], buffer_in13[i], buffer_inout[i]}, 0)); /* 0 means from host*/
    }
    OCL_CHECK(err, err = q.finish());

    for (int i = 0; i < num_cu; i++)
    {
        // Launch the kernel
        OCL_CHECK(err, err = q.enqueueTask(krnl_lu[i]));
    }
    OCL_CHECK(err, err = q.finish());

    // Copy Result from Device Global Memory to Host Local Memory
    for (int i = 0; i < num_cu; i++)
    {
        OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_inout[i]}, CL_MIGRATE_MEM_OBJECT_HOST));
    }
    OCL_CHECK(err, err = q.finish());
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
            }
            else
            {
                if (i > 10)
                    continue;
                else if (j < nrhs - 1)
                    std::cout << "x[" << i << "," << j << "] = " << b[i + n * j] << " ";
                else
                    std::cout << "x[" << i << "," << j << "] = " << b[i + n * j] << std::endl;
            }
        }
    }

    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << std::endl;
    return (match ? EXIT_SUCCESS : EXIT_FAILURE);
}
