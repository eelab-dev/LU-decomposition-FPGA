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

auto constexpr num_cu = 5;

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
    std::vector<double, aligned_allocator<double>> Ax, b_cpu;

    std::string filename, bmatrix;
	std::cout << "Left matrix file path (default - " << "./host.mtx): ";
	std::getline(std::cin, filename);
	if (filename.empty())
		filename = "./host.mtx";

	std::cout << "B matrix file path (default - " << "./host_b.mtx): ";
	std::getline(std::cin, bmatrix);
	if (bmatrix.empty())
		bmatrix = "./host_b.mtx";

    int n, nrhs;

    if (read_sparse(filename, &n, Ap, Ai, Ax))
    {
        std::cout << "Error reading matrix file" << std::endl;
        return 1;
    }

    if (read_bmatrix(bmatrix, b_cpu, &nrhs))
    {
        std::cout << "Error reading bmatrix file" << std::endl;
        return 1;
    }

    std::vector<double, aligned_allocator<double>> b2(b_cpu.begin(), b_cpu.end());
    std::cout << "nrhs: " << nrhs << ", bnumber: " << b_cpu.size() << std::endl;

    klu_symbolic Symbolic;
    klu_numeric Numeric;
    klu_common Common;
    klu_defaults(&Common);
    Symbolic = *klu_analyze(n, Ap.data(), Ai.data(), &Common);

    Numeric.Xwork = (double *)malloc(n * nrhs * sizeof(double));

    klu_factor(Ap.data(), Ai.data(), Ax.data(), &Symbolic, &Numeric, &Common);
    klu_solve(&Symbolic, &Numeric, n, nrhs, b_cpu.data(), &Common);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs - 1; j++)
            std::cout << "x [" << i << "," << j << "] = " << b_cpu[i + n * j] << " ";

        std::cout << "x [" << i << "," << nrhs - 1 << "] = " << b_cpu[i + n * (nrhs - 1)] << std::endl;
        if (i > 10)
            break;
    }

    // I/O data vectors
    std::vector<int, aligned_allocator<int>> int_lu;
    std::vector<double, aligned_allocator<double>> double_lu;

    int_lu.reserve(n * 9 + Symbolic.nzoff + 2);
    int_lu.insert(int_lu.end(), Symbolic.Q, Symbolic.Q + n);
    int_lu.insert(int_lu.end(), Numeric.Pnum, Numeric.Pnum + n);
    int_lu.insert(int_lu.end(), Numeric.Lip, Numeric.Lip + n);
    int_lu.insert(int_lu.end(), Numeric.Llen, Numeric.Llen + n);
    int_lu.insert(int_lu.end(), Numeric.LUsize, Numeric.LUsize + n);
    int_lu.insert(int_lu.end(), Numeric.Uip, Numeric.Uip + n);
    int_lu.insert(int_lu.end(), Numeric.Ulen, Numeric.Ulen + n);
    int_lu.insert(int_lu.end(), Symbolic.R, Symbolic.R + n + 1);
    int_lu.insert(int_lu.end(), Numeric.Offp, Numeric.Offp + n + 1);
    int_lu.insert(int_lu.end(), Numeric.Offi, Numeric.Offi + Symbolic.nzoff);

    double_lu.reserve(n * 2 + Numeric.lusize_sum + Symbolic.nzoff);
    double_lu.insert(double_lu.end(), Numeric.Rs, Numeric.Rs + n);
    double_lu.insert(double_lu.end(), Numeric.Udiag, Numeric.Udiag + n);
    double_lu.insert(double_lu.end(), Numeric.LUbx, Numeric.LUbx + Numeric.lusize_sum);
    double_lu.insert(double_lu.end(), Numeric.Offx, Numeric.Offx + Symbolic.nzoff);

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

    OCL_CHECK(err, cl::Buffer buffer_in0(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int) * int_lu.size(), int_lu.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in1(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * double_lu.size(), double_lu.data(), &err));

    std::vector<cl::Buffer> buffer_inout(num_cu);

    // Allocate Buffer in Global Memory
    // Buffers are allocated using CL_MEM_USE_HOST_PTR for efficient memory and
    // Device-to-host communication

    bool single_krnl;
    int chunk_size;
    std::vector<double, aligned_allocator<double>> b[num_cu];
    if (nrhs >= num_cu)
    {
        chunk_size = nrhs / num_cu;
        single_krnl = false;
        for (int i = 0; i < num_cu; i++)
        {
            if (i == num_cu - 1)
            {
                b[i].resize((nrhs - i * chunk_size) * n);
                std::copy(b2.begin() + i * chunk_size * n, b2.end(), b[i].begin());
            }
            else
            {
                b[i].resize(chunk_size * n);
                std::copy(b2.begin() + i * chunk_size * n, b2.begin() + (i + 1) * chunk_size * n, b[i].begin());
            }
        }
    }
    else
    {
        chunk_size = nrhs;
        single_krnl = true;
        b[0].resize(nrhs * n);
        std::copy(b2.begin(), b2.end(), b[0].begin());
    }

    std::cout << "chunk_size: " << chunk_size << std::endl;

    for (int i = 0; i < num_cu; i++)
    {
        if (!single_krnl)
        {
            OCL_CHECK(err, buffer_inout[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * b[i].size(), b[i].data(), &err));

            if (i == num_cu - 1)
            {
                OCL_CHECK(err, err = krnl_lu[i].setArg(5, nrhs - i * chunk_size));
            }
            else
            {
                OCL_CHECK(err, err = krnl_lu[i].setArg(5, chunk_size));
            }
        }
        else
        {
            OCL_CHECK(err, buffer_inout[i] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * b[i].size(), b[i].data(), &err));
            OCL_CHECK(err, err = krnl_lu[i].setArg(5, nrhs));
            break;
        }
    }

    for (int i = 0; i < num_cu; i++)
    {
        OCL_CHECK(err, err = krnl_lu[i].setArg(0, buffer_in0));
        OCL_CHECK(err, err = krnl_lu[i].setArg(1, buffer_in1));

        OCL_CHECK(err, err = krnl_lu[i].setArg(2, n));
        OCL_CHECK(err, err = krnl_lu[i].setArg(3, Numeric.lusize_sum));
        OCL_CHECK(err, err = krnl_lu[i].setArg(4, Symbolic.nblocks));

        OCL_CHECK(err, err = krnl_lu[i].setArg(6, buffer_inout[i]));

        // Copy input data to device global memory
        OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_inout[i]}, 0)); /* 0 means from host*/

        if (single_krnl)
            break;
    }
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in0, buffer_in1}, 0));
    OCL_CHECK(err, err = q.finish());

    for (int i = 0; i < num_cu; i++)
    {
        // Launch the kernel
        OCL_CHECK(err, err = q.enqueueTask(krnl_lu[i]));
        if (single_krnl)
            break;
    }
    OCL_CHECK(err, err = q.finish());

    // Copy Result from Device Global Memory to Host Local Memory
    for (int i = 0; i < num_cu; i++)
    {
        OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_inout[i]}, CL_MIGRATE_MEM_OBJECT_HOST));
        if (single_krnl)
            break;
    }
    OCL_CHECK(err, err = q.finish());
    // OPENCL HOST CODE AREA END

    // Compare the results of the Device to the simulation
    bool match = true;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nrhs; j++)
        {
            int idx = j / chunk_size < num_cu ? j / chunk_size : num_cu - 1;
            int idy = i + n * (j / chunk_size < num_cu ? j % chunk_size : j - chunk_size * (num_cu - 1));
            if (std::abs(b_cpu[i + n * j] - b[idx][idy]) > 1e-5)
            {
                std::cout << "Mismatched result x[" << i << "][" << j << "]: CPU x[" << i << "][" << j << "]=" << b_cpu[i + n * j] << ", FPGA x[" << idx << "][" << idy << "]=" << b[idx][idy] << std::endl;
                match = false;
            }
            else
            {
                if (i > 10)
                    continue;
                else if (j < nrhs - 1)
                    std::cout << "x[" << i << "," << j << "] = " << b[idx][idy] << " ";
                else
                    std::cout << "x[" << i << "," << j << "] = " << b[idx][idy] << std::endl;
            }
        }
    }

    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << std::endl;
    return (match ? EXIT_SUCCESS : EXIT_FAILURE);
}
