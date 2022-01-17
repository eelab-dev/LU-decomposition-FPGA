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
#include <random>
#include "lu.h"
#include "klu.h"

#define SIZE 10
#define DATA_SIZE SIZE *SIZE

int gen_random()
{
    static std::default_random_engine e;
    static std::uniform_int_distribution<int> dist(1, 10);

    return dist(e);
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    int n = 5;
    int Ap[] = {0, 2, 5, 9, 10, 12};
    int Ai[] = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
    double Ax[] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
    double b[] = {8., 45., -3., 3., 19.};

    klu_symbolic *Symbolic;
        klu_numeric *Numeric;
        klu_common Common;
        int i;
        klu_defaults(&Common);
        Symbolic = klu_analyze(n, Ap, Ai, &Common);
        Numeric = klu_factor(Ap, Ai, Ax, Symbolic, &Common);
        klu_solve(Symbolic, Numeric, 5, 1, b, &Common);
        klu_free_symbolic(&Symbolic, &Common);
        klu_free_numeric(&Numeric, &Common);
        for (i = 0; i < n; i++)
            printf("x [%d] = %g\n", i, b[i]);

    std::string binaryFile = argv[1];
    size_t vector_size_bytes = sizeof(int) * DATA_SIZE;
    size_t vectorout_size_bytes = sizeof(float) * DATA_SIZE;
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
    std::vector<int, aligned_allocator<int>> host_matrix(DATA_SIZE);
    std::vector<float, aligned_allocator<float>> source_hw_l(DATA_SIZE);
    std::vector<float, aligned_allocator<float>> source_hw_u(DATA_SIZE);
    std::vector<float, aligned_allocator<float>> source_sw_l(DATA_SIZE, 0);
    std::vector<float, aligned_allocator<float>> source_sw_u(DATA_SIZE, 0);

    // Create the test data
    std::generate(host_matrix.begin(), host_matrix.end(), gen_random);
    std::cout << "A:\n";
    print(host_matrix.data(), SIZE, SIZE);
    LUdecomposition(host_matrix, source_sw_l, source_sw_u, SIZE);
    std::cout << "L:\n";
    print(source_sw_l.data(), SIZE, SIZE);
    std::cout << "U:\n";
    print(source_sw_u.data(), SIZE, SIZE);

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
            OCL_CHECK(err, krnl_vector_add = cl::Kernel(program, "lu_kernel", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device)
    {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Allocate Buffer in Global Memory
    // Buffers are allocated using CL_MEM_USE_HOST_PTR for efficient memory and
    // Device-to-host communication
    OCL_CHECK(err, cl::Buffer buffer_in1(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, vector_size_bytes,
                                         host_matrix.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output1(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, vectorout_size_bytes,
                                             source_hw_l.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output2(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, vectorout_size_bytes,
                                             source_hw_u.data(), &err));

    int size = SIZE;
    OCL_CHECK(err, err = krnl_vector_add.setArg(0, buffer_in1));
    OCL_CHECK(err, err = krnl_vector_add.setArg(1, buffer_output1));
    OCL_CHECK(err, err = krnl_vector_add.setArg(2, buffer_output2));
    OCL_CHECK(err, err = krnl_vector_add.setArg(3, size));

    // Copy input data to device global memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in1}, 0 /* 0 means from host*/));

    cl::Event event;
    uint64_t nstimestart, nstimeend;

    // Launch the Kernel
    // For HLS kernels global and local size is always (1,1,1). So, it is
    // recommended
    // to always use enqueueTask() for invoking HLS kernel
    OCL_CHECK(err, err = q.enqueueTask(krnl_vector_add, nullptr, &event));

    // Copy Result from Device Global Memory to Host Local Memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_output1, buffer_output2}, CL_MIGRATE_MEM_OBJECT_HOST));
    q.finish();
    OCL_CHECK(err, err = event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_START, &nstimestart));
    OCL_CHECK(err, err = event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_END, &nstimeend));
    // OPENCL HOST CODE AREA END
    auto lu_time = nstimeend - nstimestart;

    // Compare the results of the Device to the simulation
    std::cout << "L:\n";
    print(source_hw_l.data(), SIZE, SIZE);
    std::cout << "U:\n";
    print(source_hw_u.data(), SIZE, SIZE);
    bool match = true;
    for (int i = 0; i < DATA_SIZE; i++)
    {
        if (std::abs(source_hw_l[i] - source_sw_l[i]) > 1e-5)
        {
            std::cout << "Error: Result mismatch L" << std::endl;
            std::cout << "i = " << i << " CPU result = " << std::setprecision(10) << source_sw_l[i]
                      << " Device result = " << source_hw_l[i] << std::endl;
            match = false;
            break;
        }
        else if (std::abs(source_hw_u[i] - source_sw_u[i]) > 1e-5)
        {
            std::cout << "Error: Result mismatch U" << std::endl;
            std::cout << "i = " << i << " CPU result = " << std::setprecision(10) << source_sw_u[i]
                      << " Device result = " << source_hw_u[i] << std::endl;
            match = false;
            break;
        }
    }

    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << std::endl;

    std::cout << "| " << std::left << std::setw(24) << "lu time: "
                  << "|" << std::right << std::setw(24) << lu_time << " |\n";

    return (match ? EXIT_SUCCESS : EXIT_FAILURE);
}
