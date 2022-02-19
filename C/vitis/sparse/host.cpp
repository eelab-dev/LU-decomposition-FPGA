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
#include <vector>

void analyse(int *Ap, int *Ai, int *Up, int *Ui, int *Lp, int *Li, double *Lx, int *lnz, int *unz, int n);

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
    std::vector<int, aligned_allocator<int>> Ap{0, 2, 4, 7, 9, 13, 14, 18, 21, 23, 24};
    std::vector<int, aligned_allocator<int>> Ai{0, 7, 1, 5, 0, 2, 7, 3, 8, 3, 4, 5, 8, 5, 0, 6, 7, 8, 4, 5, 7, 2, 8, 9};
    std::vector<double, aligned_allocator<double>> Ax{8, 8, 2, 4, 10, 3, 10, 1, 5, 2, 1, 4, 10, 2, 3, 5, 3, 15, 4, 16, 3, 7, 2, 9};

    int n = Ap.size();
    int data_size = Ai.size();
    size_t p_size_bytes = sizeof(int) * n;
    size_t i_size_bytes = sizeof(int) * data_size;
    size_t x_size_bytes = sizeof(double) * data_size;

    std::vector<int, aligned_allocator<int>> Lp(n);
    std::vector<int, aligned_allocator<int>> Li(data_size);
    std::vector<double, aligned_allocator<double>> Lx(data_size);
    std::vector<int, aligned_allocator<int>> Up(n);
    std::vector<int, aligned_allocator<int>> Ui(data_size);
    std::vector<double, aligned_allocator<double>> Ux(data_size);

    int lnz = 0, unz = 0;
    analyse(Ap.data(), Ai.data(), Up.data(), Ui.data(), Lp.data(), Li.data(), Lx.data(), &lnz, &unz, n - 1);

    for (int i = 0; i < n; i++)
        printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);
    for (int i = n + 1; i < lnz; i++)
        printf("Li[%d]=%d\tLx[%d]=%lf\n", i, Li[i], i, Lx[i]);
    for (int i = n + 1; i < unz; i++)
        printf("Ui[%d]=%d\tUx[%d]=%lf\n", i, Ui[i], i, Ux[i]);

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
            OCL_CHECK(err, krnl_vector_add = cl::Kernel(program, "sparselu", &err));
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
    OCL_CHECK(err, cl::Buffer buffer_in1(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, p_size_bytes, Ap.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in2(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, i_size_bytes, Ai.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_in3(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, x_size_bytes, Ax.data(), &err));

    OCL_CHECK(err, cl::Buffer buffer_output1(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, p_size_bytes, Lp.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output2(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, i_size_bytes, Li.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output3(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, x_size_bytes, Lx.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output4(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, p_size_bytes, Up.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output5(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, i_size_bytes, Ui.data(), &err));
    OCL_CHECK(err, cl::Buffer buffer_output6(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, x_size_bytes, Ux.data(), &err));

    std::cout << "Buffer set" << std::endl;

    OCL_CHECK(err, err = krnl_vector_add.setArg(0, n - 1));
    OCL_CHECK(err, err = krnl_vector_add.setArg(1, buffer_in1));
    OCL_CHECK(err, err = krnl_vector_add.setArg(2, buffer_in2));
    OCL_CHECK(err, err = krnl_vector_add.setArg(3, buffer_in3));
    OCL_CHECK(err, err = krnl_vector_add.setArg(4, buffer_output1));
    OCL_CHECK(err, err = krnl_vector_add.setArg(5, buffer_output2));
    OCL_CHECK(err, err = krnl_vector_add.setArg(6, buffer_output3));
    OCL_CHECK(err, err = krnl_vector_add.setArg(7, buffer_output4));
    OCL_CHECK(err, err = krnl_vector_add.setArg(8, buffer_output5));
    OCL_CHECK(err, err = krnl_vector_add.setArg(9, buffer_output6));

    std::cout << "Arg set" << std::endl;

    // Copy input data to device global memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in1, buffer_in2, buffer_in3}, 0 /* 0 means from host*/));

    std::cout << "Mem set" << std::endl;

    // Launch the Kernel
    // For HLS kernels global and local size is always (1,1,1). So, it is
    // recommended
    // to always use enqueueTask() for invoking HLS kernel
    OCL_CHECK(err, err = q.enqueueTask(krnl_vector_add));

    std::cout << "Task set" << std::endl;

    // Copy Result from Device Global Memory to Host Local Memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_output1, buffer_output2, buffer_output3, buffer_output4, buffer_output5, buffer_output6}, CL_MIGRATE_MEM_OBJECT_HOST));
    q.finish();

    std::cout << "Mem back" << std::endl;
    // OPENCL HOST CODE AREA END

    // Compare the results of the Device to the simulation
    bool match = true;
    // for (int i = 0; i < DATA_SIZE; i++)
    // {
    //     if (source_hw_results[i] != source_sw_results[i])
    //     {
    //         std::cout << "Error: Result mismatch" << std::endl;
    //         std::cout << "i = " << i << " CPU result = " << source_sw_results[i]
    //                   << " Device result = " << source_hw_results[i] << std::endl;
    //         match = false;
    //         break;
    //     }
    // }
    for (int i = 0; i < n; i++)
        printf("Lp[%d]=%d\tLi[%d]=%d\tLx[%d]=%lf\tUp[%d]=%d\tUi[%d]=%d\tUx[%d]=%lf\n", i, Lp[i], i, Li[i], i, Lx[i], i, Up[i], i, Ui[i], i, Ux[i]);
    for (int i = n + 1; i < lnz; i++)
        printf("Li[%d]=%d\tLx[%d]=%lf\n", i, Li[i], i, Lx[i]);
    for (int i = n + 1; i < unz; i++)
        printf("Ui[%d]=%d\tUx[%d]=%lf\n", i, Ui[i], i, Ux[i]);

    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << std::endl;
    return (match ? EXIT_SUCCESS : EXIT_FAILURE);
}

void analyse(int *Ap, int *Ai, int *Up, int *Ui, int *Lp, int *Li, double *Lx, int *lnz, int *unz, int n)
{
    for (int i = 0, count = 0; i < n; i++)
    {
        for (int j = Ap[i]; j < Ap[i + 1]; j++)
        {
            if (Ai[j] < i)
            {
                Ui[*unz] = Ai[count];
                (*unz)++;
            }
            else if (Ai[j] == i)
            {
                Li[*lnz] = i;
                Lx[(*lnz)++] = 1;
                Ui[(*unz)++] = i;
            }
            else
            {
                Li[*lnz] = Ai[count];
                (*lnz)++;
            }
            count++;
        }
        Lp[i + 1] = *lnz;
        Up[i + 1] = *unz;
    }
}