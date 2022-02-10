# 29/12/2021
[Online Tutorial](https://xilinxcustomertraining.netexam.com/?tkn=92697314-9F12-463F-8268-80D79E2ED373#:ct54538)


# 11/01/2022
## Error
```log
INFO: [IP_Flow 19-1686] Generating 'Simulation' target for IP 'lu_kernel_sitofp_32ns_32_5_no_dsp_1_ip'...
ERROR: '2201112057' is an invalid argument. Please specify an integer value.
    while executing
"rdi::set_property core_revision 2201112057 {component component_1}"
    invoked from within
"set_property core_revision $Revision $core"
    (file "run_ippack.tcl" line 1612)
INFO: [Common 17-206] Exiting Vivado at Tue Jan 11 20:57:59 2022...
ERROR: [IMPL 213-28] Failed to generate IP.
INFO: [HLS 200-111] Finished Command export_design CPU user time: 23.4 seconds. CPU system time: 1.52 seconds. Elapsed time: 48.34 seconds; current allocated memory: 7.418 MB.
```

**[Solution](https://support.xilinx.com/s/article/76960?language=en_US)**


# Adding Libraries and Library Paths to Vitis
You can add libraries and library paths for Application projects. If you have a custom library to link against, you should specify the library path and the library name to the linker.

To set properties for your Application project:

1. Right-click your Application project (Host one) and select **C/C++ Build Settings**. Alternatively, select **Properties** and navigate to **C/C++ Build > Settings**.
![Host](Resources/Add%20library.png)
2. Expand the target linker section and select the libraries to which you want to add the custom library path and library name.
![Add .so](Resources/Addlink.png)

3. We can also add the header search path from the option *GCC Host Compiler/Includes*.


# 12/01/2022
## Runtime Initialization & Device Configuration
```cpp
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
```

## Buffer Allocation
```cpp
// Allocate Buffer in Global Memory
// Buffers are allocated using CL_MEM_USE_HOST_PTR for efficient memory and
// Device-to-host communication
OCL_CHECK(err, cl::Buffer buffer_in1(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, vector_size_bytes,
                                        host_matrix.data(), &err));
OCL_CHECK(err, cl::Buffer buffer_output1(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, vectorout_size_bytes,
                                            source_hw_l.data(), &err));
OCL_CHECK(err, cl::Buffer buffer_output2(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, vectorout_size_bytes,
                                            source_hw_u.data(), &err));
```

## Setup Kernel Argument
```cpp
OCL_CHECK(err, err = krnl_vector_add.setArg(0, buffer_in1));
OCL_CHECK(err, err = krnl_vector_add.setArg(1, buffer_output1));
OCL_CHECK(err, err = krnl_vector_add.setArg(2, buffer_output2));
OCL_CHECK(err, err = krnl_vector_add.setArg(3, size));
```

## Writing Buffers to FPGA Memory
```cpp
// Copy input data to device global memory
OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_in1}, 0 /* 0 means from host*/));
```

## Running the Accelerators
```cpp
cl::Event event;
uint64_t nstimestart, nstimeend;

// Launch the Kernel
// For HLS kernels global and local size is always (1,1,1). So, it is
// recommended
// to always use enqueueTask() for invoking HLS kernel
OCL_CHECK(err, err = q.enqueueTask(krnl_vector_add, nullptr, &event));
```
### Running the Accelerators (Out-of-order)
```cpp
// Alternative out-of-order queue events control
cl::Event event_sp;
q.enqueueTask(krnl_vector_add, NULL, &event_sp);
clWaitForEvents(1, (const cl_event *) &event_sp);
```

## Reading Buffers from FPGA Mem
```cpp
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
```

## FPGA Acceleration
See Host Code [here](C/vitis/host.cpp).

See Kernel Code [here](C/vitis/lu_kerenl.cpp).

# 13/01/2022
## Gilbert-Peierls' algorithm
Gilbert-Peierls left-looking algorithm factors the matrix column-by-column from left to right.
```cpp
L=I; // I=identity matrix
for k=1:N
b = A(:,k); // kth column of A
x = L \ b; // \ is Lx=b solve
U(1:k) = x(1:k);
L(k+1:N) = x(k+1:N) / U(k,k);
end;
// return L and U as result
```

For Sparse matrix, we can implement the algorithm with the following code. See [here](C/sparse.cpp).
```cpp
void SparseLU(int *Ap, int *Ai, double *Ax, int *Lp, int *Li, double *Lx, int *Up, int *Ui, double *Ux)
{
	int countL = 1, countU = 0;

	for (int i = 0; i < squareSize; i++)
	{
		Lp[i] = i;
		Li[i] = i;
		Lx[i] = 1;
	}

	for (i = 0; i < squareSize; i++)
	{
		double b[squareSize], x[squareSize];
		for (int j = 0; j < squareSize; j++)
		{
			b[j] = 0;
			x[j] = 0;
		}

		std::cout << i << std::endl;
		if (Ap[i] != Ap[i + 1])
		{
			for (int j = Ap[i]; j < Ap[i + 1]; j++)
				b[Ai[j]] = Ax[j];
		}
		else
		{
			// Lp[i] = 0;
			Up[i] = 0;
			continue;
		}
		for (int j = 0; j < squareSize; j++)
		{
			for (int k = 0; k <= j; k++)
			{
				double sum = 0;
				for (int t = Lp[k]; t < Lp[k + 1]; t++)
					if (Li[t] == j)
						sum += Lx[Li[t]] * x[k];
			}
			x[j] = b[j] - sum;
		}

		for (int j = 0; j <= i; j++)
		{
			if (std::abs(x[j]) < 1e-6)
				continue;
			else
			{
				Ui[countU] = j;
				Ux[countU++] = x[j];
			}
		}
		Up[i + 1] = countU;

		for (int j = i + 1; j < squareSize; j++)
		{
			if (std::abs(x[j]) < 1e-6)
				continue;
			else
			{
				Li[countL] = j;
				Lx[countL++] = x[j] / Ux[Up[i + 1] - 1];
			}
		}
		Lp[i + 1] = countL;
		if (i < squareSize - 1)
		{
			Li[countL] = i + 1;
			Lx[countL++] = 1;
		}
	}
}
```



## Block LU Decomposition
$$
\left[ \begin{matrix}
	A_{11}&		A_{12}\\
	A_{21}&		A_{22}\\
\end{matrix} \right] =\left[ \begin{matrix}
	L_{11}&		0\\
	L_{21}&		L_{22}\\
\end{matrix} \right] \times \left[ \begin{matrix}
	U_{11}&		U_{12}\\
	0&		U_{22}\\
\end{matrix} \right]
$$


# 15/01/2021
## Improve Gilbert-Peierls' algorithm
As the elements in the **L** and **U** is fixed, we can preallocate the indexes using the function below.

See [here](C/sparse2.cpp) for more detail.
```cpp
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
```

## Implement in Vitis
[Host code](C/vitis/sparse/host.cpp).

[Kernel code](C/vitis/sparse/sparse.cpp).



# 21/01/2022
## Benchmark
**[File here](C/sparse3.c)**

First, we use a small  5x5 sparse [matrix](Matrix_Sample/test.mtx) with 12 nonzeros to test.
Run for ten times and it takes about **0.9us** to decompose on average.

Then, we use a 10x10 sparse [matrix](Matrix_Sample/test4.mtx) with 24 to test.
Run for ten times and it takes about **2.9us** to decompose on average.

Then, we try a little large matrix [HB/bcspwr01](https://sparse.tamu.edu/HB/bcspwr01). It is a 39x39 matrix with 131 nonzeros.
Run for ten times and it takes about **61.5us** to decompose on average. It is much longer that the 5x5 sparse matrix.

![bcspwr01](https://suitesparse-collection-website.herokuapp.com/files/HB/bcspwr01.png)

**[With KLU](C/KLU/kluldemo.c)**

KLU itself is more stable. For the 5x5 sparse matrix, it takes about **2.3us** to decompose on average.

For the 10x10 sparse matrix, it takes about **2.9us** to decompose on average.

For the *HB/bcspwr01* matrix, it takes about **3us**. Only 0.3us longer that the size of 5.


# 22/01/2022
## Emulation on FPGA
The reason that *the simulation exited unexpectedly* previously is that the memory migration error.
As the kernel requires 10 arguments in total, shown below, not only the first four arguments are initialized, part of the other parameters are also initialized, which also are needed to transfer to the kernel.
```c
void sparselu(int squareSize, int *AP, int *AI, double *AX, int *LP, int *LI, double *LX, int *UP, int *UI, double *UX)
```

Therefore, the memory migration should be written as:
```cpp
OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_output1, buffer_output2, buffer_output3, buffer_output4, buffer_output5, buffer_output6}, CL_MIGRATE_MEM_OBJECT_HOST));
```

## Emulation Result
![Kernel](Resources/emulation10.png)
![Kernel](Resources/emulation10kerneldata.png)
![Kernel](Resources/emulation10hls.png)

# 25/01/2022
## Improved - Reachability
The main reason that what I have written is much slower is that it does not know the data index in a row. Therefore, when implementing the Gilbert-Peierls, especially when solving
$$
x=L\backslash A(:,k)
$$
Much of the time is spent on searching for data when may exist in a certain row.

Therefore, to solve this problem, we need to have a knowledge about how the data is distributed.

![Reach](Resources/klureach.png)

For example, if we have a matrix like this

![Matrix](Resources/klumatrix.jpg)

We can derive its reachability

![Reach Diagram](Resources/klureachdiagram.png)

# 28/01/2022
## Improved timing method
Using `<chrono>` timing method in C++ instead of `<time.h>`
is C, Benchmark the code.

For the [HB/ash85](https://sparse.tamu.edu/HB/ash85) matrix, which is 85x85 with 523 entries, it will take my code about 30us to finish.

![HB/ash85](https://suitesparse-collection-website.herokuapp.com/files/HB/ash85.png)

For KLU, it will take about 3.2us to finish.

This result is better than the timing result is `C`.

# 29/01/2022
## Unroll loop in Kernel
Loop Unrolling
The compiler can also unroll a loop, either partially or completely to perform multiple loop
iterations in parallel. This is done using the pragma HLS unroll. Unrolling a loop can lead to a very
fast design, with significant parallelism. However, because all the operations of the loop
iterations are executed in parallel, a large amount of programmable logic resource are required to
implement the hardware. As a result, the compiler can face challenges dealing with such a large
number of resources and can face capacity problems that slow down the kernel compilation
process. It is a good guideline to unroll loops that have a small loop body, or a small number of
iterations.

```c
vadd: for(int i = 0; i < 20; i++) {
#pragma HLS UNROLL
c[i] = a[i] + b[i];
}
```

In the preceding example, you can see `pragma HLS UNROLL` has been inserted into the body of
the loop to instruct the compiler to unroll the loop completely. All 20 iterations of the loop are
executed in parallel if that is permitted by any data dependency.

Completely unrolling a loop can consume significant device resources, while partially unrolling the
loop provides some performance improvement while using fewer hardware resources.

**My code**

In my code, I apply the `unroll` pragma during the data reading from the global memory as well as the writing to memory process.

In addition, I also apply `unroll` for the last two steps of the Gilbert-Peierls' algorithm, i.e.
```cpp
L=I; // I=identity matrix
for k=1:N
b = A(:,k); // kth column of A
x = L \ b; // \ is Lx=b solve
U(1:k) = x(1:k); 	<--
L(k+1:N) = x(k+1:N) / U(k,k); <--
end;
// return L and U as result
```

With unroll, we can see that the kernel time decreases from 71us to 50us.

![With Unroll](Resources/emulation10%20unroll.png)

# 02/02/2022
## Further deep in KLU
To implement KLU factorization in FPGA, I find that the function
```c
static void factor2
(
    /* inputs, not modified */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Entry Ax [ ],
    KLU_symbolic *Symbolic,

    /* inputs, modified on output: */
    KLU_numeric *Numeric,
    KLU_common *Common
)
```

mainly performs the factorization process in KLU.

This function mainly takes 6 parameters, where the first three is the input matrix, `KLU_symbolic *Symbolic` is the result from symbolic analysis, `KLU_numeric *Numeric` is the factorization result, and ` KLU_numeric *Numeric` stores the control parameters of the factorization.


# 06/02/2022
[PLRAM](https://xilinx.github.io/Vitis_Accel_Examples/2020.2/html/plram_access.html)