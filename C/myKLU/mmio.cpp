#include <iostream>
#include "include/mmio.h"

template <typename T>
struct aligned_allocator
{
    using value_type = T;

    aligned_allocator() {}

    aligned_allocator(const aligned_allocator &) {}

    template <typename U>
    aligned_allocator(const aligned_allocator<U> &) {}

    T *allocate(std::size_t num)
    {
        void *ptr = nullptr;

#if defined(_WINDOWS)
        {
            ptr = _aligned_malloc(num * sizeof(T), 4096);
            if (ptr == nullptr)
            {
                std::cout << "Failed to allocate memory" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
#else
        {
            if (posix_memalign(&ptr, 4096, num * sizeof(T)))
                throw std::bad_alloc();
        }
#endif
        return reinterpret_cast<T *>(ptr);
    }
    void deallocate(T *p, std::size_t num)
    {
#if defined(_WINDOWS)
        _aligned_free(p);
#else
        free(p);
#endif
    }
};

int main(void)
{
    char filename[] = "../../Matrix_Sample/host.mtx";
    char bname[] = "../../Matrix_Sample/host_b.mtx";

    std::vector<int, aligned_allocator<int>> Ap, Ai;
    std::vector<double, aligned_allocator<double>> Ax, b;
    int n;
    if (read_sparse(filename, &n, Ap, Ai, Ax))
        return 1;

    if (read_bmatrix(bname, b))
    {
        printf("Error bmatrix\n");
        return 1;
    }

    std::vector<double> b2(b.data(), b.data() + n);

    for (int i = 0; i < Ap.size(); i++)
        std::cout << "Ap[" << i << "]=" << Ap[i] << std::endl;
    for (int i = 0; i < Ai.size(); i++)
        std::cout << "Ai[" << i << "]=" << Ai[i] << "\tAx[" << i << "]=" << Ax[i] << std::endl;

    // for (int i = 0; i < n; i++)
    // {
    //     if (std::abs(b2[i] - b[i]) > 1e-5)
    //     {
    //         std::cout << "Mismatched result x[" << i << "]! CPU x[" << i << "]=" << b2[i] << ", FPGA x[" << i << "]=" << b[i] << std::endl;
    //         break;
    //     }
    //     else
    //         std::cout << "x[" << i << "]=" << b[i] << std::endl;
    // }

    return 0;
}