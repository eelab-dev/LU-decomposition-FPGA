#include <iostream>
#include <vector>

void fa(int *a)
{
    for (int i = 0; i < 10; i++)
        a[i] = 2 * i;
}

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
    int data[] = {1, 2, 3, 4, 5, 6};
    // std::vector<int, aligned_allocator<int>> a(data, data + 6);
    std::vector<int> a(data, data + 6);

    a.insert(a.begin() + 4, 13);
    a.insert(a.begin() + 2, 15);

    a.reserve(10);
    a[8] = 38;

    for (int i = 0; i < a.capacity(); i++)
        printf("a[%d]=%d\n", i, a[i]);

    std::cout << "size: " << a.size() << std::endl;
    a.resize(10);
    std::cout << "size: " << a.size() << std::endl;

    double LU[] = {0, 1, 1, 0, 0, -1};

    int i = 0;
    double *xp = LU + i;
    int *Li;
    Li = (int *)xp;
    std::cout << "lu[" << i << "]=" << LU[0] << std::endl;
    std::cout << "xp[" << i << "]=" << xp[0] << std::endl;
    std::cout << "Li[" << i << "]=" << Li[2] << std::endl;

    return 0;
}