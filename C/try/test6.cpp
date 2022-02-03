#include <iostream>
#include <malloc.h>

typedef struct
{
    int *a;
    float *b;
} klu;

klu *func(int size_a, int size_b);

int main(void)
{
    klu *test;
    int size_a = 3, size_b = 4;
    test = func(size_a, size_b);
    for (int i = 0; i < size_a; i++)
        std::cout << "a[" << i << "]=" << test->a[i] << std::endl;
    for (int i = 0; i < size_b; i++)
        std::cout << "b[" << i << "]=" << test->b[i] << std::endl;

    return 0;
}

klu *func(int size_a, int size_b)
{
    int c[size_a];
    float d[size_b];
    klu *test;
    test->a = (int *)malloc(size_a * sizeof(int));
    test->b = (float *)malloc(size_b * sizeof(double));
    for (int i = 0; i < size_a; i++)
        test->a[i] = 2 * i;
    for (int i = 0; i < size_b; i++)
        test->b[i] = 1.5 * i;

    // test->a = c;
    // test->b = d;
    return test;
}