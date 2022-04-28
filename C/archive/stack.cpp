#include <iostream>
#include <vector>
#include <chrono>

void teststack(int size);
int main()
{
    std::string i;
    std::cin >> i;
    int size = std::stoi(i);
    teststack(size);
    // double a[size], b[size], c[size], d[size], e[size], f[size], g[size];
    return 0;
}

void teststack(int size)
{
    std::chrono::steady_clock::time_point begin, stop;
    begin = std::chrono::steady_clock::now();
    std::vector<int> a;
    a.reserve(size);
    std::vector<int> b;
    b.reserve(size);
    std::vector<int> c;
    c.reserve(size);
    std::vector<int> d;
    d.reserve(size);
    std::vector<int> e;
    e.reserve(size);
    std::vector<int> f;
    f.reserve(size);
    std::vector<int> g;
    g.reserve(size);
    stop = std::chrono::steady_clock::now();

    std::cout << "Creation time: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - begin).count() << "µs" << std::endl;

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < size; i++)
    {
        // a[i] = i;
        // b[i] = i + 1;
        // c[i] = i + 2;
        // d[i] = i + 3;
        // e[i] = i + 4;
        // f[i] = i + 5;
        // g[i] = i + 6;
        a.push_back(i);
        b.push_back(i + 1);
        c.push_back(i + 2);
        d.push_back(i + 3);
        e.push_back(i + 4);
        f.push_back(i + 5);
        g.push_back(i + 6);
    }
    stop = std::chrono::steady_clock::now();

    printf("Size of a: %d\n", a.size());

    std::cout << "Running time: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - begin).count() << "µs" << std::endl;

    double *aa, *bb, *cc, *dd, *ee, *ff, *gg;
    begin = std::chrono::steady_clock::now();
    aa = (double *)malloc(size * sizeof(double));
    bb = (double *)malloc(size * sizeof(double));
    cc = (double *)malloc(size * sizeof(double));
    dd = (double *)malloc(size * sizeof(double));
    ee = (double *)malloc(size * sizeof(double));
    ff = (double *)malloc(size * sizeof(double));
    gg = (double *)malloc(size * sizeof(double));
    stop = std::chrono::steady_clock::now();

    std::cout << "Creation time: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - begin).count() << "µs" << std::endl;

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < size; i++)
    {
        aa[i] = i;
        bb[i] = i + 1;
        cc[i] = i + 2;
        dd[i] = i + 3;
        ee[i] = i + 4;
        ff[i] = i + 5;
        gg[i] = i + 6;
    }
    stop = std::chrono::steady_clock::now();

    std::cout << "Running time: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - begin).count() << "µs" << std::endl;
}