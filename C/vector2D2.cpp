// Unit Lower matrix

#include "LU/lu.h"
#include <chrono>
#include <ctime>
#include <random>

const int SIZE = 39;

int main()
{
    int seed = 2021;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 10);

    std::vector<std::vector<double>> matrix(SIZE, std::vector<double>(SIZE)),
        l(SIZE, std::vector<double>(SIZE)),
        u(SIZE, std::vector<double>(SIZE));

    for (int i{}; i < SIZE; i++)
        for (int j{}; j < SIZE; j++)
            matrix[i][j] = distrib(gen);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    LUdecomposition(matrix, l, u);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // display(matrix);
    // display(l, "l");
    // display(u, "u");

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;

    return 0;
}