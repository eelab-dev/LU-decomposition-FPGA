#include "lu.h"
#include <chrono>
#include <ctime>
#include <random>

int main(int argc, char *argv[])
{
    int seed = std::time(0);
    // int seed=2021;
    int SIZE = 10;

    if (argc > 2)
        if (std::string{"--size"}.compare(argv[1]) == 0 || std::string{"-t"}.compare(argv[1]) == 0)
            SIZE = std::stoi(argv[2]);
    std::cout << "argc:" << argc << " size:" << SIZE << std::endl;

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
    LUPivot(matrix, l, u);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // display(matrix, "Original");
    // display(l, "l");
    // display(u, "u");
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;

    return 0;
}