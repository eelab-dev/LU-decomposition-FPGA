#include "LU/lu.h"
#include <chrono>
#include <fstream>
#include <random>

const int TIME = 10; // Sample time

int main(int argc, char *argv[])
{
    std::string filename("data_temp.csv");
    int seed = std::time(0);
    // int seed=2021;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    std::ofstream data("data/" + filename);

    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 10);

    std::vector<std::vector<double>> matrix, l, u;

    for (int size = 10; size <= 1000; size *= 10)
    // for (int size : {10, 50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000})
    {
        matrix.clear();
        l.clear();
        u.clear();
        matrix.resize(size, std::vector<double>(size));
        l.resize(size, std::vector<double>(size));
        u.resize(size, std::vector<double>(size));
        for (int i = 0; i < TIME; i++)
        {
            for (int j{}; j < size; j++)
                for (int k{}; k < size; k++)
                    matrix[j][k] = distrib(gen);

            begin = std::chrono::steady_clock::now();
            LUPivot(matrix, l, u);
            end = std::chrono::steady_clock::now();

            std::cout << "size: " << size << ",Time difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
            data << "size: " << size << ",Time difference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
        }
    }

    data.close();
    return 0;
}