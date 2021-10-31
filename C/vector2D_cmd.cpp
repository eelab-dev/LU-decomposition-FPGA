// Unit upper matrix

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

void LUdecomposition(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u);
void display(std::vector<std::vector<double>> &vect, std::string name = "a");

int main(int argc, char *argv[])
{
    // int seed = 2021;
    int seed = std::time(0);
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
    LUdecomposition(matrix, l, u);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // display(matrix);
    // display(l, "l");
    // display(u, "u");

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;

    return 0;
}

void display(std::vector<std::vector<double>> &vect, std::string name)
{
    std::cout << name << "=" << std::endl;
    for (int i{}; i < vect.size(); i++)
    {
        for (int j{}; j < vect.size(); j++)
        {
            std::cout << std::setprecision(2) << std::fixed << std::setw(8) << std::left << vect[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void LUdecomposition(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u)
{
    int i = 0, j = 0, k = 0, n = a.size();
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}