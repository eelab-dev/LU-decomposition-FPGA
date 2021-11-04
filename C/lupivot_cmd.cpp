#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

void LUPivot(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u);
void display(std::vector<std::vector<double>> &vect, std::string name = "a");

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v)
{

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0); // sequentially increasing from 0

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2)
                { return v[i1][0] > v[i2][0]; });

    return idx;
}

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

void display(std::vector<std::vector<double>> &vect, std::string name)
{
    std::cout << name << "=" << std::endl;
    for (int i{}; i < vect.size(); i++)
    {
        for (int j{}; j < vect.size(); j++)
        {
            std::cout << std::setprecision(6) << std::fixed << std::setw(10) << std::left << vect[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void LUPivot(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u)
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    int n = a.size();
    int i = 0;
    std::vector<size_t> sorted = sort_indexes(a);

    for (int sort_i : sorted)
    {
        // Upper Triangular
        for (int k = i; k < n; k++)
        {
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (l[i][j] * u[j][k]);

            // Evaluating U(i, k)
            u[i][k] = a[sort_i][k] - sum;
        }

        // Lower Triangular
        for (int k = i; k < n; k++)
        {
            if (i == k)
                l[i][i] = 1; // Diagonal as 1
            else
            {
                // Summation of L(k, j) * U(j, i)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (l[k][j] * u[j][i]);

                // Evaluating L(k, i)
                l[k][i] = (a[sorted[k]][i] - sum) / u[i][i];
            }
        }
        i++;
    }
}