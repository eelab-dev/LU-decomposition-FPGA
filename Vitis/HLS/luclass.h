#include <algorithm>
#include <ctime>
#include <iomanip>
// #include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <hls_vector.h>

const int SIZE=10;

class Matrix
{
private:
    int size;
    hls::vector<hls::vector<double,SIZE>,SIZE> matrix;
    hls::vector<int,SIZE> sort_indexes(const hls::vector<hls::vector<double,SIZE>,SIZE> &matrix)
    {

        // initialize original index locations
        hls::vector<int,SIZE> idx(SIZE);
        std::iota(idx.begin(), idx.end(), 0); // sequentially increasing from 0

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values
        std::stable_sort(idx.begin(), idx.end(),
                    [&matrix](size_t i1, size_t i2)
                    { return matrix[i1][0] > matrix[i2][0]; });

        return idx;
    }

public:
    Matrix(){}

    void set(double value, int row, int column) { matrix[row][column] = value; }

    double get(int row, int column) { return matrix[row][column]; }

    // void resize(int Size)
    // {
    //     size = Size;
    //     matrix.clear();
    //     matrix.resize(size, hls::vector<double>(size));
    // }

    void randfill(int lower = 1, int upper = 10, int seed = std::time(0))
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<unsigned> distrib(lower, upper);
        for (int i{}; i < SIZE; i++)
            for (int j{}; j < SIZE; j++)
                matrix[i][j] = distrib(gen);
    }

    void LUdecomposition(hls::vector<hls::vector<double,SIZE>,SIZE> &l, hls::vector<hls::vector<double,SIZE>,SIZE> &u)
    {
        // Decomposing matrix into Upper and Lower
        // triangular matrix
        for (int i = 0; i < SIZE; i++)
        {
            // Upper Triangular
            for (int k = i; k < SIZE; k++)
            {
                // Summation of L(i, j) * U(j, k)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (l[i][j] * u[j][k]);

                // Evaluating U(i, k)
                u[i][k] = matrix[i][k] - sum;
            }

            // Lower Triangular
            for (int k = i; k < SIZE; k++)
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
                    l[k][i] = (matrix[k][i] - sum) / u[i][i];
                }
            }
        }
    }

    void LUPivot(hls::vector<hls::vector<double,SIZE>,SIZE> &l, hls::vector<hls::vector<double,SIZE>,SIZE> &u)
    {
        int i = 0;
        hls::vector<int,SIZE> sorted = sort_indexes(matrix);

        for (int sort_i : sorted)
        {
            // Upper Triangular
            for (int k = i; k < SIZE; k++)
            {
                // Summation of L(i, j) * U(j, k)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (l[i][j] * u[j][k]);

                // Evaluating U(i, k)
                u[i][k] = matrix[sort_i][k] - sum;
            }

            // Lower Triangular
            for (int k = i; k < SIZE; k++)
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
                    l[k][i] = (matrix[sorted[k]][i] - sum) / u[i][i];
                }
            }
            i++;
        }
    }

    // void display(std::string name = "a")
    // {
    //     std::cout << name << "=" << std::endl;
    //     for (int i{}; i < matrix.size(); i++)
    //     {
    //         for (int j{}; j < matrix.size(); j++)
    //         {
    //             std::cout << std::setprecision(2) << std::fixed << std::setw(8) << std::left << matrix[i][j];
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // friend std::ostream &operator<<(std::ostream &stream, const Matrix &matrix)
    // {
    //     for (int i{}; i < matrix.size; i++)
    //     {
    //         for (int j{}; j < matrix.size; j++)
    //         {
    //             stream << std::setprecision(2) << std::fixed << std::setw(8) << std::left << matrix.matrix[i][j];
    //         }
    //         stream << std::endl;
    //     }
    //     stream << std::endl;
    //     return stream;
    // }
};

// void display(hls::vector<hls::vector<double>> &vect, std::string name = "a")
// {
//     std::cout << name << "=" << std::endl;
//     for (int i{}; i < vect.size(); i++)
//     {
//         for (int j{}; j < vect.size(); j++)
//         {
//             std::cout << std::setprecision(2) << std::fixed << std::setw(8) << std::left << vect[i][j];
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }
