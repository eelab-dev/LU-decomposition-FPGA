#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

class Matrix
{
private:
    int size;
    std::vector<std::vector<double>> matrix;
    std::vector<int> sort_indexes(const std::vector<std::vector<double>> &matrix)
    {

        // initialize original index locations
        std::vector<int> idx(size);
        std::iota(idx.begin(), idx.end(), 0); // sequentially increasing from 0

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values
        stable_sort(idx.begin(), idx.end(),
                    [&matrix](size_t i1, size_t i2)
                    { return matrix[i1][0] > matrix[i2][0]; });

        return idx;
    }

public:
    Matrix(int Size)
    {
        size = Size;
        std::cout << "Constructor called" << std::endl;
        // for (int i = 0; i < size; i++)
        // {
        //     std::vector<double> temp(size);
        //     matrix.push_back(temp);
        // }
        matrix.resize(size, std::vector<double>(size));
    }

    void set(double value, int row, int column) { matrix[row][column] = value; }

    double get(int row, int column) { return matrix[row][column]; }

    void resize(int Size)
    {
        size = Size;
        matrix.clear();
        matrix.resize(size, std::vector<double>(size));
    }

    void randfill(int lower = 1, int upper = 10, int seed = std::time(0))
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<unsigned> distrib(lower, upper);
        for (int i{}; i < size; i++)
            for (int j{}; j < size; j++)
                matrix[i][j] = distrib(gen);
    }

    void LUdecomposition(std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u)
    {
        // Decomposing matrix into Upper and Lower
        // triangular matrix
        for (int i = 0; i < size; i++)
        {
            // Upper Triangular
            for (int k = i; k < size; k++)
            {
                // Summation of L(i, j) * U(j, k)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (l[i][j] * u[j][k]);

                // Evaluating U(i, k)
                u[i][k] = matrix[i][k] - sum;
            }

            // Lower Triangular
            for (int k = i; k < size; k++)
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

    void LUPivot(std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u)
    {
        int i = 0;
        std::vector<int> sorted = sort_indexes(matrix);

        for (int sort_i : sorted)
        {
            // Upper Triangular
            for (int k = i; k < size; k++)
            {
                // Summation of L(i, j) * U(j, k)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (l[i][j] * u[j][k]);

                // Evaluating U(i, k)
                u[i][k] = matrix[sort_i][k] - sum;
            }

            // Lower Triangular
            for (int k = i; k < size; k++)
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

    void display(std::string name = "a")
    {
        std::cout << name << "=" << std::endl;
        for (int i{}; i < matrix.size(); i++)
        {
            for (int j{}; j < matrix.size(); j++)
            {
                std::cout << std::setprecision(2) << std::fixed << std::setw(8) << std::left << matrix[i][j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    friend std::ostream &operator<<(std::ostream &stream, const Matrix &matrix)
    {
        for (int i{}; i < matrix.size; i++)
        {
            for (int j{}; j < matrix.size; j++)
            {
                stream << std::setprecision(2) << std::fixed << std::setw(8) << std::left << matrix.matrix[i][j];
            }
            stream << std::endl;
        }
        stream << std::endl;
        return stream;
    }
};

void display(std::vector<std::vector<double>> &vect, std::string name = "a")
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