#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<std::vector<T>> &v)
{

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0); // sequentially increasing from 0

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort(idx.begin(), idx.end(),
                     [&v](size_t i1, size_t i2)
                     { return v[i1][0] > v[i2][0]; });

    return idx;
}

void LUdecomposition(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &l, std::vector<std::vector<double>> &u)
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    int n = a.size();
    for (int i = 0; i < n; i++)
    {
        // Upper Triangular
        for (int k = i; k < n; k++)
        {
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (l[i][j] * u[j][k]);

            // Evaluating U(i, k)
            u[i][k] = a[i][k] - sum;
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
                l[k][i] = (a[k][i] - sum) / u[i][i];
            }
        }
    }
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
