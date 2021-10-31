#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;
const int SIZE = 10;
void LUdecomposition(vector<vector<double>> &a, vector<vector<double>> &l, vector<vector<double>> &u);
void LUPivot(vector<double> **a, vector<vector<double>> &l, vector<vector<double>> &u);
void display(vector<vector<double>> &vect, string name = "a");
void display2(vector<double> **a);

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v)
{

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2)
                { return v[i1] > v[i2]; });

    return idx;
}

int main()
{
    int seed = 2021;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 10);

    vector<vector<double>> matrix(SIZE, vector<double>(SIZE)),
        l(SIZE, vector<double>(SIZE)),
        u(SIZE, vector<double>(SIZE));

    for (int i{}; i < SIZE; i++)
        for (int j{}; j < SIZE; j++)
            matrix[i][j] = distrib(gen);

    vector<double> *pivot[SIZE];

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    int j{};
    for (int i : sort_indexes(matrix))
        pivot[j++] = &matrix[i];

    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    LUPivot(pivot, l, u);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    display(matrix, "Original");
    display2(pivot);
    display(l, "l");
    display(u, "u");

    cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << endl;

    return 0;
}

void display(vector<vector<double>> &vect, string name)
{
    cout << name << "=" << endl;
    for (int i{}; i < vect.size(); i++)
    {
        for (int j{}; j < vect.size(); j++)
        {
            cout << setprecision(6) << fixed << setw(10) << left << vect[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void LUdecomposition(vector<vector<double>> &a, vector<vector<double>> &l, vector<vector<double>> &u)
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

void display2(vector<double> **a)
{
    int n = (*a[0]).size();
    // cout << setprecision(2) << fixed << setw(8) << left << a[0][1] << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << setprecision(6) << fixed << setw(10) << left << (*a[i])[j];
        cout << endl;
    }
}

void LUPivot(vector<double> **a, vector<vector<double>> &l, vector<vector<double>> &u)
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    int n = (*a[0]).size();
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
            u[i][k] = (*a[i])[k] - sum;
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
                l[k][i] = ((*a[k])[i] - sum) / u[i][i];
            }
        }
    }
}