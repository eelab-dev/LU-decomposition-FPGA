#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;
const int SIZE = 1000;
void LUdecomposition(vector<vector<double>> &a, vector<vector<double>> &l, vector<vector<double>> &u);
void display(vector<vector<double>> &vect, string name = "a");

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

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    LUdecomposition(matrix, l, u);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // display(matrix);
    // display(l, "l");
    // display(u, "u");

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
            cout << setprecision(2) << fixed << setw(8) << left << vect[i][j];
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