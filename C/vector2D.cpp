#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;
const int SIZE = 100;
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

    display(matrix);
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
            cout << setprecision(2) << fixed << setw(8) << left << vect[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void LUdecomposition(vector<vector<double>> &a, vector<vector<double>> &l, vector<vector<double>> &u)
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