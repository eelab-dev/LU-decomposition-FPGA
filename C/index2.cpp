#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

const int SIZE = 4;

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
    int seed = 2022;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 100);

    vector<vector<double>> matrix(SIZE, vector<double>(SIZE));
    vector<double> col1(SIZE);
    vector<double> *p[SIZE];

    for (int i{}; i < SIZE; i++)
        for (int j{}; j < SIZE; j++)
        {
            matrix[i][j] = distrib(gen);
            cout << setw(4) << left << matrix[i][j];
        }
    cout << endl;

    for (int i{}; i < SIZE; i++)
        col1[i] = matrix[i][0];

    int j{};
    for (int i : sort_indexes(matrix))
    {
        p[j++] = &matrix[i];
    }

    display(matrix);
    for (int i{}; i < SIZE; i++)
    {
        for (int j{}; j < SIZE; j++)
            cout << setprecision(2) << fixed << setw(8) << left << (*p[i])[j];
        cout << endl;
    }
    cout << endl;

    display2(p);

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

void display2(vector<double> **a)
{
    int n = (*a[0]).size();
    // cout << setprecision(2) << fixed << setw(8) << left << a[0][1] << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << setprecision(2) << fixed << setw(8) << left << (*a[i])[j];
        cout << endl;
    }
}