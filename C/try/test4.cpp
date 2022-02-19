#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;

template <typename T>
void sort_indexes(const vector<T> &v)
{
    bool x = (v[0] < v[1]);
    std::cout << x;
    // for (int i = 0; i < v.size() - 1; i++)
    //     cout << v[i] > v[i + 1] << " ";
    cout << endl;
}

template <typename T>
class CompareIndicesByAnotherVectorValues
{
    std::vector<T> *_values;

public:
    CompareIndicesByAnotherVectorValues(std::vector<T> *values) : _values(values) {}

public:
    bool operator()(const int &a, const int &b) const { return (_values)[a] > (_values)[b]; }
};

int main()
{
    vector<vector<double>> x(3);
    x = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    for (int i{}; i < x.size(); i++)
    {
        for (int j{}; j < x.size(); j++)
        {
            cout << setprecision(2) << fixed << setw(8) << left << x[i][j];
        }
        cout << endl;
    }
    sort_indexes(x);
    return 0;
}