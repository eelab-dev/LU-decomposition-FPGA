#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

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
                { return v[i1][0] > v[i2][0]; });

    return idx;
}

int main()
{
    int seed = 2021;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 100);

    vector<double> x(10);
    double *p[10];

    for (int i{}; i < x.size(); i++)
    {
        x[i] = distrib(gen);
        cout << setw(4) << left << x[i];
    }
    cout << endl;

    int j{};
    for (int i : sort_indexes(x))
    {
        p[j++] = &x[i];
        cout << i << "," << setw(4) << left << x[i];
    }
    cout << endl;

    for (int i; i < x.size(); i++)
        cout << setw(4) << left << *p[i];
    cout << endl;

    return 0;
}