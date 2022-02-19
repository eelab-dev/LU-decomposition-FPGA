#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;
template <typename T>
vector<size_t> sort_indexes(const vector<vector<T>> &v)
{

    // initialize original index locations
    vector<size_t> idx(v.size());
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

const int SIZE = 10000;

int main()
{
    int seed = 2021;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 10);

    vector<vector<double>> matrix(SIZE, vector<double>(SIZE));

    for (int i{}; i < SIZE; i++)
        for (int j{}; j < SIZE; j++)
            matrix[i][j] = distrib(gen);

    for (int i{}; i < 10; i++)
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        sort_indexes(matrix);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << endl;
    }

    return 0;
}