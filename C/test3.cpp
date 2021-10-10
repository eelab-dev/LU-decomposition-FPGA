#include <chrono>
#include <iostream>
#include <random>
using namespace std;

class Matrix
{
    float *array;
    int m_width;

public:
    Matrix(int w, int h) : m_width(w), array(new float[w * h]) {}
    ~Matrix() { delete[] array; }
    float at(int y, int x) const { return array[index(y, x)]; }
    void set(float t, int y, int x) { array[index(y, x)] = t; }
    int getwidth() { return m_width; }
    void display()
    {
        for (int i = 0; i < m_width; i++)
        {
            for (int j = 0; j < m_width; j++)
                cout << array[index(i, j)] << " ";
            cout << endl;
        }
    }

protected:
    int index(int y, int x) const { return x + m_width * y; }
};

void LUdecomposition(Matrix a, Matrix l, Matrix u);

int main(void)
{
    int n = 100;
    int seed = 2021;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 10);

    // cout << "Enter size of square matrix : " << endl;
    // cin >> n;

    Matrix a(n, n);
    Matrix l(n, n);
    Matrix u(n, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a.set(distrib(gen), i, j);

    cout << "A" << endl;
    a.display();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    LUdecomposition(a, l, u);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    cout << "L" << endl;
    // l.display();
    cout << "U" << endl;
    // u.display();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
    return 0;
}

void LUdecomposition(Matrix a, Matrix l, Matrix u)
{
    int n = a.getwidth();
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                l.set(0, j, i);
            else
            {
                l.set(a.at(j, i), j, i);
                for (k = 0; k < i; k++)
                    l.set(l.at(j, i) - l.at(j, k) * u.at(k, i), j, i);
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                u.set(0, i, j);
            else if (j == i)
                u.set(1, i, j);
            else
            {
                u.set(a.at(i, j) / l.at(i, i), i, j);
                for (k = 0; k < i; k++)
                    u.set(u.at(i, j) - l.at(i, k) * u.at(k, j) / l.at(i, i), i, j);
            }
        }
    }
}