#include <iostream>
using namespace std;

class Matrix
{
    float *array;
    int m_width;

public:
    Matrix(int w, int h) : m_width(w), array(new float[w * h]) {}
    ~Matrix() { delete[] array; }
    float at(int y, int x) const { return array[index(y, x)]; }
    void set(int t, int y, int x) { array[index(y, x)] = t; }
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
    int n, i, j, temp;
    cout << "Enter size of square matrix : " << endl;
    cin >> n;

    Matrix a(n, n);
    Matrix l(n, n);
    Matrix u(n, n);

    // a.set(5, 1, 1);
    // cout << a.at(0, 0) << " " << a.at(1, 1);

    cout << "Enter matrix values: " << endl;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            printf("a[%d][%d]=", i, j);
            cin >> temp;
            a.set(temp, i, j);
        }

    LUdecomposition(a, l, u);

    cout << "A" << endl;
    a.display();
    cout << "L" << endl;
    l.display();
    cout << "U" << endl;
    u.display();

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