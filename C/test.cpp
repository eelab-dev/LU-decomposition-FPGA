#include <iostream>
using namespace std;

void LUdecomposition(float **a, float **l, float **u, int n);

int main(void)
{
    int n, i, j;
    cout << "Enter size of square matrix : " << endl;
    cin >> n;

    int **a = new int *[n];
    for (i = 0; i < n; i++)
        a[i] = new int[n];

    cout << "Enter matrix values: " << endl;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            printf("a[%d][%d]=", i, j);
            cin >> a[i][j];
        }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            printf("a[%d][%d]=%d ", i, j, a[i][j]);

    for (i = 0; i < n; i++)
        delete[] a[i];
    delete[] a;

    return 0;
}

void LUdecomposition(float **a, float **l, float **u, int n)
{
    int i = 0, j = 0, k = 0;
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