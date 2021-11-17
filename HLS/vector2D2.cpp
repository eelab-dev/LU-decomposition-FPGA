// Unit Lower matrix

#include <hls_vector.h>

const int SIZE = 3;

void LUdecomposition(hls::vector<hls::vector<double,SIZE>,SIZE> &a, hls::vector<hls::vector<double,SIZE>,SIZE> &l, hls::vector<hls::vector<double,SIZE>,SIZE> &u)
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    int n = SIZE;
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

hls::vector<hls::vector<double,SIZE>,SIZE> lu()
{
//    int seed =2021;
//
//    std::mt19937 gen(seed);
//    std::uniform_int_distribution<int> distrib(1, 10);

    hls::vector<hls::vector<double,SIZE>,SIZE> matrix,l,u;

    for (int i{}; i < SIZE; i++)
        for (int j{}; j < SIZE; j++)
            matrix[i][j] = i+j;

    LUdecomposition(matrix, l, u);
    return l;

}
