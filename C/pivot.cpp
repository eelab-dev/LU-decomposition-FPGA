#include <cmath>
#include <iomanip>
#include <iostream>

double myround(double x, int n)
{
    if (n < 0)
        return x;
    else
        return round(x * pow(10, n)) / pow(10, n);
}

int main()
{
    // std::cout << myround(10.234, 0) << std::endl;
    double x[2][3] = {{0.003, 59.14, 59.17}, {5.291, -6.13, 46.78}};
    double y[2][3] = {{5.291, -6.13, 46.78}, {0.003, 59.14, 59.17}};
    double x1[2][3] = {{0.003, 59.14, 59.17}, {5.291, -6.13, 46.78}};
    double y1[2][3] = {{5.291, -6.13, 46.78}, {0.003, 59.14, 59.17}};
    // double x[2][3] = {{1, 1, 2}, {10000, 1, 10000}};
    // double y[2][3] = {{10000, 1, 10000}, {1, 1, 2}};
    // double x1[2][3] = {{1, 1, 2}, {10000, 1, 10000}};
    // double y1[2][3] = {{10000, 1, 10000}, {1, 1, 2}};
    double xd[2], yd[2];

    for (int i = 0; i < 3; i++)
    {
        x1[1][i] = x1[1][i] - x[0][i] * x[1][0] / x[0][0];

        y1[1][i] = y1[1][i] - y[0][i] * y[1][0] / y[0][0];
    }

    for (int i = 0; i < 2; i++)
    {
        std::cout << std::setprecision(20) << std::fixed << x1[i][0] << " " << x1[i][1] << " " << x1[i][2] << " " << std::endl;
    }

    xd[1] = x1[1][2] / x1[1][1];
    xd[0] = (x[0][2] - x[0][1] * xd[1]) / x[0][0];

    yd[1] = y1[1][2] / y1[1][1];
    yd[0] = (y[0][2] - y[0][1] * yd[1]) / y[0][0];

    std::cout << std::setprecision(15) << std::fixed << xd[0] << ", " << xd[1] << std::endl;
    std::cout << std::setprecision(15) << std::fixed << yd[0] << ", " << yd[1] << std::endl;

    return 0;
}