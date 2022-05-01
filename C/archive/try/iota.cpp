#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

int main(void)
{
    std::vector<int> x(10);
    std::iota(x.begin(), x.end(), 0);
    std::transform(x.begin(), x.end(), x.begin(), [&](double x)
                   { return (x + 1.0); });

    for (int i = 0; i < x.size(); i++)
    {
        std::cout << "x[" << i << "]=" << x[i] << std::endl;
    }

    std::vector<std::string> filename = {"host.mtx", "rajat14.mtx"};
    std::string prefix = "../../Matrix_Sample/Bench/";

    std::transform(filename.begin(), filename.end(), filename.begin(), [&](std::string x)
                   { return (prefix + x); });

    for (int i = 0; i < filename.size(); i++)
        std::cout << filename[i] << std::endl;

    return 0;
}