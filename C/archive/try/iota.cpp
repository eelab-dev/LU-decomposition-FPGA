#include <iostream>
#include <vector>
#include <numeric>

int main(void)
{
    int i = 10;
    std::vector<int> x(10);
    std::iota(x.begin(), x.end(), 0);

    for (int i = 0; i < 10; i++)
    {
        std::cout << "x[" << i << "]=" << x[i] << std::endl;
    }
    return 0;
}