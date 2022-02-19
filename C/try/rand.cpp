#include <iostream>
#include <random>
int main()
{
    std::random_device rd;
    int seed = 2021;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 6);

    std::cout << sizeof(gen) << std::endl;

    int x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    for (unsigned long j = 0; j < 10; ++j)
    {
        std::cout << distrib(gen) << ' ';
    }

    std::cout << '\n';
    return 0;
}