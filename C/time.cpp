#include <chrono>
#include <iostream>

using namespace std;
int main()
{
    double x = 1;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 1; i < 1000; i++)
        x += i;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[μs]" << std::endl;
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
    cout << x;
    return 0;
}
