#include <iostream>

int main(int argc, char *argv[])
{
    int size = 10;
    if (argc > 2)
        if (std::string{"--size"}.compare(argv[1]) == 0 || std::string{"-t"}.compare(argv[1]) == 0)
            size = std::stoi(argv[2]);
    std::cout << "argc:" << argc << " size:" << size << std::endl;
    return 0;
}