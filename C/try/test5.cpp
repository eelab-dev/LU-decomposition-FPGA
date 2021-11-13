#include <filesystem>
#include <iostream>
int main()
{
    std::string outfile = "data_lupivot_cmd_O3.csv";
    try
    {
        if (std::filesystem::remove(outfile))
            std::cout << "file " << outfile << " deleted.\n";
        else
            std::cout << "file " << outfile << " not found.\n";
    }
    catch (const std::filesystem::filesystem_error &err)
    {
        std::cout << "filesystem error: " << err.what() << '\n';
    }
    return 0;
}