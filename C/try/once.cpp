#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>

const int TIME = 10; // Sample time

int main()
{
    std::string filename("lupivot2");
    std::string command;

    std::array<char, 128> buffer;
    std::string result;   // Stdout
    int returnCode, size; // size refers to the size of the matrix

    FILE *pipe;
    std::ofstream data;
    data.open("data_" + filename + ".csv");
    data << "Size,1,2,3,4,5,6,7,8,9,10,Average/ns\n";

    // regex expression for pattern to be searched
    std::regex regexp("Time difference = (\\d+)\\[ns\\]");

    // flag type for determining the matching behavior (in this case on string objects)
    std::smatch sm;

    int sizerange[] = {10, 100};
    unsigned long long sumtime = 0;

    std::array<unsigned long long, TIME> time;

    data << size << ",";
    for (int i = 0; i < TIME; i++)
    {
        command = "./" + filename;
        // std::cout << "Opening reading pipe" << std::endl;
        pipe = popen(command.c_str(), "r");
        if (!pipe)
        {
            std::cerr << "Couldn't start command." << std::endl;
            return 0;
        }
        while (fgets(buffer.data(), 128, pipe) != NULL)
        {
            // std::cout << "Reading..." << std::endl;
            result += buffer.data();
        }
        returnCode = pclose(pipe);
        // std::cout << result << std::endl;

        regex_search(result, sm, regexp); // Search for running time
        time[i] = stoull(sm[1]);
        data << time[i] << ",";
        sumtime += time[i];

        std::cout << i + 1 << ": " << time[i] << " ns" << std::endl;
        // std::cout << returnCode << std::endl;

        result.clear(); // clear the string
    }
    std::cout << "Size:" << size << ", Average Time: " << sumtime / double(TIME) << " ns\n";
    data << std::fixed << std::setprecision(0) << sumtime / double(TIME) << "\n";
    sumtime = 0;

    data.close();
    return 0;
}