#include <array>
#include <iostream>
#include <regex>
#include <string>

using namespace std;
int main()
{
    array<float, 2> time;
    //string to be searched
    string str = "Time difference = 4400[ns]\nTime difference = 5600[ns]";

    // regex expression for pattern to be searched
    regex regexp("=\\s(\\d+)\\[ns\\]");

    // flag type for determining the matching behavior (in this case on string objects)
    smatch sm;

    // regex_search that searches pattern regexp in the string mystr
    string::const_iterator searchStart(str.cbegin());
    int i = 0;
    while (regex_search(searchStart, str.cend(), sm, regexp))
    {
        time[i++] = stoi(sm[1]);
        cout << sm[1] << " ";
        searchStart = sm.suffix().first; // Locate to the last character of the matched string
    }
    cout << endl;
    double sum = 0;
    for (i = 0; i < time.size(); i++)
        sum += time[i];
    cout << "Average Time:" << sum / time.size() << endl;
    return 0;
}