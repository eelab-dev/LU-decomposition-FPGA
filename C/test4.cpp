#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;
int main()
{
    vector<vector<double>> x(3);
    x = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    for (int i{}; i < x.size(); i++)
    {
        for (int j{}; j < x.size(); j++)
        {
            cout << setprecision(2) << fixed << setw(8) << left << x[i][j];
        }
        cout << endl;
    }
    return 0;
}