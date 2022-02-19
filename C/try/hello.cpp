#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main()
{
    int n, i;
    cout << "Enter size of square matrix : " << endl;
    cin >> n;
    int *a;
    a = new int[n];
    for (i = 0; i <= n * 2; i++)
    {
        printf("a[%d]=", i);
        cin >> a[i];
    }
    for (i = 0; i <= n * 2; i++)
        cout << "a[" << i << "]=" << a[i] << " ";
    delete a;
}