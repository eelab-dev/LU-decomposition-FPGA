#include <iostream>
typedef struct
{
    int n;
    int Lip[5];
    double x[5];
} numeric;

int main()
{

    numeric *Numeric;
    Numeric->Lip[0] = 5;
    printf("Lip[0]=%d\n", Numeric->Lip[0]);

    for (int head = 0; head >= 0;)
    {
        head++;
        if (head > 5)
            head = -1;
        printf("head=%d\n", head);
    }
    return 0;
}