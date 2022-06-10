#include <iostream>

int main(void)
{
    int(*a)[8] = (int **)malloc(sizeof(int[4][8]));

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            a[i][j] = i + j;
        }
    }

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            printf("%d ", a[i][j]);
        }
        printf("\n");
    }

    return 0;
}