#include <iostream>
#include "../C/myKLU/include/mmio.h"

int main()
{
    std::string homeDir = getenv("HOME");
    std::string bmatrix = homeDir + "/beng-project/Matrix_Sample/host_b.mtx";

    std::vector<double> b;
    int nrhs;
    read_bmatrix(bmatrix, b, &nrhs);

    return 0;
}