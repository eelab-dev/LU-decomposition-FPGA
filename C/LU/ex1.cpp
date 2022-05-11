// https://people.math.sc.edu/Burkardt/cpp_src/umfpack/umfpack_simple.cpp

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <vector>

#include <suitesparse/umfpack.h>

int main();

//****************************************************************************80

int main()

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for UMFPACK_SIMPLE.
//
//  Discussion:
//
//    This program uses UMFPACK to solve the 5x5 linear system A*X=B:
//
//        2  3  0  0  0        1.0         8.0
//        3  0  4  0  6        2.0        45.0
//    A = 0 -1 -3  2  0    X = 3.0    B = -3.0
//        0  0  1  0  0        4.0         3.0
//        0  4  2  0  1        5.0        10.0
//
//    The matrix contains 12 nonzero values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2014
//
//  Reference:
//
//    Timothy Davis,
//    UMFPACK User Guide,
//    Version 5.6.2, 25 April 2013
//    http://suitesparse.com
//
{

    const int N = 5;
    const int NCC = 12;
    std::vector<int> Ai = {
        0, 1,
        0, 2, 4,
        1, 2, 3, 4,
        2,
        1, 4};
    std::vector<int> Ap = {0, 2, 5, 9, 10, 12};
    std::vector<double> Ax = {
        2.0, 3.0,
        3.0, -1.0, 4.0,
        4.0, -3.0, 1.0, 2.0,
        2.0,
        6.0, 1.0};
    std::vector<double> b = {8.0, 45.0, -3.0, 3.0, 19.0};
    int i;
    int n = 5;
    double *null = (double *)NULL;
    void *Numeric;
    int status;
    void *Symbolic;
    // double x[N];
    std::vector<double> x(N);

    std::cout << "\n";
    std::cout << "UMFPACK_SIMPLE:\n";
    std::cout << "  C++ version\n";
    std::cout << "  Use UMFPACK to solve the sparse linear system A*x=b.\n";
    //
    //  Carry out the symbolic factorization.
    //
    status = umfpack_di_symbolic(n, n, Ap.data(), Ai.data(), Ax.data(), &Symbolic, null, null);
    //
    //  Use the symbolic factorization to carry out the numeric factorization.
    //
    status = umfpack_di_numeric(Ap.data(), Ai.data(), Ax.data(), Symbolic, &Numeric, null, null);
    //
    //  Free the memory associated with the symbolic factorization.
    //
    umfpack_di_free_symbolic(&Symbolic);
    //
    //  Solve the linear system.
    //
    status = umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Ax.data(), x.data(), b.data(), Numeric, null, null);
    //
    //  Free the memory associated with the numeric factorization.
    //
    umfpack_di_free_numeric(&Numeric);
    //
    //  Print the solution.
    //
    std::cout << "\n";
    std::cout << "  Computed solution:\n";
    std::cout << "\n";
    for (i = 0; i < n; i++)
    {
        std::cout << "  x[" << i << "] = " << x[i] << "\n";
    }
    //
    //  Terminate.
    //
    std::cout << "\n";
    std::cout << "UMFPACK_SIMPLE:\n";
    std::cout << "  Normal end of execution.\n";
    std::cout << "\n";

    return 0;
}