#include <chrono>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <random>

using Eigen::FullPivLU;
using Eigen::MatrixXd;
using Eigen::PartialPivLU;
using Eigen::UpLoType;
using namespace std;

const int SIZE = 4;

int main()
{
    int seed = 2021;

    std::random_device rd;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<unsigned> distrib(1, 10);

    MatrixXd P(SIZE, SIZE);
    for (int i{}; i < SIZE; i++)
        for (int j{}; j < SIZE; j++)
            P(i, j) = distrib(gen);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(P);
    MatrixXd U = lu.matrixLU().triangularView<UpLoType::Upper>();
    MatrixXd L = lu.matrixLU().triangularView<UpLoType::UnitLower>();
    // MatrixXd perm = lu.permutationP().transpose();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    cout << P << endl;
    cout << "Upper=" << endl
         << U << endl;
    cout << "Lower=" << endl
         << L << endl;
    cout << "Reconstruct=" << endl
         << L * U << endl;
    // cout << "Perm=" << endl << perm << endl;

    cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << endl;

    return 0;
}