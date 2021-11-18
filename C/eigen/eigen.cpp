#include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/LU>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::PartialPivLU;
using Eigen::UpLoType;
using namespace std;

int main()
{
    MatrixXd P(3, 3);
    P << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    MatrixXd U(3, 1);
    U << 0, 1, 2;

    PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(P);
    MatrixXd J = lu.matrixLU().triangularView<UpLoType::Upper>();
    MatrixXd F = lu.matrixLU().triangularView<UpLoType::UnitLower>();
    MatrixXd perm = lu.permutationP().transpose();
    MatrixXd F1 = perm * F;

    MatrixXd Jlamda = F.lu().solve(U);
    MatrixXd l = J.lu().solve(Jlamda);

    cout << perm << endl;
    cout << endl;
    cout << F1 << endl;
    cout << J << endl;
    cout << F << endl;

    cout << typeid(F(1, 0)).name() << endl;
    return 0;
}