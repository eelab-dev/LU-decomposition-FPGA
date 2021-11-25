#include "SparseMatrix.cpp"

int main()
{
    SparseMatrix::SparseMatrix<int> matrix(3, 4);
    matrix.set(2, 1, 1);
    matrix.set(7, 1, 3);
    matrix.set(4, 2, 2);
    matrix.set(1, 3, 4);

    std::cout << matrix << std::endl
              << std::endl;

    SparseMatrix::SparseMatrix<int> matrixA(2, 3);
    SparseMatrix::SparseMatrix<int> matrixB(3, 4);

    SparseMatrix::SparseMatrix<int> product(2, 4); // will be of size 2×4
    product = matrixA.multiply(matrixB);           // method
    product = matrixA * matrixB;                   // operator

    std::cout << product << std::endl;
    return 0;
}