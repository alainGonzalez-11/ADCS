#include <iostream>
#include <D:\Documentos\GitHub\ADCS\C\Eigen\Dense>
#include <math.h>

using namespace Eigen;

Vector3d ESOQ2 ( VectorXd a_i, MatrixXd v_i, MatrixXd s_i)
{
    Matrix3d B = a_i.cwiseProduct(s_i) * v_i.transpose();
    Matrix3d S = B + B.transpose();
    double sigma = B.trace();
    Vector3d Z;
    Z << B(1,2) - B(2,1),
        B(2,0) - B(0,2),
        B(0,1) - B(1,0);
    Matrix3d B_k = S -  MatrixXd::Identity(3, 3) * sigma;
    Matrix4d K;
    K << sigma, Z.transpose(),
        Z, B_k;
    double A = a_i.sum();
    double error = 1;
    while (error > 1e-10){
        double A_old = A;
        Matrix4d I;
        I = MatrixXd::Identity(4, 4);
        double r = (K - A * I).determinant();
        double r1 = (K - A * I).determinant() * (((K - A * I).inverse()) * (- I)).trace();
        A = A - r/r1;
        error = abs(A - A_old);
    }
    Vector4d q = (((A + sigma) * MatrixXd::Identity(3, 3) - S).inverse())* Z;
    return q;
}

int main()
{
    MatrixXd v_i;
    v_i.resize(3,2);
    v_i << 1,     0,
           0,     0,
           0,     1;
    MatrixXd s_i;
    v_i.resize(3,2);
    v_i << 0.0100,    0.0050,
           0.9999,   -0.0050,
           0.0050,    0.9999;
    VectorXd weight;
    weight.resize(2);
    weight << 0.5, 0.5;
    ESOQ2(weight, v_i, s_i);


}
