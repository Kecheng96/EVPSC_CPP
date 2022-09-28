#ifndef TOOLBOX_H
#define TOOLBOX_H

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <iostream>
#include <Eigen/Dense>
#include <cmath>

#define SQR2 1.41421356237309
#define RSQ2 0.70710678118654744
#define RSQ3 0.57735026918962584
#define RSQ6 0.40824829046386304

using namespace std;
using namespace Eigen;

//
#define Intn 10 //number of integral points in Eshelby calculation
#define Mtr 1 //number of Multithread
//
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double,5,5> Matrix5d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<int, 6, 1> Vector6i;
typedef Matrix<double, 10, 1> Vector10d;
typedef Matrix<double,3,Intn*Intn> Integralpoint3;
typedef Matrix<double,6,Intn*Intn> Integralpoint6;
typedef Matrix<double,Intn*Intn,1> Integralpoint1;

//////////
struct Gausspoint
{
    //integral points and weights
    Integralpoint3 Gpalpha, Gpaww;
    Integralpoint6 Gpaa6, Gpaaww6; //coordinate and weigts in Fourier space 
    Integralpoint1 Gpww;
};

/////////////////////////
//choose the gausspoint set according to the length of ellisoid 
int Eshelby_case(Vector3d);
/////////////////////////

/////////////////////////
//Identity tensor
double Iij(int, int);
/////////////////////////

/////////////////////////
//transform into the deviatoric tensor
Matrix3d devia(Matrix3d);
Vector6d devia(Vector6d);
Matrix6d devia(Matrix6d);
/////////////////////////

/////////////////////////
//A(I,K)=B(I,J,K,L)*C(J,L) using voigt's notation
Vector6d Mult_voigt(Matrix6d, Vector6d);
/////////////////////////

/////////////////////////
//calculate the inverse of a forth order tensor A
void Inv_voigt(double A[3][3][3][3], double Ainv[3][3][3][3]);
//parameters: double A[3][3][3][3], double Ainv[3][3][3][3]
//Ainv = A^-1 as a output
/////////////////////////

/////////////////////////
//calculate the error between matrix A and B
double Errorcal(Matrix6d, Matrix6d);
double Errorcal(Vector6d, Vector6d);
double Errorcal(Matrix3d, Matrix3d);
double Errorcal(Matrix5d, Matrix5d);
/////////////////////////

/////////////////////////
//Sort the eigen values from largest to smallest 
// and change the order of the eigen vectors
void Eigsrt(Matrix3d &, Vector3d &);
/////////////////////////

/////////////////////////
//calculate C(i) = A(i)*B(i)
Vector6d mult_dot(Vector6d, Vector6i);
//calculate Cijkl = Aijmn * Bmnkl
void mult_4th(double A[3][3][3][3],double B[3][3][3][3],double C[3][3][3][3]);
//calculate Cij = Aijkl * Bkl
Matrix3d mult_4th(double A[3][3][3][3],Matrix3d B);
/////////////////////////

/////////////////////////
//THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION MATRIX
Matrix3d Euler_trans(Vector3d);
Vector3d Euler_trans(Matrix3d); 
/////////////////////////

/////////////////////////
//Chg_basis and its overload function
//with BASIS TENSOR Basis[3][3][6]
Matrix3d Chg_basis(VectorXd);
void Chg_basis(MatrixXd, double M3333[3][3][3][3]);
///
Vector6d Chg_basis6(Matrix3d);
Matrix6d Chg_basis6(Matrix6d);
Matrix6d Chg_basis6(double M3333[3][3][3][3]);
///
Vector5d Chg_basis5(Matrix3d);
Matrix5d Chg_basis5(Matrix6d);
Matrix5d Chg_basis5(double M3333[3][3][3][3]);
/////////////////////////

/////////////////////////
//change a B-basis tensor into voigt notation
Matrix6d Btovoigt(MatrixXd);
Matrix5d voigttoB5(Matrix6d);
Matrix6d voigttoB6(Matrix6d);
/////////////////////////

/////////////////////////
//Voigt and its overload function
//to realize the 3X3X3X3 to 6X6 and the reverse transformation
//and realize the 3X3 to 6X1 and the reverse transformation
//11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6 
Matrix3d voigt(Vector6d);  // 6X1 to 3X3
Vector6d voigt(Matrix3d);  // 3X3 to 6X1
Matrix6d voigt(double M3333[3][3][3][3]);    // 3X3X3X3 to 6X6
void voigt(Matrix6d, double M3333[3][3][3][3]); // 6X6 to 3X3X3X3
Vector10d voigt(Matrix4d);  // 4X4 to 10X1
Matrix4d voigt(Vector10d);  // 10X1 to 4X4
Vector6i voigt(Matrix3i); //int version

Matrix5d voigt6to5(Matrix6d);
Matrix6d voigt5to6(Matrix5d);
/////////////////////////

/////////////////////////
Matrix6d Bbasisadd(Matrix6d,Matrix5d);
Vector6d Bbasisadd(Vector6d,Vector5d);
/////////////////////////

/////////////////////////
//Gauss-Legendre points
//double x1, double x2, Vector x,w,int n
//calculate the points and weights (VectorXd x,w)
//according the Integral interval (x1, x2) and number of points (n)
void gau_leg(double, double, VectorXd &, VectorXd &, int);
/////////////////////////

/////////////////////////
//Rotate the stiffness or compliance in 6*6 matrix
Matrix6d rotate_C66(Matrix6d C66, Matrix3d M);
Matrix6d rot_stress(Matrix3d);
void rot_4th(double A[3][3][3][3], Matrix3d, double B[3][3][3][3]);
/////////////////////////

/////////////////////////
//Refer to Appendix D: Crystal rotation (Rodrigues) and misorientation
Matrix3d Rodrigues(Matrix3d);
/////////////////////////

int Jacobi(Matrix3d A,Vector3d &D,Matrix3d &V);

#endif