#include "Toolbox.h"
using namespace std;
using namespace Eigen;

/////////////////////////
//Chg_basis and its overload function
//with BASIS TENSOR Basis[3][3][6]
const double Basis[3][3][6] = {
     //1,1
     -RSQ2, -RSQ6, 0, 0, 0, RSQ3,
     //1,2
     0, 0, 0, 0, RSQ2, 0,
     //1,3
     0, 0, 0, RSQ2, 0, 0,
     //2,1
     0, 0, 0, 0, RSQ2, 0,
     //2,2
     RSQ2, -RSQ6, 0, 0, 0, RSQ3,
     //2,3 
     0, 0, RSQ2, 0, 0, 0,
     //3,1
     0, 0, 0, RSQ2, 0, 0,
     //3,2
     0, 0, RSQ2, 0, 0, 0,
     //3,3
     0, 2.0*RSQ6, 0, 0, 0, RSQ3
     };

     
/////////////////////////

/////////////////////////
//Voigt and its overload function
//to realize the 3X3X3X3 to 6X6 and the reverse transformation
//and realize the 3X3 to 6X1 and the reverse transformation
const int IJV[6][2]= {0,0,1,1,2,2,1,2,0,2,0,1};
/////////////////////////

///
int Eshelby_case(Vector3d Axis){
    double RATIO1=Axis(1)/Axis(2);
    double RATIO2=Axis(0)/Axis(2);
    Matrix<double, 11, 1> DTE;
    DTE << 0.0,
    -0.7*RATIO1+7,
    -RATIO1+17,
    -RATIO1+23,
    -RATIO1+26,
    -RATIO1+29.3,
    -RATIO1+32,
    -RATIO1+34.85,
    -RATIO1+37,
    -RATIO1+41.9,
    -RATIO1+44.5;
    int case_c = 11;
    for(int i = 1; i <= 10; i++)
        if(RATIO2 >= DTE(i-1)&& RATIO2 < DTE(i)) case_c = i;
    return case_c-1;
}


/// function definition of gau_leg
void gau_leg(double x1, double x2, VectorXd &x, VectorXd &w, int n)
{
	double eps = 1e-7;

	double m = n / 2;
	double xm = 0.5*(x1 + x2);
	double xl = 0.5*(x2 - x1);
	double xn = n;

	for (int i = 0; i < m; i++)
	{
		double xi = i + 1;
		double z = cos(M_PI*(xi - 0.25) / (xn + 0.5));
		double p1, p2, p3, xj, pp, z1;

		z1 = z + 1;

		while (abs(z1 - z) > eps)
		{
		p1 = 1;
		p2 = 0;
		for (int j = 0; j < n; j++)
		{
			xj = j + 1;
			p3 = p2;
			p2 = p1;
			p1 = ((2 * j + 1)*z*p2 - (xj - 1)*p3) / xj;
		}

		pp = n*(z*p1 - p2) / (z*z - 1);
		z1 = z;
		z = z1 - p1 / pp;
		}

		x(i) = xm - xl*z;
		x(n - 1 - i) = xm + xl*z;
		w(i) = 2 * xl / ((1 - z*z)*pp*pp);
		w(n - 1 - i) = w(i);
	}
}

/// function definition of voigt
Matrix3d voigt(Vector6d Vin)
{
    Matrix3d Mout;
    int I1,I2;
    for(int i = 0; i < 6; i++)
    {
        I1 = IJV[i][0];
        I2 = IJV[i][1];
        Mout(I1,I2) = Vin(i);
        Mout(I2,I1) = Vin(i);
    }
    return Mout;
}

Vector6d voigt(Matrix3d Min)
{
    Vector6d Vout(6);
    int I1,I2;
    for(int i = 0; i < 6; i++)
        {
        I1 = IJV[i][0];
        I2 = IJV[i][1];
        Vout(i) = Min(I1,I2);
        }
    return Vout;
}

Vector6i voigt(Matrix3i Min)
{
    Vector6i Vout(6);
    int I1,I2;
    for(int i = 0; i < 6; i++)
        {
        I1 = IJV[i][0];
        I2 = IJV[i][1];
        Vout(i) = Min(I1,I2);
        }
    return Vout;
}

Matrix6d voigt(double M3333[3][3][3][3])
{
    Matrix6d Mout(6,6);
    int I1,I2,J1,J2;
    for(int i = 0; i < 6; i++)
    {
        I1 = IJV[i][0];
        I2 = IJV[i][1];
        for(int j = 0; j < 6; j++)
        {
            J1 = IJV[j][0];
            J2 = IJV[j][1];
            Mout(i,j) = M3333[I1][I2][J1][J2];
        }
    }
    return Mout;
}

void voigt(Matrix6d Min, double M3333[3][3][3][3])
{
    int I1,I2,J1,J2;
    for(int i = 0; i < 6; i++)
    {
        I1 = IJV[i][0];
        I2 = IJV[i][1];
        for(int j = 0; j < 6; j++)
        {
            J1 = IJV[j][0];
            J2 = IJV[j][1];
            M3333[I1][I2][J1][J2] = Min(i,j);
            M3333[I2][I1][J1][J2] = Min(i,j);
            M3333[I1][I2][J2][J1] = Min(i,j);
            M3333[I2][I1][J2][J1] = Min(i,j);
        }
    }
}

//11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6 \
  14-->7, 24-->8, 34-->9, 44-->10
Vector10d voigt(Matrix4d Min)
{
    Vector10d Vout;
    Vout(0) = Min(0,0);
    Vout(1) = Min(1,1); 
    Vout(2) = Min(2,2); 
    Vout(3) = Min(1,2); 
    Vout(4) = Min(2,0); 
    Vout(5) = Min(0,1); 
    Vout(6) = Min(0,3); 
    Vout(7) = Min(1,3); 
    Vout(8) = Min(2,3); 
    Vout(9) = Min(3,3);  
    return Vout;
}

Matrix4d voigt(Vector10d Vin)
{
    Matrix4d Mout;
    Mout(0,0) = Vin(0);
    Mout(1,1) = Vin(1);
    Mout(2,2) = Vin(2);
    Mout(1,2) = Vin(3); Mout(2,1) = Vin(3);
    Mout(2,0) = Vin(4); Mout(0,2) = Vin(4);
    Mout(0,1) = Vin(5); Mout(1,0) = Vin(5);
    Mout(0,3) = Vin(6); Mout(3,0) = Vin(6);
    Mout(1,3) = Vin(7); Mout(3,1) = Vin(7);
    Mout(2,3) = Vin(8); Mout(3,2) = Vin(8);
    Mout(3,3) = Vin(9);
    return Mout;
}

Matrix5d voigt6to5(Matrix6d M66)
{
    Matrix<double,5,5> M55;
	for (int i = 0; i < 2; i++)
	{
		M55(0,i) = M66(0,i);
		M55(1,i) = M66(1,i);
		M55(2,i) = M66(3,i);
		M55(3,i) = M66(4,i);
		M55(4,i) = M66(5,i);
	}
	for (int i = 2; i < 5; i++)
	{
		M55(0,i) = M66(0,i+1);
		M55(1,i) = M66(1,i+1);
		M55(2,i) = M66(3,i+1);
		M55(3,i) = M66(4,i+1);
		M55(4,i) = M66(5,i+1);
	}

	for (int i = 0; i < 2; i++)
	{
		M55(0,i) = M55(0,i) - M66(0,2);
		M55(1,i) = M55(1,i) - M66(1,2);
		M55(2,i) = M55(2,i) - M66(3,2);
		M55(3,i) = M55(3,i) - M66(4,2);
		M55(4,i) = M55(4,i) - M66(5,2);
	}
    return M55;
}

Matrix6d voigt5to6(Matrix5d M55)
{
    Matrix6d M66;
    for (int i = 0; i < 2; i++)
	{
		M66(0,i) = M55(0,i);
		M66(1,i) = M55(1,i);
		M66(3,i) = M55(2,i);
		M66(4,i) = M55(3,i);
		M66(5,i) = M55(4,i);
    }
	for (int i = 3; i < 6; i++)
	{
		M66(0,i) = M55(0,i-1);
		M66(1,i) = M55(1,i-1);
		M66(3,i) = M55(2,i-1);
		M66(4,i) = M55(3,i-1);
		M66(5,i) = M55(4,i-1);

		M66(2,i) = - M66(0,i) - M66(1,i);
		M66(i,2) = M66(2,i);
	}

	M66(2,2) = 0;
	for (int i = 3; i < 6; i++)
	{
		M66(i,0) = M66(i,0) + M66(i,2);
		M66(i,1) = M66(i,1) + M66(i,2);
	}


	for (int i = 0; i < 2; i++)
	{
		M66(2,i) = -M66(0,i) - M66(1,i);
		M66(i,2) = M66(2,i);
	}

	for (int i = 0; i < 2; i++)
	{
		M66(0,i) = M66(0,i) + M66(0,2);
		M66(1,i) = M66(1,i) + M66(1,2);
	}
    return M66;
}

/// function definition of Chg_basis
Matrix3d Chg_basis(VectorXd Vin)
{
    Matrix3d Mout = Matrix3d::Zero();
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int n = 0; n < Vin.size(); n++)
                Mout(i,j) += Vin(n) * Basis[i][j][n];
    return Mout;
}

void Chg_basis(MatrixXd Min, double M3333[3][3][3][3])
{
    int KDIM = Min.cols();
    double dummy = 0; 
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
    {
        dummy = 0;
        for(int n = 0; n < KDIM; n++)
            for(int m = 0; m < KDIM; m++)
                dummy += Min(n,m)*Basis[i][j][n]*Basis[k][l][m];
        M3333[i][j][k][l] = dummy;
    }
}

Vector6d Chg_basis6(Matrix3d Min)
{
    Vector6d Vout = Vector6d::Zero();
    for(int n = 0; n < 6; n++)
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
    Vout(n) += Min(i,j)*Basis[i][j][n];
    return Vout;
}

Matrix6d Chg_basis6(Matrix6d Min)
{
    double M3333[3][3][3][3] ={0};
    voigt(Min,M3333);
    return Chg_basis6(M3333);
}

Matrix6d Chg_basis6(double M3333[3][3][3][3])
{
    Matrix6d Mout = Matrix6d::Zero();
    for(int n = 0; n < 6; n++)
        for(int m = 0; m < 6; m++)
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                    for(int k = 0; k < 3; k++)
                        for(int l = 0; l < 3; l++)
    Mout(n,m) += M3333[i][j][k][l]*Basis[i][j][n]*Basis[k][l][m];

    return Mout;
}

Vector5d Chg_basis5(Matrix3d Min)
{
    Vector5d Vout = Vector5d::Zero();
    for(int n = 0; n < 5; n++)
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
    Vout(n) += Min(i,j)*Basis[i][j][n];
    return Vout;
}

Matrix5d Chg_basis5(Matrix6d Min)
{
    double M3333[3][3][3][3] ={0};
    voigt(Min,M3333);
    return Chg_basis5(M3333);
}

Matrix5d Chg_basis5(double M3333[3][3][3][3])
{
    Matrix5d Mout = Matrix5d::Zero();
    for(int n = 0; n < 5; n++)
        for(int m = 0; m < 5; m++)
            for(int i = 0; i < 3; i++)
                for(int j = 0; j < 3; j++)
                    for(int k = 0; k < 3; k++)
                        for(int l = 0; l < 3; l++)
    Mout(n,m) += M3333[i][j][k][l]*Basis[i][j][n]*Basis[k][l][m];

    return Mout;
}

/// function definition of Euler
Matrix3d Euler_trans(Vector3d Vin)
{
    Matrix3d Mout;
    double PH,TH,TM;
    Vin = Vin / 180 * M_PI; // in radians
    PH = Vin(0);
    TH = Vin(1);
    TM = Vin(2);

    double SPH,STH,STM;
    double CPH,CTH,CTM;
    SPH=sin(PH); CPH=cos(PH);
    STH=sin(TH); CTH=cos(TH);
    STM=sin(TM); CTM=cos(TM);
    Mout(0,0)=CTM*CPH-SPH*STM*CTH;
    Mout(1,0)=-STM*CPH-SPH*CTM*CTH;
    Mout(2,0)=SPH*STH;
    Mout(0,1)=CTM*SPH+CPH*STM*CTH;
    Mout(1,1)=-SPH*STM+CPH*CTM*CTH;
    Mout(2,1)=-STH*CPH;
    Mout(0,2)=STH*STM;
    Mout(1,2)=CTM*STH;
    Mout(2,2)=CTH;

    return Mout;
}

Vector3d Euler_trans(Matrix3d Min)
{
    Vector3d Vout;
    double PH,TH,TM;
    double STH;

    TH = acos(Min(2,2));
    if(abs(Min(2,2))>= 0.999)
    {
        TM = 0;
        PH = atan2(Min(0,1),Min(0,0));
    }
    else
    {
        STH = sin(TH);
        TM = atan2(Min(0,2)/STH,Min(1,2)/STH);
        PH = atan2(Min(2,0)/STH,-Min(2,1)/STH); 
    }
    Vout << PH,TH,TM;
    Vout = Vout * 180 / M_PI; // in degree
    return Vout;
}

/// function definition of Identity tensor
double Iij(int i, int j)
{
    if(i == j) return 1;
    return 0;
}

// function definition of Eigsrt
void Eigsrt(Matrix3d &V, Vector3d &D)
{
    int n = 3;
    double P,t;
    int K;
    for(int i = 0; i < n-1; i++)
    {
        K = i;
        P = D(i);
        for(int j = i+1; j < n; j++)
            if(D(j) >= P) 
            {   K = j; P = D(j);}

        if(K != i)
            {
                D(K) = D(i);
                D(i) = P;
                for(int m = 0; m < n; m++)
                {
                    t = V(m,i);
                    V(m,i) = V(m,K);
                    V(m,K) = t;
                }
            }
    }

    // int count = 0;
    // for(int j = 0; j < 3; j++)
    // {
    //     count = 0;
    //     for(int i = 0; i < 3; i++)
    //         if( V(i,j) < 0) count+=1;
    //     if(count > 2) V.col(j) = -V.col(j);
    // }
}

//fuction definition of Mult_voigt
Vector6d Mult_voigt(Matrix6d B, Vector6d C)
{
    Vector6d A;
    A(0)=B(0,0)*C(0)+B(5,5)*C(1)+B(4,4)*C(2)\
         +2*(B(4,5)*C(3)+B(0,4)*C(4)+B(0,5)*C(5));
    
    A(1)=B(5,5)*C(0)+B(1,1)*C(1)+B(3,3)*C(2)\
         +2*(B(1,3)*C(3)+B(3,5)*C(4)+B(1,5)*C(5));
    
    A(2)=B(4,4)*C(0)+B(3,3)*C(1)+B(2,2)*C(2)\
         +2*(B(2,3)*C(3)+B(2,4)*C(4)+B(3,4)*C(5));
    
    A(3)=B(4,5)*C(0)+B(1,3)*C(1)+B(2,3)*C(2)\
         +(B(1,2)+B(3,3))*C(3)\
         +(B(2,5)+B(3,4))*C(4)\
         +(B(3,5)+B(1,4))*C(5);
    
    A(4)=B(0,4)*C(0)+B(3,5)*C(1)+B(2,4)*C(2)\
         +(B(2,5)+B(3,4))*C(3)\
         +(B(0,2)+B(4,4))*C(4)\
         +(B(0,3)+B(4,5))*C(5);
    
    A(5)=B(0,5)*C(0)+B(1,5)*C(1)+B(3,4)*C(2)\
        +(B(3,5)+B(1,4))*C(3)\
        +(B(0,3)+B(4,5))*C(4)\
        +(B(0,1)+B(5,5))*C(5);
    return A;
}

//fuction definition of Inv_voigt
void Inv_voigt(double A[3][3][3][3], double Ainv[3][3][3][3])
{
    Matrix6d Ainv66;
    Ainv66 = voigt(A).inverse();
    voigt(Ainv66, Ainv);
}

//fuction definition of Errorcal
double Errorcal(Matrix6d A, Matrix6d B)
{
    double err;
    Matrix6d C1,C2;
    C1 = A - B;
    C2 = 0.5*(A + B);
    err = C1.norm()/C2.norm();
    return err;
}

double Errorcal(Vector6d A, Vector6d B)
{
    double err;
    double A1 = A.norm();
    double B1 = B.norm();
    double Tmax = max(A1, B1);
    double Tmin = min(A1, B1);
    if(Tmax==0)
    {err = 0;}
    else if(Tmin/Tmax < 1e-3)
    {err = Tmax - Tmin;}
    else
    {
        Vector6d C1,C2;
        C1 = A - B;
        C2 = 0.5*(A + B);
        err = C1.norm()/C2.norm();
    }
    return err;
}

double Errorcal(Matrix3d A, Matrix3d B)
{
    double err;
    Matrix3d C1,C2;
    C1 = A - B;
    C2 = 0.5*(A + B);
    err = C1.norm()/C2.norm();
    return err;
}

double Errorcal(Matrix5d A, Matrix5d B)
{
    double err;
    Matrix5d C1,C2;
    C1 = A - B;
    C2 = 0.5*(A + B);
    err = C1.norm()/C2.norm();
    return err;
}

//function definition of mult_dot
Vector6d mult_dot(Vector6d Vin, Vector6i Iv)
{
    Vector6d Vout;
    int n = 6;
    for(int i = 0; i < n; i++)
        Vout(i) = Vin(i)*Iv(i);
    return Vout;
}

void mult_4th(double A[3][3][3][3],double B[3][3][3][3],double C[3][3][3][3])
{
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
                    for(int m = 0; m < 3; m++)
                        for(int n = 0; n < 3; n++)
    C[i][j][k][l] = A[i][j][m][n] * B[m][n][k][l];        
}

Matrix3d mult_4th(double A[3][3][3][3],Matrix3d B)
{
    Matrix3d C = Matrix3d::Zero();
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
    C(i,j) = A[i][j][k][l] * B(k,l);
    return C;
}

//function definition of devia
Matrix3d devia(Matrix3d X)
{
    Matrix3d Mout = X;
    //transform into the deviatoric tensor
    double pr = (X(0,0) + X(1,1) + X(2,2))/3;
    Mout(0,0) = X(0,0) - pr;
    Mout(1,1) = X(1,1) - pr;
    Mout(2,2) = X(2,2) - pr;
    return Mout;
}

Vector6d devia(Vector6d X)
{
    Vector6d Vout = X;
    //transform into the deviatoric tensor
    double pr = (X(0) + X(1) + X(2))/3;
    Vout(0) = X(0) - pr;
    Vout(1) = X(1) - pr;
    Vout(2) = X(2) - pr;
    return Vout;
}

Matrix6d devia(Matrix6d Min)
{
    Matrix6d deviaM;
    deviaM<< 2/3, -1/3, -1/3, 0, 0, 0,
    -1/3, 2/3, -1/3, 0, 0, 0,
    -1/3, -1/3, 2/3, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1;
    return Min*deviaM;
}

//definition of rotate 6*6 tensor
Matrix6d rotate_C66(Matrix6d C66, Matrix3d E)
{
    Matrix6d E66, E66T;
    E66 = rot_stress(E);
    E66T = E66.transpose();
    return E66*C66*E66T;
}

Matrix6d rot_stress(Matrix3d M)
{
	Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = 2*xy*xz;
    M66(0,4) = 2*xx*xz;
    M66(0,5) = 2*xx*xy;

    M66(1,3) = 2*yy*yz;
    M66(1,4) = 2*yx*yz;
    M66(1,5) = 2*yx*yy;

    M66(2,3) = 2*zy*zz;
    M66(2,4) = 2*zx*zz;
    M66(2,5) = 2*zx*zy;

    M66(3,0) = yx*zx;
    M66(3,1) = yy*zy;
    M66(3,2) = yz*zz;

    M66(4,0) = xx*zx;
    M66(4,1) = xy*zy;
    M66(4,2) = xz*zz;

    M66(5,0) = xx*yx;
    M66(5,1) = xy*yy;
    M66(5,2) = xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

void rot_4th(double A[3][3][3][3], Matrix3d E, double B[3][3][3][3])
{
    double DUMMY;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int m = 0; m < 3; m++)
                for(int n = 0; n < 3; n++)
    {
    DUMMY = 0;
    for(int i1 = 0; i1 < 3; i1++)
        for(int j1 = 0; j1 < 3; j1++)
            for(int m1 = 0; m1 < 3; m1++)
                for(int n1 = 0; n1 < 3; n1++)
    DUMMY += E(i,i1)*E(j,j1)*E(m,m1)*E(n,n1)*A[i1][j1][m1][n1];
    B[i][j][m][n] = DUMMY;
    }
}

//convention between voigt C66 and Bbasis C66orC55
Matrix6d Btovoigt(MatrixXd B)
{
    double C4[3][3][3][3] = {0};
    Chg_basis(B,C4);
    return voigt(C4);
}

Matrix5d voigttoB5(Matrix6d M)
{
    double C4[3][3][3][3] = {0};
    voigt(M,C4);
    return Chg_basis5(C4);
}

Matrix6d voigttoB6(Matrix6d M)
{
    double C4[3][3][3][3] = {0};
    voigt(M,C4);
    return Chg_basis6(C4);
}

//definition of Rodrigues
Matrix3d Rodrigues(Matrix3d C)
{
    Matrix3d AROT;
//    Matrix3d Iij3; Iij(Iij3);
    Vector3d V;
    V(0)=C(2,1); V(1)=C(0,2); V(2)=C(1,0);
    double Vnorm = V.norm();
    if(Vnorm < 1e-30)
        AROT = Matrix3d::Identity();
    else
    {
        double Coef1 = sin(Vnorm)/Vnorm;
        double Coef2 = (1.0 - cos(Vnorm))/(Vnorm*Vnorm);
        AROT = Matrix3d::Identity() + Coef1*C + Coef2*C*C;
    }
    return AROT;
}

Matrix6d Bbasisadd(Matrix6d M66,Matrix5d M55)
{
    Matrix6d Mout;
    Mout = M66;
    Mout.block(0,0,5,5) += M55;
    return Mout;
}

Vector6d Bbasisadd(Vector6d V6,Vector5d V5)
{
    Vector6d Vout;
    Vout = V6;
    Vout.head(5) += V5;
    return Vout;
}

int Jacobi(Matrix3d A,Vector3d &D,Matrix3d &V)
{
    int n = 3;
    V = Matrix3d::Identity();
    Vector3d B,Z;
    for(int i = 0; i < n; i++)
    {
        B(i) = A(i,i);
        D(i) = B(i);
        Z(i) = 0;
    }
    int Nrot = 0;

    double Sm, Tresh, G, H, T, theta;
    double C, S, Tau;

    for(int i = 0; i < 50; i++)
    {
        Sm = 0;
        for(int ip = 0; ip < n-1; ip++)
            for(int iq = ip+1; iq < n; iq++)
                Sm += abs(A(ip,iq));
        if(Sm == 0) return 0;

        Tresh = 0;
        if(i < 4) Tresh = 0.2*Sm/n/n;

        for(int ip = 0; ip < n-1; ip++)
            for(int iq = ip+1; iq < n; iq++)
            {
                G=100*abs(A(ip,iq));
                if((i > 4)&&(abs(D(ip))+G == abs(D(ip)))&&(abs(D(iq))+G == abs(D(iq))))
                {A(ip,iq) = 0;}    
                else if(abs(A(ip,iq))>Tresh)
                {
                    H = D(iq)-D(ip);
                    if(abs(H)+G == abs(H))
                    {T=A(ip,iq)/H;}
                    else
                    {
                    theta = 0.5*H/A(ip,iq);
                    T=1.0/(abs(theta)+sqrt(1.0+theta*theta));
                    if(theta < 0) T=-T;
                    }
                    C=1.0/sqrt(1+T*T);
                    S=T*C;
                    Tau=S/(1.0+C);
                    H=T*A(ip,iq);
                    Z(ip)=Z(ip)-H;
                    Z(iq)=Z(iq)+H;
                    D(ip)=D(ip)-H;
                    D(iq)=D(iq)+H;
                    A(ip,iq)=0.0;
                    for(int j = 0; j < ip-1; j++)
                    {
                        G=A(j,ip);
                        H=A(j,iq);
                        A(j,ip)=G-S*(H+G*Tau);
                        A(j,iq)=H+S*(G-H*Tau);
                    }
                    for(int j = ip+1; j < iq-1; j++)
                    {
                        G=A(ip,j);
                        H=A(j,iq);
                        A(ip,j)=G-S*(H+G*Tau);
                        A(j,iq)=H+S*(G-H*Tau);
                    }
                    for(int j = iq+1; j < n; j++)
                    {
                        G=A(ip,j);
                        H=A(iq,j);
                        A(ip,j)=G-S*(H+G*Tau);
                        A(iq,j)=H+S*(G-H*Tau);
                    }
                    for(int j = 0; j < n; j++)
                    {
                        G=V(j,ip);
                        H=V(j,iq);
                        V(j,ip)=G-S*(H+G*Tau);
                        V(j,iq)=H+S*(G-H*Tau);
                    }
                    Nrot += 1;
                }
            }   
        for(int ik = 0; ik < n; ik++) 
        {
            B(ik)=B(ik)+Z(ik);
            D(ik)=B(ik);
            Z(ik)=0;
        } 
    }
    return 1;
}