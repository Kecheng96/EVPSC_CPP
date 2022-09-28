#ifndef POLYCRYSTALS_H
#define POLYCRYSTALS_H

#include <omp.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

#include "Grains.h"
#include "Toolbox.h"

namespace Polycs{

//const int  max_grain_n = 3000; 

class polycrystal
{
    private:
        Grains::grain* g  = NULL; 
        int grains_num = 0;

        double GZ = 30; //grain size (um)

        Gausspoint* Gpsets  = NULL; 

        string crysym; //crystal symmetry
        MatrixXd Trans_Miller;// the conversion matrix of Miller indices
        int Miller_n = 3;// the number of the Miller indices (3 or 4(hcp))
        Matrix3d Mabc; 

        Vector3d Cdim; //lattice constant
        Vector3d Cang; //lattice constant
        Vector6d therm; //Thermal expansion coefficients [K^(-1)]
        double Temp; // /(K^-1)

        //the shape of ellipsoid
        //Vector3d ell_axis_o;
        Vector3d ell_axis; //axial length of ellipsoid/grain (_o means the original value)
        Matrix3d ell_axisb;
        double ell_crit_shape = 25;
        Vector3d ellip_ang; //the rotate angle of the ellipsoid
        bool Iflat = 0; //flag of stretch the ellipsoid (0: yes; 1: no)
        bool Ishape = 0; //flag of individual ellipsoid for grain (1: yes; 0: no)
        ///////

        double Cijkl[3][3][3][3]; //elastic constants 3X3X3X3
        Matrix6d Cij6; //elastic constants 6X6 of crystal
        //ESC
        Matrix6d CSC; //The elastic consistent stiffness
        Matrix6d COLD;//The initial guessed elastic consistent stiffness
        Matrix6d SSC; //CSC^-1
        //VPSC
        Matrix5d C_VP_SC; //The visco-plastic stiffness 
        Matrix5d M_VP_SC; //The visco-plastic compliance C_VP_SC^-1
        Vector6d D0; //the macro back-extrapolated term (follow Equ[5-41b])

        Vector5d DVP_AV;

        Matrix3d Fij_m; //the macro deformation tensor in grain
        Matrix3d Udot_m; //the macro velocity gradient;
        Matrix3i IUdot; //the flag of known (control by Udot_m) and unknown (calculated by EVPSC) velocity components
        Vector6i IDdot; //the flag of known and unknown strain rate components
        Vector6i ISdot; //the flag of known and unknown stress componets
    
        Matrix3d Dij_m; //the macro strain rate tensor
        Matrix3d Wij_m; //the macro rotation rate tensor
        Matrix3d Udot_AV;
        // Udot_m = Dij_m + Wij_m
        Matrix3d Dij_AV; //the average strain rate tensor of all garins
        Matrix3d Dije_AV; //the elastic part
        Matrix3d Dijp_AV; //the vp part

        Matrix3d Eps_m; //macro strain
        Matrix3d Sig_m; //macro stress 
        Matrix3d Sig_AV; //the average stress tensor of all garins
        Matrix3d Sig_m_old; //macro stress in last increment

        //some parameters for error control and iteration control
        double SC_err_m = 0.01; //the error limit of Self-consistent compliance or stiffness
        int SC_iter_m = 20; //the max iteration number of SC
        double errD_m = 0.01; //error limit between the input macro strain rate and output at each iteration 
        double errS_m = 0.01; //error limit of the macro stress
        double err_g_AV = 0.01; //error limit of the average grain stress
        
        Matrix3d Sig_in;
        Matrix3d Dij_in;
        Matrix3d sig_in_AV;

        Matrix6d Msup;

    public:
        polycrystal();

        void ini_Udot_m(Matrix3d);
        void ini_Sig_m(Matrix3d);
        void set_IUdot(Matrix3i);
        void set_ISdot(Vector6i);

        //input the number of grains
        int grains_n(int);  
        int check_grains_n();
        
        //input the euler angles and weights
        //needs loop over grains
        void ini_euler(Vector4d, int); 
        void Norm_weight();

        int ini_cry(string, VectorXd);
        //input the crystal constant
        void check_cry();
        int get_Millern();
        //output the need number of Miller indice 

        void ini_Cij6(MatrixXd);
        //input the elastic constant in 6X6 tensor
        int check_Cij6();

        int ini_therm(VectorXd);
        //input the thermal expansion coefficients [K^(-1)]
        int check_therm();

        int ini_gmode(int n);
        //input the deformation modes
        //needs loop over grains
        int check_gmode();
             
        //input the normal and Burgers vector of deformation system in a mode
        //input parameters:
        //(MatrixXd) sn, (int) flag of twin or slip, (int) number of systems,(int) mode label
        //needs loop over grains
        int ini_sn(MatrixXd, int, int, int);
        int check_sn();

        int ini_GZ(double);
        //input the grain size;

        int ini_hardening(double, VectorXd, VectorXd, int);
        //input the hardening parameters
        //input parameters:
        //double nrsx_in; VectorXd CRSS_p_in; VectorXd hst_in; int modei
        int check_hardening();

        int Update_Fij(double);
        int Update_shape();
        Vector3d get_ell_axis(); 
        Vector3d get_ellip_ang(); 

        //the singular step according to a certain process
        int EVPSC(int, double, bool, bool, bool);

        //calculate the macro&grain elastic compliance
        int Selfconsistent_E(int, double, int); 

        //calculate the macro&grain VP compliance
        int Selfconsistent_P(int, double, int); 

        // calculate the Sig_m and Dij_m
        void Cal_Sig_m(double);  

        // calculate the Sig_g and Dij_g including the elastic and vp part
        void Cal_Sig_g(double);

        void Update_AV(); //update the volume average value

        //output
        Vector6d get_Sig_m();
        Vector6d get_Eps_m();
        void get_euler(fstream &);

        
};


}
#endif