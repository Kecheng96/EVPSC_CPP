#ifndef SYSTEMS_VOCE_H
#define SYSTEMS_VOCE_H

#include "Toolbox.h"

namespace DSystems_Voce{
//deformation system from which slip&twinning system can be derived 
class dsystem
{
    protected:
        double gamma0 = 1e-3; //reference shear strain rate
        double gammarate = 0;

        double tau0,tau1;//the initial CRSS and the back-extrapolated CRSS(tau0+tau1)
        double theta0,theta1;//the initial hardening rate, the asymptotic hardening rate
        double aratio = 0; // = theta0/tau1;
        double CRSS;
        double CRSS_delta;

        //double nrsx; //rate sensitivity component
        Vector3d s,n; //the Burgers vector(s) and normal(n) of plane
        Matrix3d Pij;
        Matrix3d Rij;
    
    public:

        //initial the sn and calculate the pij
        int ini_sn_s(VectorXd);
        int check_sn_s();

        //input the harderning parameter tau0, tau1, theta0, theta1;
        int ini_hardening_s(VectorXd);
        int check_hardening_s();

        double cal_RSSx(Matrix3d);

        double get_gammarate();
        
        //calculate the gamma_rate * Pij to get the vp strain rate Dijp_g
        //Parametes:
        //Matrix3d sig, double nrsx
        Matrix3d cal_dijp_a(Matrix3d, double);

        //calculate the gamma_rate * Rij
        Matrix3d cal_rotslip_a();

        //update CRSS of system
        //parameters:
        //double Hstpart, double gamma_delta, double gamma_total
        void Update_CRSS_a(double, double, double);
        
        Matrix6d get_Fgrads(Matrix3d, double);
};

class slip : public dsystem
{
};

class twin : public dsystem
{
};


}
#endif