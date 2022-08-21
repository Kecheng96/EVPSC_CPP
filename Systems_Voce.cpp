#include "Systems_Voce.h"

using namespace DSystems_Voce;

//
int dsystem::ini_sn_s(VectorXd vin)
{   
    n = vin(seq(0,2)); //normal 
    s = vin(seq(3,5)); //Burgers
    Pij=0.5*(s*n.transpose()+n*s.transpose());
    Rij=0.5*(s*n.transpose()-n*s.transpose());
    return 0;
}

int dsystem::check_sn_s()
{
    cout << n.transpose() << "\t" << s.transpose() << endl;
    cout <<  "Pij\n" << Pij << endl;
    return 0;
}

int dsystem::ini_hardening_s(VectorXd CRSS_p_in)
{
    tau0 = CRSS_p_in(0);
    tau1 = CRSS_p_in(1);
    theta0 = CRSS_p_in(2);
    theta1 = CRSS_p_in(3);

    aratio = theta0/tau1;
    CRSS = tau0;
    return 0;
}

int dsystem::check_hardening_s()
{
    Vector4d vout(tau0, tau1, theta0, theta1);
    cout << vout.transpose() << endl;
    return 0;
}

double dsystem::cal_RSSx(Matrix3d sig)
{
    double RSS = 0;
    for(int i = 0; i < 3; i++ )
        for(int j = 0; j < 3; j++ )
            RSS += Pij(i,j) * sig(i,j);

    return RSS/CRSS;
}

Matrix3d dsystem::cal_dijp_a(Matrix3d sig, double nrsx)
{
    double RSSx = cal_RSSx(sig);
    gammarate = gamma0 * pow(abs(RSSx), nrsx-1) * RSSx;
    return gammarate * Pij;
}

Matrix3d dsystem::cal_rotslip_a()
{
    return gammarate * Rij;
}

Matrix6d dsystem::get_Fgrads(Matrix3d sig, double nrsx)
{
    double RSSx = cal_RSSx(sig);

    Vector6d Pijv = voigt(Pij);
    Matrix6d Mout = Pijv * Pijv.transpose();
    
    return Mout * gamma0 * pow(abs(RSSx), nrsx-1) / CRSS * nrsx;
}

double dsystem::get_gammarate(){return abs(gammarate);}

void dsystem::Update_CRSS_a(double Hstpart,double gamma_delta, double gamma_total)
{
    double exp1 = exp(-aratio*gamma_total);
    double exp2 = exp(-aratio*gamma_delta);
    CRSS_delta = Hstpart * ( \
                 theta1*gamma_delta \
               - (aratio*tau1-theta1)/aratio*exp1*(exp2-1.0) \
               - theta1/aratio*exp1*(exp2*(1+aratio*(gamma_total+gamma_delta))\
               - (1.0 + aratio*gamma_total)) \
               );
    CRSS += CRSS_delta;
}