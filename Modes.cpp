#include "Modes.h"

using namespace Modes;

int mode::ini_sn_mode(MatrixXd Min, int flag, int system_n)
{
    mtype = flag;
    system_num = system_n;
    if(mtype == 1)
    {
        sl = new DSystems_Voce::slip[system_num];
        for(int i = 0; i < system_num; i++)
            sl[i].ini_sn_s(Min.row(i));
    }
    if(mtype == 0)
    {
        tw = new DSystems_Voce::twin[system_num];
        for(int i = 0; i < system_num; i++)
            tw[i].ini_sn_s(Min.row(i));
    }
    return 0;
}

int mode::check_sn_mode()
{
    if(mtype == 1)
    {
        cout << "Slip systems\n";
        for(int i = 0; i < system_num; i++)
            sl[i].check_sn_s();
    }
    if(mtype == 0)
    {
        cout << "Twinning systems\n";
        for(int i = 0; i < system_num; i++)
            tw[i].check_sn_s();
    }
    return 0;
}

int mode::ini_hardening_mode(double nrsx_in, VectorXd CRSS_p_in, VectorXd hst_in)
{
    nrsx = nrsx_in;
    hst = hst_in;
    if(mtype == 1)
    {
        for(int i = 0; i < system_num; i++)
            sl[i].ini_hardening_s(CRSS_p_in);
    }
    if(mtype == 0)
    {
        for(int i = 0; i < system_num; i++)
            tw[i].ini_hardening_s(CRSS_p_in);
    }    
    return 0;
}

int mode::check_hardening_mode()
{
    if(mtype == 1)
    {
        cout << "Slip systems\n";
        for(int i = 0; i < system_num; i++)
            sl[i].check_hardening_s();
    }
    if(mtype == 0)
    {
        cout << "Twinning systems\n";
        for(int i = 0; i < system_num; i++)
            tw[i].check_hardening_s();
    }
    return 0;    
}

double mode::get_gamma0(){return gamma0;}
double mode::get_nrsx(){return nrsx;}

double mode::get_RSSxM(Matrix3d sig)
{
    double temp;
    double RSSx = 0;
    for(int i = 0; i < system_num; i++)
    {
        if(mtype == 1) temp = abs(sl[i].cal_RSSx(sig));
        if(mtype == 0) temp = abs(tw[i].cal_RSSx(sig));
        if(RSSx < temp) RSSx = temp;
    }
    return RSSx;
}

Matrix3d mode::cal_dijpmode(Matrix3d sig)
{
    Matrix3d Pij = Matrix3d::Zero();
    for(int i = 0; i < system_num; i++)
    {
        if(mtype == 1) Pij += sl[i].cal_dijp_a(sig, nrsx);
        if(mtype == 0) Pij += tw[i].cal_dijp_a(sig, nrsx);
    }
    //cout << "Pij total\t" << Pij << endl;
    return Pij;
}

Matrix3d mode::cal_rotslip_m()
{
    Matrix3d wij = Matrix3d::Zero();
    for(int i = 0; i < system_num; i++)
    {
        if(mtype == 1) wij += sl[i].cal_rotslip_a();
        if(mtype == 0) wij += tw[i].cal_rotslip_a();
    }
    return wij;
}

Matrix6d mode::get_Fgradm(Matrix3d sig)
{
    Matrix6d Pij = Matrix6d::Zero();
    for(int i = 0; i < system_num; i++)
    {
        if(mtype == 1) Pij += sl[i].get_Fgrads(sig, nrsx);
        if(mtype == 0) Pij += tw[i].get_Fgrads(sig, nrsx);
    }
    return Pij;
}

double mode::Update_shear_strain_m()
{
    gamma_rate_abs_m = 0;
    for(int i = 0; i < system_num; i++)
    {
        if(mtype == 1) gamma_rate_abs_m += sl[i].get_gammarate();
        if(mtype == 0) gamma_rate_abs_m += tw[i].get_gammarate();
    }
    return gamma_rate_abs_m;
}

void mode::Update_CRSS_m(double Tincr, double gamma_total, double gamma_delta, double gamma_delta_gmode[], int mode_num)
{
    //refer to Eq[C-12]
    double Hstpart = 0;
    for(int i = 0; i < mode_num; i++)
        Hstpart += hst(i)*gamma_delta_gmode[i];
    Hstpart = Hstpart/gamma_delta;
    for(int i = 0; i < system_num; i++)
    {
        if(mtype == 1) sl[i].Update_CRSS_a(Hstpart, gamma_delta, gamma_total);
        if(mtype == 0) tw[i].Update_CRSS_a(Hstpart, gamma_delta, gamma_total);
    }
}