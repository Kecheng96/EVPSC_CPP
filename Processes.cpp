#include "Processes.h"
using namespace Procs;

Process::Process()
{
}

void Process::load_ctrl(Vector4d Vin)
{
    Nsteps = int(Vin(0));
    Ictrl = int(Vin(1)) - 1;
    Eincr = Vin(2);
    Temp = Vin(3);
}
void Process::get_Udot(Matrix3d Min)
{   
    Udot_input = Min;
    Ddot_input = 0.5*(Min + Min.transpose());
    Wdot_input = 0.5*(Min - Min.transpose());
    
    //calculate Time increment Tincr
    Vector6d Vtemp = voigt(Ddot_input);
    Tincr = Eincr / Vtemp(Ictrl);
}

void Process::set_IUdot(Matrix3i Min)
{   
    IUdot = Min;
    IDdot(0)=IUdot(0,0);
    IDdot(1)=IUdot(1,1);
    IDdot(2)=IUdot(2,2);
    IDdot(3)=IUdot(1,2)*IUdot(2,1);
    IDdot(4)=IUdot(0,2)*IUdot(2,0);
    IDdot(5)=IUdot(0,1)*IUdot(1,0);
}

void Process::set_ISdot(Vector6i Min){ISdot = Min;}

void Process::loading(Polycs::polycrystal &pcrys)
{
    double SC_err_m = 0.01;
    int SC_iter_m = 20;
    double errD_m_AV = 0.01;
    double errS_m_AV = 0.01;
    double err_g_AV = 0.01;
    double errd, errs, err_g;
}