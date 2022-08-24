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
    
    //calculate Time increment Tincr
    Vector6d Vtemp = voigt(Ddot_input);
    Tincr = Eincr / Vtemp(Ictrl);
}
void Process::get_Sdot(Matrix3d Min){Sdot_input = Min;}
void Process::get_IUdot(Matrix3i Min){IUdot = Min;}
void Process::get_ISdot(Vector6i Vin){ISdot = Vin;}

void Process::loading(Polycs::polycrystal &pcrys)
{
    pcrys.ini_Udot_m(Udot_input);
    pcrys.ini_Sig_m(Sdot_input);
    pcrys.set_IUdot(IUdot);
    pcrys.set_ISdot(ISdot);
    for(int istep = 0; istep < Nsteps; ++istep)
    {
        pcrys.EVPSC(istep, Tincr);
    }

}