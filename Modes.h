#ifndef MODES_H
#define MODES_H

#include "Systems_Voce.h"
#include "Toolbox.h"

namespace Modes{

//machnical system from which slip&twinning system can be derived 
class mode
{
    private:
        int mtype = 1; //type of deformation modes (1: slip; 0: twin)
        int system_num = 0; 
        DSystems_Voce::slip* sl = NULL;
        DSystems_Voce::twin* tw = NULL;

        double gamma0 = 1e-3; //reference shear strain rate

        double gamma_rate_abs_m = 0; //shear strain of one mode

        double nrsx; //rate sensitive

        VectorXd hst; //hst
        
    public:
        int ini_sn_mode(MatrixXd, int, int);
        //input the normal and Burgers vector in each system
        //(MatrixXd) sn, (int) flag of twin or slip, (int) number of systems
        //need loop over systems
        int check_sn_mode();

        int ini_hardening_mode(double, VectorXd, VectorXd);
        //input the hardening parameters
        //input parameters:
        //double nrsx_in; VectorXd CRSS_p_in; VectorXd hst_in
        int check_hardening_mode();

        Matrix3d cal_dijpmode(Matrix3d);

        Matrix3d cal_rotslip_m();
        
        Matrix6d get_Fgradm(Matrix3d);

        double get_gamma0();
        double get_nrsx();
        double get_RSSxM(Matrix3d);

        double Update_shear_strain_m();

        //update the CRSS in each mode
        //parameters:
        //double Tincr, double gamma_total, double gamma_delta,
        //double gamma_delta_gmode[], int modes_num
        void Update_CRSS_m(double, double, double, double gamma_delta_gmode[], int);
};

}

#endif