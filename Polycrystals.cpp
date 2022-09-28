#include "Polycrystals.h"
using namespace Polycs;
using namespace std;

polycrystal::polycrystal()
{

    //initial the macro stress&strain
    Eps_m = Matrix3d::Zero();
    Sig_m = Matrix3d::Zero();

    //initial the shape of ellipsoid
    ell_axis = Vector3d::Ones();

    ellip_ang << 90,90,90;
    ell_axisb = Euler_trans(ellip_ang);

    Fij_m = Matrix3d::Identity();


    //initial the VP consistent
    M_VP_SC = Matrix5d::Identity();
    M_VP_SC = 1e-8 * M_VP_SC;
    C_VP_SC = M_VP_SC.inverse();
    D0 = Vector6d::Zero();

    Msup<<1,0,0,0,0,0,
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,2,0,0,
    0,0,0,0,2,0,
    0,0,0,0,0,2;

    Vector10d xth,xph,wth,wph;
    //integral points and weights
    Integralpoint3 alpha, aww;
    Integralpoint6 aa6, aaww6; //coordinate and weigts in Fourier space 
    Integralpoint1 ww;

    Gpsets = new Gausspoint[11];

    for(int Gpcase = 0; Gpcase < 11; Gpcase++ ){
        switch(Gpcase){
            case 0:
            xth << 4.71236594e-02,0.241774723e0,0.565131843e0,0.968887568e0,1.37937832e0,
            1.76221442e0,2.17270517e0,2.57646084e0,2.89981818e0,3.09446883e0;
            wth << 0.120191820e0,0.264987558e0,0.373805553e0,0.420841277e0,0.390970200e0,
            0.390970260e0,0.420841366e0,0.373805553e0,0.264987499e0,0.120192111e0;
            break;
            case 1:
            xth << 1.57080423e-02,0.144995824e0,0.425559640e0,0.829968274e0,1.31460333e0,
            1.82698941e0,2.31162453e0,2.71603298e0,2.99659705e0,3.12588477e0;
            wth << 5.41692823e-02,0.207461149e0,0.348739326e0,0.452716887e0,0.507709801e0,
            0.507709682e0,0.452716798e0,0.348738998e0,0.207461327e0,5.41692935e-02;
            break; 
            case 2:
            xth << 3.76990959e-02,0.198626831e0,0.483041346e0,0.871647120e0,1.32964790e0,
            1.81194484e0,2.26994562e0,2.65855122e0,2.94296598e0,3.10389376e0;
            wth << 9.68142375e-02,0.224478707e0,0.341134071e0,0.430180043e0,0.478189558e0,
            0.478189170e0, 0.430180043e0, 0.341134191e0, 0.224478647e0, 9.68143344e-02;
            break;
            case 3:
            xth << 3.45576368e-02,0.187556863e0,0.468425453e0,0.859980166e0,1.32527423e0,
            1.81631863e0,2.28161263e0,2.67316723e0,2.95403576e0,3.10703516e0;
            wth << 8.95763785e-02,0.217725381e0,0.341026783e0,0.435772508e0,0.486694932e0,
            0.486695170e0,0.435772508e0,0.341026902e0,0.217725128e0,8.95764604e-02;
            break;
            case 4:
            xth << 3.14158052e-02,0.177928671e0,0.457155794e0,0.851592362e0,1.32222414e0,
            1.81936860e0,2.29000044e0,2.68443704e0,2.96366405e0,3.1101768e0;
            wth << 8.26927349e-02,0.213228315e0,0.342008322e0,0.440196186e0,0.492670894e0,
            0.492670983e0,0.440195888e0,0.342008322e0, 0.213227972e0, 8.26930404e-02;
            break;
            case 5:
            xth << 2.98452154e-02,0.173592165e0,0.452448040e0,0.848216832e0,1.32101476e0,
            1.82057810e0,2.29337597e0,2.68914461e0,2.96800065e0,3.11174774e0;
            wth << 7.93928578e-02,0.211627841e0,0.342669785e0,0.442057431e0,0.495048553e0,
            0.495048642e0,0.442057490e0,0.342670023e0,0.211627468e0,7.93929026e-02;
            break;
            case 6:
            xth << 2.67036632e-02,0.165752888e0,0.444431901e0,0.842614472e0,1.31902647e0,
            1.82256627e0,2.29897833e0,2.69716072e0,2.97583985e0,3.11488938e0;
            wth << 7.30879456e-02,0.209402516e0,0.344104946e0,0.445234656e0,0.498966068e0,
            0.498966306e0,0.445234746e0, 0.344104946e0,0.209402665e0,7.30878562e-02;
            break;
            case 7:
            xth << 2.67036632e-02,0.165752888e0,0.444431901e0,0.842614472e0,1.31902647e0,
            1.82256627e0,2.29897833e0,2.69716072e0,2.97583985e0,3.11488938e0;
            wth << 7.30879456e-02,0.209402516e0,0.344104946e0,0.445234656e0,0.498966068e0,
            0.498966306e0,0.445234746e0,0.344104946e0,0.209402665e0,7.30878562e-02;
            break;
            case 8:
            xth <<2.43473575e-02,0.160516247e0,0.439386278e0,0.839168847e0,1.31781363e0,
            1.82377899e0,2.30242372e0,2.70220637e0,2.98107672e0,3.11724544e0;
            wth << 6.86219111e-02,0.208388865e0,0.345189095e0,0.447236270e0,0.501360059e0,
            0.501359940e0,0.447236151e0,0.345189214e0,0.208388969e0,6.86219335e-02;
            break;
            case 9:
            xth << 2.19910536e-02,0.155757755e0,0.434985727e0,0.836206555e0,1.31677616e0,
            1.82481658e0,2.30538607e0,2.70660710e0,2.98583508e0,3.11960149e0;
            wth << 6.43825606e-02,0.207786217e0,0.346235514e0,0.448981822e0,0.503410578e0,
            0.503410578e0,0.448981792e0,0.346235693e0,0.207785636e0,6.43827692e-02;
            break;
            case 11:
            xth << 2.04204638e-02,0.152822554e0,0.432348520e0,0.834448099e0,1.31616223e0,
            1.82543063e0,2.30714464e0,2.70924401e0,2.98877001e0,3.12117243e0;
            wth << 6.16818815e-02,0.207559645e0,0.346902698e0,0.450027168e0,0.504624724e0,
            0.504624426e0,0.450027317e0,0.346902847e0,0.207559645e0,6.16819337e-02;
            break;
        }
        wph = wth;
        xph = xth;	
        double sinth, costh, simbtet;
        Matrix3d aa, aaww;
        int ny;
        for (int ith = 0; ith < Intn; ith++){
            sinth = sin(xth(ith));
            costh = cos(xth(ith));
            simbtet = wth(ith) * sinth / (2.0*M_PI);

            for (int iph = 0; iph < Intn; iph++){
                ny = iph + ith * Intn;
                ww(ny) = simbtet*wph(iph);
                alpha(0,ny) = sinth*cos(xph(iph));
                alpha(1,ny) = sinth*sin(xph(iph));
                alpha(2,ny) = costh;
                //
                for(int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++){
                        aa(i,j) = alpha(i,ny) * alpha(j,ny);
                        aaww(i,j) = aa(i,j) * ww(ny);
                    }
                aa6.col(ny) = voigt(aa);
                aaww6.col(ny) = voigt(aaww);
                for(int i = 0; i < 3; i++)
                    aww(i,ny) = alpha(i,ny) * ww(ny);
            }
        }
        Gpsets[Gpcase].Gpaa6 = aa6;
        Gpsets[Gpcase].Gpaaww6 = aaww6;
        Gpsets[Gpcase].Gpalpha = alpha;
        Gpsets[Gpcase].Gpaww = aww;
        Gpsets[Gpcase].Gpww = ww;
    }
}


void polycrystal::ini_Udot_m(Matrix3d Udot_input)
{   
    Udot_m = Udot_input;
    Dij_m = 0.5*(Udot_m + Udot_m.transpose());
    Wij_m = 0.5*(Udot_m - Udot_m.transpose());
    
    Dij_AV = Dij_m;
    Dije_AV = Dij_m;
    Dijp_AV = Matrix3d::Zero();
}

void polycrystal::ini_Sig_m(Matrix3d Min){Sig_m = Min;}

void polycrystal::set_IUdot(Matrix3i Min)
{   
    IUdot = Min;
    IDdot(0)=IUdot(0,0);
    IDdot(1)=IUdot(1,1);
    IDdot(2)=IUdot(2,2);
    IDdot(3)=IUdot(1,2)*IUdot(2,1);
    IDdot(4)=IUdot(0,2)*IUdot(2,0);
    IDdot(5)=IUdot(0,1)*IUdot(1,0);
}
void polycrystal::set_ISdot(Vector6i Min){ISdot = Min;}


int polycrystal::grains_n(int n)
{
    grains_num = n;
    g = new Grains::grain[n];
    //change the number of grains
    for(int i = 0; i < n; i++)  g[i].grain_i = i;
    return 0;
}

int polycrystal::check_grains_n()
{
    cout << "the number of grains:\n" << grains_num << endl;
    return grains_num;
    //print the number of grains
}

void polycrystal::ini_euler(Vector4d vin, int i){g[i].ini_euler_g(vin);} 

void polycrystal::Norm_weight()
{
    double total_w = 0;
    for(int i = 0; i < grains_num; i++)
        total_w += g[i].get_weight_g();
    
    for(int i = 0; i < grains_num; i++)
        g[i].set_weight_g(g[i].get_weight_g()/total_w);

}

int polycrystal::ini_cry(string strin, VectorXd vin)
{
    //transform str to lower case
    for (int i = 0; i < strin.size(); i++)
		strin[i] = tolower(strin[i]);
    crysym = strin;    
    Cdim = vin(seq(0,2));
    Cang = vin(seq(3,5)) / 180 * M_PI; //Converting degrees to radians
    
    //calculate conversion matrix of Miller indices according to the crysym
    if(!crysym.compare("hexag"))
    {
        Miller_n = 4;
        MatrixXd Mtemp(3,4);
        Trans_Miller = Mtemp;
        Trans_Miller <<
        1, 0, -1, 0,
        0, 1, -1, 0,
        0, 0,  0, 1;
    }
    else if(!crysym.compare("cubic")) 
    {
        Miller_n = 3;
        MatrixXd Mtemp(3,3);
        Trans_Miller = Mtemp;
        Trans_Miller <<
        1, 0,  0,
        0, 1,  0,
        0, 0,  1;
    }

    Mabc(0,0)=sin(Cang(1));
    Mabc(1,0)=0.;
    Mabc(2,0)=cos(Cang(1));
    Mabc(0,1)=(cos(Cang(2))-cos(Cang(0))*cos(Cang(1)))/sin(Cang(1));
    Mabc(2,1)=cos(Cang(0));
    Mabc(1,1)=sqrt(1.0-pow(Mabc(0,1),2)-pow(Mabc(2,1),2));
    Mabc(0,2)=0.;
    Mabc(1,2)=0.;
    Mabc(2,2)=1.;

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            Mabc(i,j) = Cdim(j) * Mabc(i,j);

    return 0;
}

int polycrystal::get_Millern(){return Miller_n;}

void polycrystal::check_cry()
{
    cout << crysym << endl;
    cout << "the cdim and cang:" << endl;
    cout << Cdim.transpose() << endl;
    cout << Cang.transpose() << endl;
}

void polycrystal::ini_Cij6(MatrixXd Min)
{   
    Cij6 = Min;
    voigt(Cij6,Cijkl);
}

int polycrystal::check_Cij6()
{
    cout << "elastic constant:\n" << Cij6 << "\n";
    return 0;
}

int polycrystal::ini_therm(VectorXd vin)
{
    therm = vin;
    return 0;
}

int polycrystal::check_therm()
{
    cout << "Thermal coefficient:\n" << therm.transpose() << "\n";
    return 0;
}

int polycrystal::ini_gmode(int n)
{
    for(int i = 0; i < grains_num; i++)
    {
        g[i].ini_gmode_g(n);
    }
    return 0;
}

int polycrystal::check_gmode()
{
    for(int i = 0; i < grains_num; i++)
    {
        cout << "the number of modes in Grain " << i << ":\n";
        cout  << g[i].check_gmode_g() << endl;
    }
    return 0;
}

int polycrystal::ini_sn(MatrixXd Min, int flag, int system_n, int modei)
{
    MatrixXd Min_s, Min_n;

    Min_n = Min(all,seq(0,Miller_n-1))*Trans_Miller.transpose();
    Min_s = Min(all,seq(Miller_n,2*Miller_n-1))*Trans_Miller.transpose(); 
    //calculate the coordinate in Cartesian system
    Min_n = Min_n*Mabc.transpose();
    Min_s = Min_s*Mabc.transpose();

    //normalization
    for(int i = 0; i < system_n; i++)
    {
        Min_n.row(i) = Min_n.row(i).normalized();
        Min_s.row(i) = Min_s.row(i).normalized();
    }

    MatrixXd Min_ns(system_n,6);
    Min_ns.block(0,0,system_n,3) = Min_n;
    Min_ns.block(0,3,system_n,3) = Min_s;

    for(int i = 0; i < system_n; i++)
        for(int j = 0; j < 6; j++)
            if(abs(Min_ns(i,j)) <= 1e-3 ) Min_ns(i,j) = 0.0;

    for(int i = 0; i < grains_num; i++)
    {
        g[i].ini_sn_g(Min_ns, flag, system_n, modei);
    }
    return 0;
}

int polycrystal::check_sn()
{
    for(int i = 0; i < grains_num; i++)
    {
        cout << "the deformation system in Grain " << i << ":\n";
        g[i].check_sn_g();
    }
    return 0;
}

int polycrystal::ini_GZ(double x)
{
    GZ = x;
    return 0;
}

int polycrystal::ini_hardening(double nrsx_in, VectorXd CRSS_p_in, VectorXd hst_in, int modei)
{
    for(int i = 0; i < grains_num; i++)
    {
        g[i].ini_hardening_g(nrsx_in, CRSS_p_in, hst_in, modei);
    }
    return 0;
}

int polycrystal::check_hardening()
{
    for(int i = 0; i < grains_num; i++)
    {
        cout << "the hardening parameters in Grain " << i << ":\n";
        g[i].check_hardening_g();
    }    
    return 0;
}

int polycrystal::Selfconsistent_E(int Istep, double ERRM, int ITMAX)
{
    //-1 Calculate the Elastic stiffness in Jaumann rate in all grains
    //(Rotate from crystal to Sample axes according to the euler angle)
    // and sum with the weight
    // the result is CUB
    Matrix6d CUB;
    CUB = Matrix6d::Zero();
    double C4SA[3][3][3][3];  // Elastic stiffness Rotate from crystal to Sample axes
    double C4SAS[3][3][3][3]; // ...in Jaumann rate
    Matrix3d sig_g;

    //cout << "\nCijkl:\n" << Cij6 << endl;

    double DUMMY = 0; //temporate variable in sum calculation;

    for(int G_n = 0; G_n < grains_num; G_n++)
    {
        // -2  Rotate the tensor Cijkl in grain to Sample Axes
        Matrix3d Euler_M = g[G_n].get_Euler_M_g();
        Matrix3d ET = Euler_M.transpose();
        voigt(rotate_C66(Cij6, ET), C4SA);
        g[G_n].Update_Cij6_SA_g(voigt(C4SA));
        // -2  Rotate the tensor Cijkl in grain to Sample Axes

        sig_g = g[G_n].get_stress_g();

        // -3 ( Eq[2.14] in Wang et al., 2010)
        // the elastic stiffness invovling Jaumann rate
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < 3; k++)
                    for(int l = 0; l < 3; l++)
        {
            C4SAS[i][j][k][l] = C4SA[i][j][k][l] - sig_g(i,j)*Iij(k,l);
        //                        +0.5*Iij(i,k)*sig_g(j,l) \
                                -0.5*Iij(i,l)*sig_g(j,l) \
                                -0.5*Iij(j,l)*sig_g(i,k) \
                                +0.5*Iij(j,k)*sig_g(i,l);
        }
        // -3 ( Eq[2.14] in Wang et al., 2010)
        // the elastic stiffness invovling Jaumann rate 

        g[G_n].Update_Mij6_J_g(Chg_basis6(C4SAS).inverse());
        //store the Jaumann rate elastic stiffness in grains

        CUB += Chg_basis6(C4SAS) * g[G_n].get_weight_g();

        // CUB is the volume average Elastic stiffness of all grains
    }
    //-1 Calculate the Elastic stiffness in Jaumann rate in all grains
    //(Rotate to Sample axes according to the euler angle)
    // and sum with the weight 
    // the result is CUB

    if(Istep == 0)  
        CSC = CUB;  //first step, use the volume average
    
    //-5
    //loop to make the guessed elastic stiffness CSC
    // to the Eshelby calculated CNEW 
    Matrix6d SSC_new;

    int IT = 0; //loop flag
    double RER = 2*ERRM;
    while((RER >= ERRM) && (IT < ITMAX))
    {
        IT++;
        SSC_new = Matrix6d::Zero();

        Chg_basis(CSC, C4SA);
        
        Matrix6d Ctilde;
        Matrix6d S66; //the Sijkl Equ[5-33] in sample axes in Manual 7d
        double R4_SA[3][3][3][3]; //PIijkl
        double RSinv_SA[3][3][3][3];
        Matrix6d S66inv; 
        Matrix3d axisb_t;
        Vector3d axis_t;
        Matrix6d Me_g; 

        //solve the eshelby tensor in the common ellipsoid
        if(Ishape == 0)
        {
            //-4
            //rotate the macro Elastic stiffness 
            //from Sample axes to the ellipsoid axes;
            // if Ishape = 0, which means the grains share the same ellipsoid axes 
            axisb_t = ell_axisb; 
            axis_t = ell_axis; 
            //Rotate the self consistent Elastic stiffness (CSC)
            //from sample axes to ellipsoid axes
            Matrix6d C66 = rotate_C66(voigt(C4SA), axisb_t.transpose());
            //-4

            //-9
            //calculate the Eshelby tensor
            double S4_EA[3][3][3][3] = {0}; 
            double R4_EA[3][3][3][3] = {0};
            
            int case_c = Eshelby_case(axis_t);
            g[0].Eshelby_E(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,Gpsets[case_c].Gpalpha);
            //-9

            //-10
            //rotate the eshelby tensor
            //back to the sample axes
            S66 = voigttoB6(rotate_C66(voigt(S4_EA), axisb_t));
            rot_4th(R4_EA, axisb_t, R4_SA);
            //-10

            //-11
            //Calculate Ctilde
            //refer to Equ[5-33] in Manual 7d
            //M~=(I-S)^-1 * S * M
            //but here we used C~ = (M~)^-1
            //that is 
            // C~ = C * (S^-1 - I)
            //-12 calculate (S)^-1
            S66inv = S66.inverse(); // 6x6 of S^-1
            Matrix6d S66inv_I; // 6x6 of (S^-1 - I)
            //for(int i = 0; i < 6; i++)
            //    for(int j = 0; j < 6; j++)
            //        S66inv_I(i,j) = S66inv(i,j) - Iij(i,j);
            S66inv_I = S66inv - Matrix6d::Identity();
            //-12
            // C~ = C * (S^-1 - I)
            Ctilde = CSC * S66inv_I;
            //-11
        }
        //end of Ishape == 0

        for(int G_n = 0; G_n < grains_num; G_n++)
        {
            //solve the eshelby tensor in all the ellipsoid of grain 
            if(Ishape == 1)
            {
                //-7
                //rotate the macro Elastic stiffness
                //from Sample axes to the ellipsoid axes;
                // if Ishape = 1, which means the ellipsoid axes varies in grains,
                axis_t = g[G_n].get_ell_axis_g();
                axisb_t = g[G_n].get_ell_axisb_g();
                Matrix6d C66 = rotate_C66(CSC, axisb_t.transpose());

                //-9
                //calculate the Eshelby tensor
                double S4_EA[3][3][3][3] = {0}; 
                double R4_EA[3][3][3][3] = {0};
                int case_c = Eshelby_case(axis_t);
                g[G_n].Eshelby_E(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,Gpsets[case_c].Gpalpha);
                //-9

                //-10
                //rotate the eshelby tensor
                //back to the sample axes
                S66 = rotate_C66(voigt(S4_EA), axisb_t);
                rot_4th(R4_EA, axisb_t, R4_SA);
                //-10
            
                //-11
                //Calculate Ctilde
                //refer to Equ[5-33] in Manual 7d
                //M~=(I-S)^-1 * S * M
                //but here we used C~ = (M~)^-1
                //that is 
                // C~ = C * (S^-1 - I)
                //-12 calculate (S)^-1
                S66inv = S66.inverse(); // 6x6 of S^-1
                Matrix6d S66inv_I; // 6x6 of (S^-1 - I)
                //for(int i = 0; i < 6; i++)
                //    for(int j = 0; j < 6; j++)
                //        S66inv_I(i,j) = S66inv(i,j) - Iij(i,j);
                S66inv_I = S66inv - Matrix6d::Identity();
                //-12
                // C~ = C * (S^-1 - I)
                Ctilde = CSC * S66inv_I;
                //-11
            }
            //-7 end of Ishape == 1;


            //-13 store some matrix into grain
            //store the C~
            g[G_n].Update_Metilde_g(Ctilde.inverse());
            //store the PI*(S^-1)
            double S66inv4th[3][3][3][3] = {0};
            Chg_basis(S66inv,S66inv4th);
            mult_4th(R4_SA,S66inv4th,RSinv_SA);
            g[G_n].Update_RSinv_C_g(RSinv_SA);
            //-13

            //-14 
            //Calculate the localization tensor B_g 
            // _g means the value depends on the grain
            //refer to Equ[5-35] in Manual 7d
            // B_g = (M_g + M~)^-1 * (M_ + M~)
            Matrix6d Metilde = Ctilde.inverse();
            Me_g = g[G_n].get_Mij6_J_g();
            Matrix6d Part1 = Me_g + Metilde;
            Matrix6d Part1_inv = Part1.inverse();
            Matrix6d Part2 = SSC + Metilde; 
            Matrix6d B_g = Part1_inv * Part2;
            //-14

            ///////
            //-15
            //Calculate the New elastic consistent stiffness CNEW
            //refer to Equ[5-40a] in Manual 7d
            SSC_new += Me_g * B_g * g[G_n].get_weight_g();
            //-15           
        } //loop over grains

        //cout << "\nThe CNEW:\n" << CNEW << endl;
        //-16
        //error between CSC and CNEW
        RER=Errorcal(SSC,SSC_new);
        SSC = 0.5*(SSC_new+SSC_new.transpose());
        CSC = SSC.inverse();
        //-16
        cout << "**Error in  ESC iteration_" << IT << ":\t" << RER << endl;

    } //while loop

    //SSC = Msup * Btovoigt(CSC.inverse());

    return 0;
}

int polycrystal::Selfconsistent_P(int Istep, double ERRM, int ITMAX)
{
    Matrix5d MNEW; //VP compliance updated in every do-while
    Vector5d D0_new; //the back-extrapolated term updated in every do-while
    Matrix5d B_g_ave; // <B_g> in Equ[5-41a] average of B_g
    Vector5d b_g_ave; // <b_g> in Equ[5-41b] average of b_g

    double DUMMY;

    int IT = 0; //loop flag
    double RER = 2*ERRM;
    while((RER >= ERRM) && (IT < ITMAX))
    {   
        IT++;

        MNEW = Matrix5d::Zero();
        B_g_ave = Matrix5d::Zero();
        D0_new = Vector5d::Zero();
        b_g_ave = Vector5d::Zero();
        //
        Matrix5d Mtilde;
        Matrix5d S55;//the Sijkl Equ[5-33] in sample axes in Manual 7d
        double R4_SA[3][3][3][3]; //PIijkl
        double RSinv_SA[3][3][3][3];
        Matrix5d S55_inv;
        Matrix5d R55;
        Matrix3d axisb_t;
        Vector3d axis_t;
        //

        //solve the eshelby tensor in the common ellipsoid
        if(Ishape == 0)
        {
            //-1
            //rotate the macro VP stiffness (C_VP_SC)
            //from Sample axes to the ellipsoid axes;
            // if Ishape = 0, which means the ellipsoid axes 
            // keep unchanged in grains,
            // the transform matrix should be taken from
            // the polycrystal 
            axisb_t = ell_axisb; 
            axis_t = ell_axis;

            Matrix6d C66 = rotate_C66(Btovoigt(C_VP_SC), axisb_t.transpose());
            //-1
            
            //-3
            //Calculate Eshelby tensor
            double S4_EA[3][3][3][3] = {0}; 
            double R4_EA[3][3][3][3] = {0};
            
            int case_c = Eshelby_case(axis_t);
            g[0].Eshelby_P(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,\
            Gpsets[case_c].Gpalpha,Gpsets[case_c].Gpaww,Gpsets[case_c].Gpww);
            //-3

            //-4
            //rotate the eshelby tensor
            //back to the sample axes
            S55 = voigttoB5(rotate_C66(voigt(S4_EA), axisb_t));
            rot_4th(R4_EA, axisb_t, R4_SA);
            //-4

            //Calculate Mtilde
            //refer to Equ[5-33] in Manual 7d
            //M~=(I-S)^-1 * S * M
            Matrix5d I_S55, I_S55_inv; // (I-S) and (I-S)^-1
            //for(int i = 0; i < 5; i++)
            //    for(int j = 0; j < 5; j++)
            //        I_S55(i,j) = Iij(i,j) - S55(i,j);
            I_S55 = Matrix5d::Identity() - S55;
            I_S55_inv = I_S55.inverse();

            Mtilde = I_S55_inv * S55 * M_VP_SC;
            S55_inv = S55.inverse();
        }
        //-1


        DVP_AV = Vector5d::Zero();

        for(int G_n = 0; G_n < grains_num; G_n++)
        {
            //-2
            //rotate the macro VP stiffness 
            //from Sample axes to the ellipsoid axes;
            // if Ishape = 1, which means the ellipsoid axes varies in grains,
            // the transform matrix should be taken out of
            // each grain  
            if(Ishape == 1)
            {
                axis_t = g[G_n].get_ell_axis_g();
                axisb_t = g[G_n].get_ell_axisb_g();
                Matrix6d C66 = rotate_C66(Btovoigt(C_VP_SC), axisb_t.transpose());
                
                //Calculate Eshelby tensor
                double S4_EA[3][3][3][3] = {0}; 
                double R4_EA[3][3][3][3] = {0};
                int case_c = Eshelby_case(axis_t);
                g[G_n].Eshelby_P(S4_EA,R4_EA,axis_t,C66,Gpsets[case_c].Gpaa6,Gpsets[case_c].Gpaaww6,\
                Gpsets[case_c].Gpalpha,Gpsets[case_c].Gpaww,Gpsets[case_c].Gpww);
                //-3

                //-4
                //rotate the eshelby tensor
                //back to the sample axes
                S55 = voigttoB5(rotate_C66(voigt(S4_EA), axisb_t));
                rot_4th(R4_EA, axisb_t, R4_SA);
                //-4

                //-4
                //Calculate Mtilde
                //refer to Equ[5-33] in Manual 7d
                //M~=(I-S)^-1 * S * M
                Matrix5d I_S55, I_S55_inv; // (I-S) and (I-S)^-1
                //for(int i = 0; i < 5; i++)
                //    for(int j = 0; j < 5; j++)
                //        I_S55(i,j) = Iij(i,j) - S55(i,j);
                I_S55 = Matrix5d::Identity() - S55;
                I_S55_inv = I_S55.inverse();

                Mtilde = I_S55_inv * S55 * M_VP_SC;
                S55_inv = S55.inverse();
            }
            //-2

            //-5 store some matrix into grain
            //store the M~
            g[G_n].Update_Mptilde_g(Mtilde); 
            //only affine case (interaction = 1) 
            //store the PI*(S^-1)
            double S66inv4th[3][3][3][3] = {0};
            Chg_basis(S55.inverse(),S66inv4th);
            mult_4th(R4_SA,S66inv4th,RSinv_SA);
            g[G_n].Update_RSinv_VP_g(RSinv_SA);
            //06.02
            //-5

            //07.14
            //Matrix5d M_g = voigttoB5(g[G_n].get_Mpij6_g());
            //Vector5d d0_g = Chg_basis5(g[G_n].get_d0_g());
            Matrix5d M_g = g[G_n].get_Mpij6_g();
            Vector5d d0_g = g[G_n].get_d0_g(); 

            double wei = g[G_n].get_weight_g();

            //-6 
            //Calculate the localization tensor B_g 
            // _g means the value depends on the grain
            //refer to Equ[5-35] in Manual 7d
            // B_g = (M_g + M~)^-1 * (M_ + M~)
            Matrix5d Part1 = M_g + Mtilde;
            Matrix5d Part1_inv = Part1.inverse();
            Matrix5d Part2 = M_VP_SC + Mtilde; 
            Matrix5d B_g = Part1_inv * Part2;
            ///
            //refer to Equ[5-35] in Manual 7d
            // b_g = (M_g + M~)^-1 * (d- - d-_g)
            //Matrix3d dtemp =  voigt(D0) - d0_g;//(d- - d-_g)
            //Vector6d b_gv =  voigt5to6(Part1_inv) * Msupinv * voigt(dtemp);
            Vector5d b_gv =  Part1_inv * (Chg_basis5(voigt(D0))-d0_g);
            //Matrix3d b_g =  voigt(b_gv);
            //-6
            
            // MNEW = <M_g * B_g>
            // Equ[5-40a]
            MNEW += M_g * B_g * wei;
            //<B_g>
            B_g_ave += B_g * wei;
            // Equ[5-40b]
            // < M_g * b_g + d0_g>
            //Vector6d Mr_br =  voigt5to6(M_g) * Msupinv * b_gv; // M_g * b_g
            Vector5d Mr_br =  M_g * b_gv; 
            D0_new += ( d0_g +  Mr_br ) * wei; // < M_g * b_g + d0_g>
            // <b_g>
            b_g_ave += b_gv * wei;

            //06.06
            DVP_AV +=  Chg_basis5(g[G_n].get_Dijp_g()) * wei;
            //06.06

        } //loop over grains

    // Equ[5-41a]
    // <M_g * B_g> * <B_g>^-1
    MNEW = MNEW * B_g_ave.inverse();
    Matrix5d MNEW2 = 0.5*(MNEW+MNEW.transpose());
    MNEW = MNEW2;
    // Equ[5-41b]
    // < M_g * b_g + d0_g> - <M_g * B_g> * <B_g>^-1 * <b_g>
    //<M_g * B_g> * <B_g>^-1 * <b_g>
    //Vector6d M_bg =  voigt5to6(MNEW) * Msupinv  * b_g_ave;
    Vector5d M_bg =  MNEW * b_g_ave;
    D0_new = D0_new - M_bg;

    //calculate the error between 
    //the input M_VP_SC of do-while loop
    //and the output MNEW of the loop 
    RER = Errorcal(M_VP_SC, MNEW);
    M_VP_SC = MNEW;
    D0 = voigt(Chg_basis(D0_new));

    C_VP_SC = M_VP_SC.inverse();
  
    cout << "**Error in VPSC iteration_" << IT << ":\t" << RER << endl;
    }//while loop
    return 0;
}


int polycrystal::Update_Fij(double Tincr)
{
    Matrix3d Fnew;
    Fnew = Matrix3d::Zero();
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
    {
        Fnew(i,j) += (Tincr*Udot_AV(i,k)+Iij(i,k))*Fij_m(k,j);
    }
    Fij_m = Fnew;
    return 0;
}

int polycrystal::Update_shape()
{
    //-1 F.transpose()*F
    //BX = Matrix3d::Zero();
    //for(int i = 0; i < 3; i++)
    //    for(int j = 0; j < 3; j++)
    //        BX(i,j) = Fij_g.row(i)*Fij_g.col(j);
    Matrix3d BX;
    Matrix3d FT = Fij_m.transpose();
    BX =  Fij_m*FT;
    //-1 F.transpose()*F

    //-2 solve the eigen vector of BX
    //and sort the value from largest to smallest
    //EigenSolver<Matrix3d> es(BX);
    Matrix3d BX_vectors;
    Vector3d BX_value;
    Jacobi(BX,BX_value,BX_vectors);
    //BX_value = (es.eigenvalues()).real();
    //BX_vectors = (es.eigenvectors()).real();

    Eigsrt(BX_vectors, BX_value);

    //-2 solve the eigen vector of BX
    //and sort the value from largest to smallest
    
    Matrix3d B = BX_vectors;
    Vector3d W = BX_value;
    //-3 
    //redefine Axis(1) (the second) to be the largest  
    //to improve the accuracy in calculation of Eshelby tensor
    //IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
    //HANDED BY EXCHANGING 1 AND 2.
    double Sign = -1;
    double temp;
    if(B.determinant() <= 0) Sign = 1;
    for(int i = 0; i < 3; i++)
    {
        temp = B(i,0);
        B(i,0) = B(i,1);
        B(i,1) = temp * Sign;
    }
    temp = W(0); W(0) = W(1); W(1) = temp;
    //-3 

    //-4 update the stretching of ellipsoid
    double Ratmax=sqrt(W(1)/W(2));
    double Ratmin=sqrt(W(0)/W(2));

    ell_axisb = B;
    if(!Iflat) //if Iflat = 0
        for(int i = 0; i < 3; i++)
            ell_axis(i) = sqrt(W(i));
    //if Iflat = 1, the axis of ellipsoid keeps unchange
    //-4 update the stretching of ellipsoid
    
    //-5 update the ellipsoid orientation
    Matrix3d BT = B.transpose();
    ellip_ang = Euler_trans(BT);
    
    //-5

    //-6 Update the Iflat_g according to the Max axes ratio of ellipsoid
    if((!Iflat)&&(Ratmax >= ell_crit_shape))
    {
        Iflat = 1;
    }
    //-6

    //-7 Iflat_g = 1; recaculates the Fij_g in grain
    else if((Iflat)&&(Ratmax >= ell_crit_shape))
    {
        W(1) = W(1)/4;
        if(Ratmin >= 0.5 * ell_crit_shape)
            W(0) = W(0)/4;
        for(int i = 0; i < 3; i++)
            ell_axis(i) = sqrt(W(i));
        
        Matrix3d Fijx;
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                Fijx(i,j) = Iij(i,j) * ell_axis(i);
        
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
        {
            Fij_m(i,j) = 0;
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
                    Fij_m(i,j) = Fij_m(i,j) + B(i,k)*B(j,l)*Fijx(k,l);
        }
    }
    return 0;
}


int polycrystal::EVPSC(int istep, double Tincr,\
 bool Iupdate_ori,bool Iupdate_shp,bool Iupdate_CRSS)
{   
    double errd, errs, err_g;
    cout << "********STEP********\n\t" \
        << istep << endl << "********STEP********\n";

    Sig_m_old = Sig_m; //save the stress of the last step
    for(int G_n = 0; G_n < grains_num; ++G_n)
        g[G_n].save_sig_g_old();//save the grain stress of the last step

    for(int i = 0; i < 30; ++i)
    {
        //save the input for error calculation
        Sig_in = Sig_m;
        Dij_in = Dij_m;
        sig_in_AV = Matrix3d::Zero();

        for(int G_n = 0; G_n < grains_num; ++G_n) 
            sig_in_AV += g[G_n].get_stress_g() * g[G_n].get_weight_g();

        ///////////
        cout << "iteration: " << i+1 << endl;
        Selfconsistent_E(istep, SC_err_m, SC_iter_m);
        Selfconsistent_P(istep, SC_err_m, SC_iter_m);
        Cal_Sig_m(Tincr); 
        Cal_Sig_g(Tincr);
        Update_AV();
        ///////////

        errs = Errorcal(Sig_m, Sig_in);
        errd = Errorcal(Dij_m, Dij_in);
        err_g = Errorcal(Sig_AV, sig_in_AV);

        cout << "\nerr_s:\t" << errs << endl;
        cout << "err_d:\t" << errd << endl;
        cout << "err_g:\t" << err_g << endl;
        cout << "=-=-=-=-=-=-=-=-=-=\n"; 
        if((errs<errS_m)&&(errd<errD_m)&&(err_g<err_g_AV)) break;
    }

    Eps_m += Dij_m * Tincr; //update the macro strain tensor

    //update the shape of ellipsoid
    if(Ishape == 0)
    {
        Update_Fij(Tincr);
        if(Iupdate_shp) Update_shape();
    }

    //update the state in deformation systems and 
    // crystalline orientation 
     for(int G_n = 0; G_n < grains_num; ++G_n)
    {
        g[G_n].Update_shear_strain(Tincr);
        if(Iupdate_ori) g[G_n].Update_orientation(Tincr, Wij_m, Dije_AV, Dijp_AV);
        if(Iupdate_CRSS) g[G_n].Update_CRSS(Tincr);
        if(Ishape == 1)
        {
            g[G_n].Update_Fij_g(Tincr);
            if(Iupdate_shp) g[G_n].Update_shape_g();
        }
    }
    return 0;
}


void polycrystal::Cal_Sig_m(double Tincr)
{
    Wij_m = Matrix3d::Zero();
    //why not Wij = Udot - Dij ?
    //because the Udot need be update
    //some components are not imposed
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    {
        if(IUdot(i,j)==1 && IUdot(j,i)==1)
           Wij_m(i,j) = 0.5*(Udot_m(i,j) - Udot_m(j,i));
        else if(IUdot(i,j)==1)
            {
                Wij_m(i,j) = Udot_m(i,j) - Dij_m(i,j);
                Wij_m(j,i) = -Wij_m(i,j);
            }
        Udot_m(i,j) = Dij_m(i,j) + Wij_m(i,j);
    }


    //-1
    //calculate the Jaumann rate
    Matrix3d Sig_J = Wij_m * Sig_m - Sig_m * Wij_m;
    //-1

    //-2
    //calulate the De = D - Dp - D0
    Vector6d De = Bbasisadd(Chg_basis6(Dij_m), -DVP_AV);
    //-2

    //-3
    //Calculate AX = B
    //X is the macro stress
    Matrix6d A = SSC/Tincr;
    //B
    Matrix3d Mtemp = Sig_m_old / Tincr;

    Vector6d B = De + SSC * ( Chg_basis6(Mtemp) + Chg_basis6(Sig_J));
    //-3

    //-4
    //According to the IUdot and ISDOT to solve AX = B
    //transform to solve At Xt = Bt
    // Xt in the unkown set of Udot and Sig
    Vector6i Is = ISdot;
    Vector6i Id = IDdot;

    Vector6d profac;
    profac << 1,1,1,2,2,2;

    Vector6d BC_D = voigt(Chg_basis(B));
    Vector6d BC_S = voigt(Sig_m);
    Matrix6d AUX2 = Btovoigt(A);

    Vector6d AUX11;
    Matrix6d AUX21;

    for(int i = 1; i <= 6; i++){
        AUX11(i-1) = -Id(i-1)*BC_D(i-1);
        for(int j = 1; j <= 6; j++){
                AUX11(i-1) = AUX11(i-1) + AUX2(i-1,j-1)*Is(j-1)*BC_S(j-1)*profac(j-1);
                AUX21(i-1,j-1) = Is(j-1)*(i/j)*(j/i) - Id(j-1)*AUX2(i-1,j-1)*profac(j-1);
            }
    }

    Vector6d AUX6 = AUX21.inverse()*AUX11;

    //-4

    Vector6d Sig_x, B_x;
    B_x = mult_dot(BC_D,Id) + mult_dot(AUX6,Is);
    Sig_x = mult_dot(BC_S,Is) + mult_dot(AUX6,Id);

    De = Chg_basis6(voigt(B_x)) - SSC * ( Chg_basis6(Mtemp) + Chg_basis6(Sig_J));

    //calulate the D = De + Dp + D0
    Vector6d Dij_m_v = Bbasisadd(De, DVP_AV);

    Dij_m = Chg_basis(Dij_m_v);
    Sig_m = voigt(Sig_x);

}

void polycrystal::Cal_Sig_g(double Tincr)
{
    #pragma omp parallel for num_threads(Mtr)
    for(int G_n = 0; G_n < grains_num; G_n++)
        g[G_n].grain_stress(Tincr, Wij_m, Dij_m, Dije_AV, Dijp_AV, Sig_m, Sig_m_old);
}


void polycrystal::Update_AV()
{
    Dije_AV = Matrix3d::Zero();
    Dijp_AV = Matrix3d::Zero();
    Sig_AV = Matrix3d::Zero();
    Udot_AV = Matrix3d::Zero();
    double wei;
    for(int G_n = 0; G_n < grains_num; G_n++)
    {
        wei = g[G_n].get_weight_g();
        Dije_AV += g[G_n].get_Dije_g() * wei;
        Dijp_AV += g[G_n].get_Dijp_g() * wei;
        Sig_AV += g[G_n].get_stress_g() * wei;
        Udot_AV += g[G_n].get_Udot_g() * wei;
    }
    Dij_AV = Dije_AV + Dijp_AV;
}

Vector6d polycrystal::get_Sig_m(){return voigt(Sig_m);}
Vector6d polycrystal::get_Eps_m(){return voigt(Eps_m);}
void polycrystal::get_euler(fstream &texfile)
{
    IOFormat Outformat(StreamPrecision);
    for(int i = 0; i < grains_num; i++)
    {
        texfile << setprecision(4) << scientific << g[i].get_euler_g().transpose().format(Outformat) << "  ";
        texfile << setprecision(4) << scientific << g[i].get_weight_g() << endl;
    }
}

Vector3d polycrystal::get_ell_axis(){return ell_axis;}
Vector3d polycrystal::get_ellip_ang(){return ellip_ang;}