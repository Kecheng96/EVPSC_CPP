#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

#include "Polycrystals.h"
#include "Processes.h"

int EVPSCinput(string &,string &,string &, Procs::Process &); //read .in file
int sxinput(string, Polycs::polycrystal &); //read .sx file
int texinput(string, Polycs::polycrystal &); //read .tex file
int loadinput(string, Procs::Process &Proc);
VectorXd getnum(string, int);

int EVPSCinput(string &ftex,string &fsx,string &fload, Procs::Process &Proc)
{
    fstream ininp;
    ininp.open("EVPSC_CPP.in",ios::in); //open EVPSC.in
    if (ininp.is_open())
    {
        //read the file path
        string tp;
        getline(ininp, tp); //skip
        getline(ininp, ftex);
        getline(ininp, tp); //skip
        getline(ininp, fsx); 
        getline(ininp, tp); //skip
        getline(ininp, fload); 

        //read the update control
        getline(ininp, tp); //skip
        getline(ininp, tp); //skip
        getline(ininp, tp); //skip
        getline(ininp, tp); 
        VectorXd temp1 = getnum(tp, 3);
        Vector3i temp2;
        for(int i=0; i<3; i++)
            temp2(i) = int(temp1(i));
        Proc.Update_ctrl(temp2);

        //read output control
        getline(ininp, tp); //skip
        getline(ininp, tp); //skip
        getline(ininp, tp); //skip
        getline(ininp, tp); 
        VectorXd temp = getnum(tp, 1);
        Proc.Out_texset(int(temp(0)));

        ininp.close(); 
        return 0;
    }
    else
    {
        cout << "Error code 0: loading file cannot be opened.\n";
        return 1;
    }
}

int loadinput(string fname, Procs::Process &Proc)
{
    fstream loadinp;
    loadinp.open(fname,ios::in); //open load

    if (loadinp.is_open())
    {   //checking whether the file is open
        string tp;

        //1st line is the loading control option
        getline(loadinp, tp);
        Vector4d Victrl = getnum(tp, 4);
        Proc.load_ctrl(Victrl);
         
        getline(loadinp, tp);//skip one line  
        //boundary condition
        Matrix3i IUdot;
        for(int i = 0; i < 3; i++)
        {
            getline(loadinp, tp);
            Vector3d temp = getnum(tp, 3);
            for(int j = 0; j < 3; j++)
                IUdot(i,j) = int(temp(j));
        }
        Proc.get_IUdot(IUdot);

        getline(loadinp, tp);//skip one line  
        //boundary condition
        Matrix3d Udot;
        for(int i = 0; i < 3; i++)
        {
            getline(loadinp, tp);
            Udot.row(i) = getnum(tp, 3);
        }
        Proc.get_Udot(Udot);

        getline(loadinp, tp);//skip one line  
        //boundary condition
        Vector6i ISdot;
        getline(loadinp, tp);
        VectorXd temp = getnum(tp, 3);
        ISdot(0)=int(temp(0));ISdot(5)=int(temp(1));ISdot(4)=int(temp(2));
        getline(loadinp, tp);
        temp = getnum(tp, 2);
        ISdot(1)=int(temp(0));ISdot(3)=int(temp(1));
        getline(loadinp, tp);
        temp = getnum(tp, 1);
        ISdot(2)=int(temp(0));

        Proc.get_ISdot(ISdot);

        getline(loadinp, tp);//skip one line  
        //boundary condition
        Vector6d Sig_m;
        getline(loadinp, tp);
        temp = getnum(tp, 3);
        Sig_m(0)=temp(0);Sig_m(5)=temp(1);Sig_m(4)=temp(2);
        getline(loadinp, tp);
        temp = getnum(tp, 2);
        Sig_m(1)=temp(0);Sig_m(3)=temp(1);
        getline(loadinp, tp);
        temp = getnum(tp, 1);
        Sig_m(2)=temp(0);

        Proc.get_Sdot(voigt(Sig_m));

        loadinp.close(); //close the file object.
        return 0;        
    }
    else
    {
        cout << "Error code 0: loading file cannot be opened.\n";
        return 1;
    }
}

int sxinput(string fname, Polycs::polycrystal &pcrys)
{
    fstream sxinp;
    sxinp.open(fname,ios::in); //open .sx

    if (sxinp.is_open())
    {   //checking whether the file is open
        string tp;
        //skip first line; 
        getline(sxinp, tp);
        //second line: crystal symmetry (crysym)
        string crysym;
        getline(sxinp, crysym);
        crysym = crysym.substr(0,5);
        //cout << "crystal symmetry:\n" << crysym << "\n"; 
        //third line: crystal constant (cdim,cang)
        getline(sxinp, tp);
        VectorXd ccon = getnum(tp, 6);
        //cout << "crystal constant:\n" <<ccon.transpose() << "\n";
        //
        pcrys.ini_cry(crysym, ccon);
        //
        int Millern = pcrys.get_Millern();

        //Elastic constant
        getline(sxinp, tp);  //skip a line;
        MatrixXd Cij6(6,6);
        for (int i=0; i < 6; i++)
        {
            getline(sxinp, tp);
            Cij6.row(i) = getnum(tp, 6);
        }
        //cout << "elastic constant:\n" << Cij6 << "\n"; 
        pcrys.ini_Cij6(Cij6);

        //Thermal coefficients
        getline(sxinp, tp);  //skip a line;
        getline(sxinp, tp);
        VectorXd ther = getnum(tp, 6);
        //cout << "Thermal coefficients:\n" << ther.transpose() << "\n";
        //
        pcrys.ini_therm(ther);
        //

        //slip and twinning modes
        getline(sxinp, tp);  //skip a line;
        getline(sxinp, tp);
        VectorXd nmodesx = getnum(tp, 1);
            //total number of the modes in crystal(nmodesx)
        getline(sxinp, tp);
        VectorXd nmodes = getnum(tp, 1);
            //considered in current simmulation(nmodes)
        getline(sxinp, tp);
        VectorXd mode_i = getnum(tp, int(nmodes(0)));
            //the index of modes(mode_i)
        //cout << "mode_i:\n" << mode_i.transpose() << "\n";
        //
        pcrys.ini_gmode(int(nmodes(0)));
        //
        //pcrys.check_gmode();

        //the normal and direction of each system
        int modeflag = 0;
        bool f = 1;
        for(int i = 0; i < nmodesx(0); i++)
        {
            getline(sxinp, tp);  //skip a line;
            getline(sxinp, tp);
            VectorXd mode_info = getnum(tp, 4);
            //mode_info 0: the serial number
            //1: number of mechnical systems
            //2: flag of slip (0 for twin; 1 for slip)
            //3: flag of twin (1 for twin; 0 for slip)
            //cout << int(mode_info(1)) << endl;
            MatrixXd nor_dir(int(mode_info(1)),2*Millern);
            //normal and direction of slip plane

            //special for twin
            if(int(mode_info(3))) getline(sxinp, tp);

            for (int i = 0; i < int(mode_info(1)); i++)
            {            
                getline(sxinp, tp);
                nor_dir.row(i) = getnum(tp, 2*Millern);    
            }
            //cout << "normal and direction:\n" << nor_dir << "\n"; 
            //
            f = (mode_i.array() == i+1).any();
            if(f)
            {
                pcrys.ini_sn(nor_dir, int(mode_info(2)), int(mode_info(1)), modeflag);
                modeflag++;
            }                
            //
        }
        //pcrys.check_sn();

        getline(sxinp, tp);  //skip a line;
        //hardening law selection
        getline(sxinp, tp); 
        VectorXd Hlaw = getnum(tp, 1); 
        //flag of rate sensitive model (1: yes; 0: no)
        getline(sxinp, tp); 
        VectorXd Irate = getnum(tp, 1);
        //grain siez: um
        getline(sxinp, tp); 
        VectorXd GZ = getnum(tp, 1);
        //
        pcrys.ini_GZ(GZ(0));
        //
        //hardening parameters of modes
        VectorXd nrsx, CRSS_p, hst;
        for(int i = 0; i < int(nmodes(0)); i++)
        {
            getline(sxinp, tp);  //skip a line;
            //cout << tp << "\n";
            //rate sensitive
            getline(sxinp, tp);
            nrsx = getnum(tp, 1);
            //CRSS parameters
            getline(sxinp, tp);
            CRSS_p = getnum(tp, 4);
            //hst
            getline(sxinp, tp);
            hst = getnum(tp, int(nmodes(0)));
            //cout << "hst:\n" << hst.transpose() << "\n";
            pcrys.ini_hardening(nrsx(0), CRSS_p, hst, i);
        }
        sxinp.close(); //close the file object.
        return 0;
    }
    else
    {
        cout << "Error code 0: .sx file cannot be opened\n";
        return 1;
    }
}

int texinput(string fname, Polycs::polycrystal &pcrys)
{
    fstream texinp;
    texinp.open(fname,ios::in); //open .tex
    if (texinp.is_open())
    {   //checking whether the file is open
        string tp;
        //skip 3 lines;
        for(int i = 0; i < 3; i++)
        {
            getline(texinp, tp);
        }
        //number of grains 
        getline(texinp, tp);
        VectorXd Gn = getnum(tp, 1);
        //
        pcrys.grains_n(int(Gn(0)));
        // 
        //Euler angle and weighs
        Vector4d Euler_w(0,0,0,0);
        for(int i = 0; i < int(Gn(0)); i++)
        {
            getline(texinp, tp);
            Euler_w = getnum(tp, 4);
            pcrys.ini_euler(getnum(tp, 4),i);
        }
        texinp.close(); //close the file object.
        pcrys.Norm_weight();
        return 0;        
    }
    else
    {
        cout << "Error code 0: texture file cannot be opened.\n";
        return 1;
    }
}

VectorXd getnum(string strin, int num)
{
    int i = 0;
    VectorXd Vtemp(num);

    //string pattern("\\d+(\\.\\d+)?");
    string pattern("[+-]?[\\d]+([\\.][\\d]*)?([Ee][+-]?[\\d]+)?");
	regex r(pattern);
	smatch results;

	string::const_iterator iter_begin = strin.cbegin();
	string::const_iterator iter_end = strin.cend();
	while (regex_search(iter_begin, iter_end,  results,  r))
	{
        if (i >= num) break;
		Vtemp(i)=stof(results[0].str());
		iter_begin = results[0].second;	
        i++;
	}
    return Vtemp;
}

#endif