#ifndef PROCESSES_H
#define PROCESSES_H

#include <string>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

#include "Grains.h"
#include "Polycrystals.h"
#include "Toolbox.h"

namespace Procs{

class Process
{
    private:
        Matrix3d Udot_input; //the velocity gradient in process file
        Matrix3d Ddot_input; //the strain rate calculated by the velocity gradient in process file
        Matrix3d Sdot_input; //the stress tensor in process file

        Matrix3i IUdot; //the flag of known (control by Udot_input) and unknown (calculated by EVPSC) velocity components
        Vector6i IDdot; //the flag of known and unknown strain rate components
        Vector6i ISdot; //the flag of known and unknown stress componets

        double Eincr; //the increment of strain in every istep
        int Ictrl;
        int Nsteps; // total steps
        double Temp; // /(K^-1) temperature 
        double Tincr;

        //Output files
        fstream ss_out; //output of the macro stress-strain curves
        fstream tex_out; //output of the texture

        int texctrl; //print the texture every n steps(0 means only print at the end)

        //update
        bool Iupdate_ori; //update the orientation 1:yes 0:no
        bool Iupdate_shp; //update the ellipsoid shape 1:yes 0:no
        bool Iupdate_CRSS; //update the CRSS 1:yes 0:no

    public:
        Process();
        ~Process();

        //get the total steps and increment of a process from files
        //Vector4d 0: Nsteps; 1: Ictrl; 2: Eincr; 3: Temperature;
        void load_ctrl(Vector4d);

        //
        void Update_ctrl(Vector3i);

        //get the velocity gradient in a process file
        void get_Udot(Matrix3d);

        //get the Cauchy stress tensor in a process file
        void get_Sdot(Matrix3d);

        //get the known and unknown flag tensor in a process file
        void get_IUdot(Matrix3i);
        void get_ISdot(Vector6i);

        void loading(Polycs::polycrystal &);

        /////
        //Output functions:
        //output of stress&strain curves  
        void Out_sscurves(Polycs::polycrystal &);

        //output of texture
        void Out_texture(Polycs::polycrystal &, int);
        void Out_texset(int);

};

}

#endif