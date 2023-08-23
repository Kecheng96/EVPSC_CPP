#include <time.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <Eigen/Dense>

#include "Input.h"
#include "Processes.h"
#include "Polycrystals.h"

int main()
{
    Polycs::polycrystal metal; //declare object
    Procs::Process Proc1;

    string ftex, fsx, fload;
    if(EVPSCinput(ftex, fsx, fload, Proc1)) exit(0);
    if(texinput(ftex, metal)) exit(0);
    if(sxinput(fsx, metal)) exit(0);
    if(loadinput(fload, Proc1)) exit(0);
    
    double Start = clock();
    Proc1.loading(metal);
    double End = clock();
    cout << "The run time is: " <<(double)(End - Start) / CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}