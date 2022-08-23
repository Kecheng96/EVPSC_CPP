#include <time.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <Eigen/Dense>

#include "Input.h"
#include "Polycrystals.h"

int main()
{   
    //declare object
    Polycs::polycrystal metal; 

    string ftex, fsx, fload;
    if(EVPSCinput(ftex, fsx, fload)) exit(0);
    //need a new .in file
    if(texinput(ftex, metal)) exit(0);
    if(sxinput(fsx, metal)) exit(0);
    if(loadinput(fload, metal)) exit(0);
    
    double Start = clock();
    metal.EVPSC();
    double End = clock();
    cout << "The run time is: " <<(double)(End - Start) / CLOCKS_PER_SEC << " sec" << std::endl;
    metal.check_euler(0);
    metal.check_euler(1);
    return 0;
}