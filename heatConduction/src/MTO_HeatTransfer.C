//Author: Yu Minghao    Updated: May 2020

static char help[] = "topology optimization of heat conduction problem\n";
#include "fvCFD.H"
#include "simpleControl.H"
#include "MMA/MMA.h"
#include <diff.c>
int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "readThermalProperties.H" 
    #include "opt_initialization.H"
    while (simple.loop(runTime))
    {
        #include "update.H"
        #include "conduction.H"
        #include "costfunction.H"              
        #include "sensitivity.H"
    }
    #include "finalize.H"
    return 0;
}
