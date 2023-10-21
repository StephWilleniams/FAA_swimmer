/* 
Title: Three body microswimmer.
*/

// Private header files
#include "header.hpp"
#include "constants.hpp"
#include "rng.hpp"
#include "functions.hpp"

int main() 
{
    string filename = "output.txt"; // Output file name.
    ofstream outfile = initialise_file(filename); // Initialise output.

    //zero_IVP(ns,x); // All particles at (x,y)=(0,0).
    partpair_tester_IVP(ns,x);
    //no_OL_IVP(ns, x, xL, xR, yB, yT, rS, gen); // Uniform distribution, non-overlapping.
    output_pos(-1,ns,x,RS,DT,outfile); // Output particle positions.

    for(int t=0; t<nSteps; t++){
        update_pos_RKII(ns,x,vS,sigR,xL,xR,yB,yT,dT,DT,gen); // Increment particle positions.
        output_pos(t,ns,x,RS,DT,outfile); // Output particle positions.
    }

    outfile.close(); // Finalise output.
    return 0; // End script.

} 