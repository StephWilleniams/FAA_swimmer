/* 
Title: Three body microswimmer w/ passive particles
Author: Stephen Williams
*/

// Private header files
#include "header.hpp"
#include "constants.hpp"
#include "functions.hpp"
#include "functions_rng.hpp"
#include "functions_IVP.hpp"
#include "functions_overhead.hpp"

int main()
{

    // Active particle outputs
    string filename_active = "output_active.txt"; // Output file name.
    ofstream outfile_active = initialise_file(filename_active); // Initialise output.

    // Passive particle outputs
    string filename_passive = "output_passive.txt"; // Output file name.
    ofstream outfile_passive = initialise_file(filename_passive); // Initialise output.

    //zero_IVP(ns,x); // All particles at (x,y)=(0,0).
    partpair_tester_IVP(ns,xs);
    //no_OL_IVP(ns, xs, np, xp, xL, xR, yB, yT, rA, rP, gen); // Uniform distribution, non-overlapping.
    output_pos(-1,ns,xs,np,xp,RA,DT,outfile_active,outfile_passive); // Output active particle positions.

    for(int t=0; t<nSteps; t++){
        update_pos_RKII(ns,xs,rA,rSeg,sDV,polarization,np,xp,vA,sigR,sigT,xL,xR,yB,yT,dT,DT,gen); // Increment particle positions.
        if((t % 1)==0) {output_pos(t,ns,xs,np,xp,RA,DT,outfile_active,outfile_passive);}; // Output particle positions.
    }

    outfile_active.close(); // Finalise output.
    outfile_passive.close(); // Finalise output.
    return 0; // End script.

} 