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

int main(int argc, char* argv[])
{
    /* USER INPUTS */
    double dR = DT*double(atoi(argv[1]))/100;
    double sigR = sqrt(dR);
    int pol = atoi(argv[2]);
    double KickFreq = double(atoi(argv[3]))/100; // Kick poisson-process frequency.
    double kickFreq = KickFreq*DT; // ND kick-frequency.
    double KickStr = double(atoi(argv[4]))/100; // Kick rotational 'speed'.
    double kickStr = KickStr*DT; // ND kick rotation per frame.
    /* USER INPUTS */

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    cout << "Now starting, running " << T << "s." << endl;

    // Active particle outputs
    string filename_active =  string(argv[5]) + "_outputs/output_active.txt"; // Output file name.
    ofstream outfile_active = initialise_file(filename_active); // Initialise output.
    // Passive particle outputs
    string filename_passive = string(argv[5]) + "_outputs/output_passive.txt"; // Output file name.
    ofstream outfile_passive = initialise_file(filename_passive); // Initialise output.

    // Initialise random seed.
    mt19937 gen = generate_seed(atoi(argv[5]));
    //zero_IVP(ns,x); // All particles at (x,y)=(0,0).
    //partpair_tester_IVP(ns,xs);
    set_active_geometry(pol,rSeg,sDV);
    no_OL_IVP(na, xa, np, xp, xL, xR, yB, yT, rA, rP, gen); // Uniform distribution, non-overlapping.
    output_pos(-1,na,xa,np,xp,RA,DT,outfile_active,outfile_passive); // Output active particle positions.

    for(int t=0; t<nSteps; t++){
        update_pos_RKII(na,xa,rA,rSeg,sDV,np,xp,rP,vA,sigR,sigT,xL,xR,yB,yT,f0,fR,fricPar,fricPerp,kickFreq,kickStr,dT,DT,gen); // Increment particle positions.
        if((t % 50)==0) {output_pos(t,na,xa,np,xp,RA,DT,outfile_active,outfile_passive);}; // Output particle positions.
        if((t % int(T/(DT*100)))==0) {cout << 100*double(t)/double(nSteps) << " percent complete." << endl;};  // Update progress.
    }

    outfile_active.close(); // Finalise output.
    outfile_passive.close(); // Finalise output.

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0; // End script.

} 