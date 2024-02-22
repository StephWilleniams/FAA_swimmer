/* 
Title: Three body microswimmer w/ passive particles
Author: Stephen Williams
*/

// Private header files
#include "functions/header.hpp"

int main(int argc, char* argv[])
{
    /* USER INPUTS */
    double dR = DT*20/100.0;
    double sigR = sqrt(dR);
    int pol = -1;
    
    double KickFreq = 0.06667; // Kick poisson-process frequency.
    double kickFreq = KickFreq*DT; // ND kick-frequency.

    double KickStr = 0.1 ; // Kick rotational 'speed' (Rad/s).
    double kickStr = KickStr*DT; // ND kick rotation per frame.
    /* USER INPUTS */

    int run_ind = atoi(argv[1]);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    cout << "Now starting, running " << T << "s." << endl;

    // Active particle outputs
    string filename_active = "outputs_fine/" + string(argv[1]) + string("_outputs/output_active.txt"); // Output file name.
    ofstream outfile_active = initialise_file(filename_active); // Initialise output.
    // Passive particle outputs
    string filename_passive = "outputs_fine/" + string(argv[1]) + string("_outputs/output_passive.txt"); // Output file name.
    ofstream outfile_passive = initialise_file(filename_passive); // Initialise output.

    // Initialise random seed.
    mt19937 gen = generate_seed(run_ind);
    //mt19937 gen = generate_seed(run_ind);
    set_active_geometry(rA,pol,rSeg,sDV);
    //zero_IVP(ns,x); // All particles at (x,y)=(0,0).
    //partpair_tester_IVP(ns,xs);
    //partpair_tester_IVP_2(na,xa,np,xp);
    no_OL_IVP(na, xa, np, xp, xL, xR, yB, yT, rA, rP, gen); // Uniform distribution, non-overlapping.
    output_pos(-1,na,xa,np,xp,RA,DT,outfile_active,outfile_passive); // Output active particle positions.

    for(int t=0; t<nSteps; t++){
        update_pos_RKII(na,xa,rA,rSeg,sDV,np,xp,rP,vA,sigR,sigT,xL,xR,yB,yT,f0,fR,fricPar,fricPerp,kickFreq,kickStr,dT,DT,gen); // Increment particle positions.
        if((t % 1)==0) {output_pos(t,na,xa,np,xp,RA,DT,outfile_active,outfile_passive);}; // Output particle positions.
        if((t % int(T/(DT*100)))==0) {cout << 100*double(t)/double(nSteps) << " percent complete." << endl;};  // Update progress.
    }

    outfile_active.close(); // Finalise output.
    outfile_passive.close(); // Finalise output.

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0; // End script.

} 