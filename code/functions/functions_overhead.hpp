#ifndef functions_overhead
#define functions_overhead

#include "header.hpp"
#include "functions_rng.hpp"

ofstream initialise_file(string inputString){

    // Create an ofstream object to work with files
    ofstream outputFile;

    // Open a file for writing (this will create or overwrite the file)
    outputFile.open(inputString);

    if (!outputFile.is_open()) {
        cerr << "Failed to open the file for writing." << endl;
        exit(1);
    }

    return outputFile;
}

void output_pos(double t, int na, double xa[][4], int nc, double xc[][3], double NDL, double NDT, ofstream& outfile_active, ofstream& outfile_passive){
    for (int n=0; n<na; n++){ // Loop over particle index.
        outfile_active << t << ' ' << n+1;
        for (int i=0; i<4; i++){ // Loop over dimension.
            outfile_active << ' ' << xa[n][i];
        }
        outfile_active << endl;
    }

    for (int n=0; n<nc; n++){ // Loop over particle index.
        outfile_passive << t << ' ' << n+1;
        for (int i=0; i<2; i++){ // Loop over dimension.
            outfile_passive << ' ' << xc[n][i];
        }
        outfile_passive << endl;
    }

    return;
}

void output_parameters(){
    return;
}

#endif // 