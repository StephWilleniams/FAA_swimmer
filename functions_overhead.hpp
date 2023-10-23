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

void output_pos(double t, int ns, double x[][3], double NDL, double NDT, ofstream& outfile){
    for (int n=0; n<ns; n++){ // Loop over particle index.
        outfile << t << ' ' << n+1;
        for (int i=0; i<3; i++){ // Loop over dimension.
            outfile << ' ' << x[n][i];
        }
        outfile << endl;
    }
    return;
}

#endif // 