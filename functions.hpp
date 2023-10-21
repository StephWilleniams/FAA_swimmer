#ifndef functions
#define functions

#include "header.hpp"
#include "rng.hpp"

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

void zero_IVP(int ns, double x[][3]){
    for (int n=0; n<ns; n++){ // Loop over particle index.
        for (int i=0; i<3; i++){ // Loop over dimension.
            x[n][i]=0; // Set 0.
        }
    }
    return;
}

// Artificial seeder to have two particles input to collide.
void partpair_tester_IVP(int ns, double x[][3]){
    
    double pi = 3.14159265;

    // Particle 1.
    x[0][0]=0;
    x[0][1]=0;
    x[0][2]=0.5*pi;

    // Particle 2.
    x[1][0]=0;
    x[1][1]=5;
    x[1][2]=1.5*pi;

    return;
}

void no_OL_IVP(int ns, double x[][3], double xL, double xR, double yB, double yT, double rS, mt19937& gen){
    int check;
    double l;
    for (int n=0; n<ns; n++){ // Loop over particle index.
        check = 0;
        while(check != 1){
            check = 1;
            for (int i=0; i<3; i++){ // Loop over dimension.
                x[n][0] = urand(xL,xR,gen); // Set x-coordinate.
                x[n][1] = urand(yB,yT,gen); // Set x-coordinate.
            }
            for (int n2=0; n2<n;n2++){
                l = sqrt(pow(x[n][0] - x[n2][0],2) + pow(x[n][1] - x[n2][1],2));
                if( l < rS*2 ){
                    check = 0;
                }
            }
        }
    }
    return;
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

// WCA potential.
void WCA_force(int ns, double F[][3],double x[][3]){

    double sigma = 4; double eps = 0.1;
    double r;
    double dx; double dy;

    for(int n1=0; n1<ns; n1++){
        F[n1][0] = 0;
        F[n1][1] = 0;
        for(int n2=0; n2<ns; n2++){
            if(n1 != n2){
                dx = x[n2][0]-x[n1][0];
                dy = x[n2][1]-x[n1][1];
                r = sqrt(pow(dx,2) + pow(dy,2));
                if( r < pow(2,1/6)*sigma ){
                    F[n1][0] -= (24*eps*dx/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                    F[n1][1] -= (24*eps*dy/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                }
            }
        }
        cout << F[n1][0] << endl;
    }

    return;
}

// Basic implementation. Forward Euler.
void update_pos_basic(int ns,double x[][3], double vS, double sigR, double xL, double xR, double yB, double yT, double dT, mt19937& gen){
    for (int n=0; n<ns; n++){ // Loop over particle index.
        x[n][0]+= vS*dT*cos(x[n][2]); 
        x[n][1]+= vS*dT*sin(x[n][2]); 
        x[n][2]+= nrand(0,sigR,gen); 
    }
    return;
}

// Stochastic Runge-Kutta scheme, 2nd order.
void update_pos_RKII(int ns,double x[][3], double vS, double sigR, double xL, double xR, double yB, double yT, double dT, double DT, mt19937& gen){
    
    double    F1[ns][3];
    double    F2[ns][3];
    double     F[ns][3];
    double   f0 = -0.01;
    double xtemp[ns][3];

    WCA_force(ns,F,x);

    // Get F1.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        F1[n][0] = vS*dT*cos(x[n][2]) + f0*F[n][0]*pow(dT,2)*DT;
        F1[n][1] = vS*dT*sin(x[n][2]) + f0*F[n][1]*pow(dT,2)*DT;
        F1[n][2] = nrand(0,sigR,gen);

        xtemp[n][0] = x[n][0] + F1[n][0];
        xtemp[n][1] = x[n][1] + F1[n][1];
        xtemp[n][2] = x[n][2] + F1[n][2];
    }

    WCA_force(ns,F,xtemp);

    // Get F2.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        F2[n][0] = vS*dT*cos(x[n][2]+F1[n][2]) + f0*F[n][0]*DT;;
        F2[n][1] = vS*dT*sin(x[n][2]+F1[n][2]) + f0*F[n][1]*DT;;
        F2[n][2] = nrand(0,sigR,gen);
    }

    // Get x(t+\delta{t}).
    for (int n=0; n<ns; n++){ // Loop over particle index.
        x[n][0] += 0.5*(F1[n][0] + F2[n][0]);
        x[n][1] += 0.5*(F1[n][1] + F2[n][1]) ;
        x[n][2] += nrand(0,sigR,gen);
    }
}

// Even higher order solver?

#endif // 