#ifndef functions_IVP
#define functions_IVP

#include "header.hpp"
#include "functions_rng.hpp"

void set_active_geometry(double rA, double pol, double rSeg[], double sDV[]){
    
    rSeg[0] = (4*sqrt(2)/(5*sqrt(2)+2))*rA;
    rSeg[1] = (4*sqrt(2)/(5*sqrt(2)+2))*rA/sqrt(2);
    rSeg[2] = (4*sqrt(2)/(5*sqrt(2)+2))*rA/2;
    sDV[0]  = double(pol)*(rSeg[0] - rA);
    sDV[1]  = double(pol)*(rSeg[0]*2 - rA);
    sDV[2]  = double(pol)*(rSeg[0]*2 + rSeg[1] - rA);
}

void zero_IVP(int ns, double x[][4]){
    for (int n=0; n<ns; n++){ // Loop over particle index.
        for (int i=0; i<3; i++){ // Loop over dimension.
            x[n][i]=0; // Set 0.
        }
    }
    return;
}

// Artificial seeder to have two particles input to collide.
void partpair_tester_IVP_2(int ns, double xs[][4],int np, double xp[][3]){
    
    //double pi = 3.14159265;

    // Set active particle 1.
    xs[0][0]=10;
    xs[0][1]=-6.5;
    xs[0][2]=(3*M_PI/2)+0.01;

    // // Set active particle 2.
    // xs[1][0]=0;
    // xs[1][1]=6.5;
    // xs[1][2]=M_PI/2+0.1;

    // // // Set active particle 3.
    // xs[2][0]=20;
    // xs[2][1]=6.5;
    // xs[2][2]=M_PI/2+0.1;
     
    // // // Set active particle 4.
    // xs[3][0]=30;
    // xs[3][1]=6.5;
    // xs[3][2]=M_PI/2+0.1;

    // // // Set active particle 5.
    // xs[4][0]=40;
    // xs[4][1]=6.5;
    // xs[4][2]=M_PI/2+0.1;

    // Set passive particle.
    xp[0][0] = 100;
    xp[0][1] = 100;

    return;
}

void no_OL_IVP(int ns, double xs[][4], int np, double xp[][3], double xL, double xR, double yB, double yT, double rS, double rP, mt19937& gen){
    int check;
    double l;

    // Seed the active particles (full no overlap)
    for (int n=0; n<ns; n++){ // Loop over particle index.
        check = 0;
        while(check != 1){
            check = 1;
            for (int i=0; i<3; i++){ // Loop over dimension.
                xs[n][0] = urand(xL,xR,gen); // Set x-coordinate.
                xs[n][1] = urand(yB,yT,gen); // Set x-coordinate.
                xs[n][2] = urand(0,2*M_PI,gen); // Set x-coordinate.
            }
            for (int n2=0; n2<n;n2++){
                l = sqrt(pow(xs[n][0] - xs[n2][0],2) + pow(xs[n][1] - xs[n2][1],2)); // Get distance to existing particle centers.
                if( l < rS*2 ){ // If center is too close re-roll position. 
                    check = 0;
                }
            }
        }
    }

    // Seed "ghost" passive particles, no overlap with active particles.
    for (int n=0; n<np; n++){ // Loop over particle index.
        check = 0;
        while(check != 1){
            check = 1;
            for (int i=0; i<3; i++){ // Loop over dimension.
                xp[n][0] = urand(xL,xR,gen); // Set x-coordinate.
                xp[n][1] = urand(yB,yT,gen); // Set x-coordinate.
            }
            for (int n2=0; n2<ns;n2++){
                l = sqrt(pow(xp[n][0] - xs[n2][0],2) + pow(xp[n][1] - xs[n2][1],2)); // Get distance to existing particle centers.
                if( l < rS+rP ){ // If center is too close re-roll position. 
                    check = 0;
                }
            }
        }
    }

    return;
}

#endif // 