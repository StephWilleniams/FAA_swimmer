#ifndef functions_IVP
#define functions_IVP

#include "header.hpp"
#include "functions_rng.hpp"

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
    x[0][0]=2.01;
    x[0][1]=0;
    x[0][2]=pi;

    // Particle 2.
    x[1][0]=-1.01;
    x[1][1]=0.35;
    x[1][2]=0;

    return;
}

void no_OL_IVP(int ns, double xs[][3], int np, double xp[][2], double xL, double xR, double yB, double yT, double rS, double rP, mt19937& gen){
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

// void no_OL_IVP(int ns, double x[][3], double xL, double xR, double yB, double yT, double rS, mt19937& gen){
//     int check;
//     double l;
//     for (int n=0; n<ns; n++){ // Loop over particle index.
//         check = 0;
//         while(check != 1){
//             check = 1;
//             for (int i=0; i<3; i++){ // Loop over dimension.
//                 x[n][0] = urand(xL,xR,gen); // Set x-coordinate.
//                 x[n][1] = urand(yB,yT,gen); // Set x-coordinate.
//             }
//             for (int n2=0; n2<n;n2++){
//                 l = sqrt(pow(x[n][0] - x[n2][0],2) + pow(x[n][1] - x[n2][1],2)); // Get distance to existing particle centers.
//                 if( l < rS*2 ){ // If center is to close re-roll position. 
//                     check = 0;
//                 }
//             }
//         }
//     }
//     return;
// }

#endif // 