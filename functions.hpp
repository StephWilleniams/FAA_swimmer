#ifndef functions
#define functions

#include "header.hpp"
#include "functions_rng.hpp"

// WCA potential.
void WCA_force(int ns, double F[][3], double x[][3], double xR, double xL, double yB, double yT){

    double sigma = 4; double eps = 1;
    double r; double dx; double dy;

    for(int n1=0; n1<ns; n1++){
        F[n1][0] = 0;
        F[n1][1] = 0;
        for(int n2=0; n2<ns; n2++){
            if(n1 != n2){
                dx = x[n2][0]-x[n1][0];
                dy = x[n2][1]-x[n1][1];
                if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;}
                if(abs(dy) > (yT-yB)/2){dy = (dy/abs(dy))*(abs(dx) - (yT-yB)) ;}
                r = sqrt(pow(dx,2) + pow(dy,2));
                if( r < sigma ){
                    F[n1][0] -= (24*eps*dx/pow(r,2))*(2*pow(sigma/r,6)-8*pow(sigma/r,12));
                    F[n1][1] -= (24*eps*dy/pow(r,2))*(2*pow(sigma/r,6)-8*pow(sigma/r,12));
                }
            }
        }
    }

    return;
}

// Basic implementation. Forward Euler. Not maintained.
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
    double     F[ns][3]; // Forces store.
    double   f0 = 0.001; // Force coefficient.
    double xtemp[ns][3]; // Temporary position store.

    WCA_force(ns,F,x, xR,xL,yB,yT); // Current particle potentials.

    // Get F1.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        F1[n][0] = vS*dT*cos(x[n][2]) - f0*F[n][0]*pow(dT,2)*DT;
        F1[n][1] = vS*dT*sin(x[n][2]) - f0*F[n][1]*pow(dT,2)*DT;
        F1[n][2] = nrand(0,sigR,gen);

        // Get intermediate position.
        xtemp[n][0] = x[n][0] + F1[n][0];
        xtemp[n][1] = x[n][1] + F1[n][1];
        xtemp[n][2] = x[n][2] + F1[n][2];
        // Boundary conditions (inefficient).
        while(xtemp[n][0] > xR){xtemp[n][0] -= (xR-xL);}
        while(xtemp[n][0] < xL){xtemp[n][0] += (xR-xL);}
        while(xtemp[n][1] > yT){xtemp[n][1] -= (yT-yB);}
        while(xtemp[n][1] < yB){xtemp[n][1] += (yT-yB);}
    }

    WCA_force(ns,F,xtemp, xR,xL,yB,yT); // Intermediate particle potentials.

    // Get F2.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        F2[n][0] = vS*dT*cos(x[n][2]+F1[n][2]) + f0*F[n][0]*DT;;
        F2[n][1] = vS*dT*sin(x[n][2]+F1[n][2]) + f0*F[n][1]*DT;;
        F2[n][2] = nrand(0,sigR,gen);
    }

    // Get x(t+\delta{t}).
    for (int n=0; n<ns; n++){ // Loop over particle index.
        // Get RKII-displacement.
        x[n][0] += 0.5*(F1[n][0] + F2[n][0]);
        x[n][1] += 0.5*(F1[n][1] + F2[n][1]) ;
        x[n][2] += nrand(0,sigR,gen);
        // Boundary conditions (inefficient).
        while(x[n][0] > xR){x[n][0] -= (xR-xL);}
        while(x[n][0] < xL){x[n][0] += (xR-xL);}
        while(x[n][1] > yT){x[n][1] -= (yT-yB);}
        while(x[n][1] < yB){x[n][1] += (yT-yB);}
    }
}

// Even higher order solver?

#endif // 