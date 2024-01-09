#ifndef functions
#define functions

#include "header.hpp"
#include "functions_rng.hpp"

// WCA potential, 3-part active particle, includes passive.
void WCA_force_2(int ns, double Fa[][3], double xs[][3], double rA, double rSeg[], double sDV[], double polarisation, int np, double Fp[][2], double xp[][2], double xR, double xL, double yB, double yT){

    double sigma; double eps = 0.5;
    double r; double x1; double x2; double y1; double y2; double dx; double dy;

    // Calulate the forces on the active particles.
    //#pragma omp parallel for
    for(int n1=0; n1<ns; n1++){

        Fa[n1][0] = 0;
        Fa[n1][1] = 0;

        // Get comparisons of active-active forces.
        for(int n2=0; n2<ns; n2++){
            if(n1 != n2){

                // loop over the 3 segments in each swimmer
                for(int seg1=0; seg1<3; seg1++){
                    for(int seg2=0; seg2<3; seg2++){

                        sigma = pow(2,-1/6)*(rSeg[seg1]+rSeg[seg2]); // Modift this 2

                        x1 = xs[n1][0] + polarisation*sDV[seg1]*cos(xs[n1][2]);
                        x2 = xs[n2][0] + polarisation*sDV[seg2]*cos(xs[n2][2]);
                        y1 = xs[n1][1] + polarisation*sDV[seg1]*sin(xs[n1][2]);
                        y2 = xs[n2][1] + polarisation*sDV[seg2]*sin(xs[n2][2]);
                    
                        dx = x2-x1; // x-distance between centers.
                        dy = y2-y1; // y-distance between centers.

                        if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
                        if(abs(dy) > (yT-yB)/2){dy = (dy/abs(dy))*(abs(dx) - (yT-yB)) ;} // Get wrapped y

                        r = sqrt(pow(dx,2) + pow(dy,2));

                        if( r < pow(2,1/6)*sigma ){
                            Fa[n1][0] -= (24*eps*dx/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                            Fa[n1][1] -= (24*eps*dy/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                        }
                    }
                }
            }
        }

        // Get comparisons of active-passive forces (back-force).

    }

    // Calulate the forces on the passive particles.
    for(int n1=0; n1<np; n1++){

        Fp[n1][0] = 0;
        Fp[n1][1] = 0;

        // Get comparisons of passive-active forces.
        for(int n2=0; n2<ns; n2++){
            if(n1 != n2){

                // loop over the 3 segments in each swimmer

                dx = xs[n2][0]-xp[n1][0]; // x-distance between centers.
                dy = xs[n2][1]-xp[n1][1]; // y-distance between centers.

                if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
                if(abs(dy) > (yT-yB)/2){dy = (dy/abs(dy))*(abs(dx) - (yT-yB)) ;} // Get wrapped y

                r = sqrt(pow(dx,2) + pow(dy,2));

                if( r < sigma ){
                    Fa[n1][0] -= (24*eps*dx/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12));
                    Fa[n1][1] -= (24*eps*dy/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12));
                }

            }
        }

        // Get comparisons of passive-passive forces.

    }

    return;
}

// Stochastic Runge-Kutta scheme, 2nd order.
void update_pos_RKII(int ns, double xs[][3], double rA, double rSeg[], double sDV[], double polarisation, int np, double xp[][2], double vS, double sigR, double sigT, double xL, double xR, double yB, double yT, double dT, double DT, mt19937& gen){
    
    // Active forces
    double    F1a[ns][3];
    double    F2a[ns][3];
    double     Fa[ns][3]; // Forces store.

    // Passive forces
    double    F1p[ns][2];
    double    F2p[ns][2];
    double     Fp[ns][2]; // Forces store.

    // Force constants
    double   f0 = 0.001; // Force coefficient.

    // Position stores
    double xtemp_a[ns][3]; // Temporary position store.
    double xtemp_p[np][2]; // Temporary position store.

    WCA_force_2(ns,Fa,xs,rA,rSeg,sDV,polarisation,np,Fp,xp,xR,xL,yB,yT); // Current particle potentials.

    // Get F1.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        F1a[n][0] = vS*dT*cos(xs[n][2]) - f0*Fa[n][0]*pow(dT,2);
        F1a[n][1] = vS*dT*sin(xs[n][2]) - f0*Fa[n][1]*pow(dT,2);
        F1a[n][2] = nrand(0,sigR,gen);

        // Get intermediate position.
        xtemp_a[n][0] = xs[n][0] + F1a[n][0];
        xtemp_a[n][1] = xs[n][1] + F1a[n][1];
        xtemp_a[n][2] = xs[n][2] + F1a[n][2];

        // Boundary conditions (inefficient).
        while(xtemp_a[n][0] > xR){xtemp_a[n][0] -= (xR-xL);}
        while(xtemp_a[n][0] < xL){xtemp_a[n][0] += (xR-xL);}
        while(xtemp_a[n][1] > yT){xtemp_a[n][1] -= (yT-yB);}
        while(xtemp_a[n][1] < yB){xtemp_a[n][1] += (yT-yB);}
    }

    WCA_force_2(ns,Fa,xtemp_a,rA,rSeg,sDV,polarisation,np,Fp,xtemp_p,xR,xL,yB,yT); // Intermediate particle potentials.

    // Get F2.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        F2a[n][0] = vS*dT*cos(xs[n][2]+F1a[n][2]) - f0*Fa[n][0]*pow(dT,2);
        F2a[n][1] = vS*dT*sin(xs[n][2]+F1a[n][2]) - f0*Fa[n][1]*pow(dT,2);
        F2a[n][2] = nrand(0,sigR,gen);
    }

    // Get x(t+\delta{t}).
    for (int n=0; n<ns; n++){ // Loop over particle index.
        // Get RKII-displacement.
        xs[n][0] += 0.5*(F1a[n][0] + F2a[n][0]);
        xs[n][1] += 0.5*(F1a[n][1] + F2a[n][1]) ;
        xs[n][2] += nrand(0,sigR,gen);
        
        // Boundary conditions (inefficient).
        while(xs[n][0] > xR){xs[n][0] -= (xR-xL);}
        while(xs[n][0] < xL){xs[n][0] += (xR-xL);}
        while(xs[n][1] > yT){xs[n][1] -= (yT-yB);}
        while(xs[n][1] < yB){xs[n][1] += (yT-yB);}
    }
}

// Even higher order solver?

/* /////////////// */
/* Older functions */
/* /////////////// */

// Basic implementation. Forward Euler. Not maintained.
void update_pos_basic(int ns,double x[][3], double vS, double sigR, double xL, double xR, double yB, double yT, double dT, mt19937& gen){
    for (int n=0; n<ns; n++){ // Loop over particle index.
        x[n][0]+= vS*dT*cos(x[n][2]); 
        x[n][1]+= vS*dT*sin(x[n][2]); 
        x[n][2]+= nrand(0,sigR,gen); 
    }
    return;
}

// WCA potential, 1-part active particle, no passive, old.
void WCA_force(int ns, double F[][3], double x[][3], double xR, double xL, double yB, double yT){

    double sigma = 1; double eps = 1;
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
                    F[n1][0] -= (24*eps*dx/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12));
                    F[n1][1] -= (24*eps*dy/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12));
                }
            }
        }
    }

    return;
}

#endif // 