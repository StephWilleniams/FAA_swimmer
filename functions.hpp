#ifndef functions
#define functions

#include "header.hpp"
#include "functions_rng.hpp"

// WCA potential, 3-part active particle, includes passive.
void WCA_force_2(int ns, double Fa[][3], double xs[][3], double rA, double rSeg[], double sDV[], double polarisation, int np, double Fp[][2], double xp[][2], double rP, double xR, double xL, double yB, double yT){

    double sigma; double eps = 0.05;
    double r; double x1; double x2; double y1; double y2; double y3; double dx; double dy;

    // Calulate the forces on the active particles.
    //#pragma omp // parallel for
    for(int n1=0; n1<ns; n1++){

        // Reset force on current particle.
        Fa[n1][0] = 0;
        Fa[n1][1] = 0;

        // Get comparisons of active-active forces.
        for(int n2=0; n2<ns; n2++){
            if(n1 != n2){

                // loop over the 3 segments in each swimmer
                for(int seg1=0; seg1<3; seg1++){
                    for(int seg2=0; seg2<3; seg2++){

                        sigma = pow(2,-1/6)*(rSeg[seg1]+rSeg[seg2]); // Get the interaction parameter

                        // Get the coordinates of the relative segments.
                        x1 = xs[n1][0] + polarisation*sDV[seg1]*cos(xs[n1][2]);
                        x2 = xs[n2][0] + polarisation*sDV[seg2]*cos(xs[n2][2]);
                        y1 = xs[n1][1] + polarisation*sDV[seg1]*sin(xs[n1][2]);
                        y2 = xs[n2][1] + polarisation*sDV[seg2]*sin(xs[n2][2]);
                    
                        dx = x2-x1; // x-distance between centers.
                        dy = y2-y1; // y-distance between centers.

                        if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
                        if(abs(dy) > (yT-yB)/2){dy = (dy/abs(dy))*(abs(dx) - (yT-yB)) ;} // Get wrapped y

                        r = sqrt(pow(dx,2) + pow(dy,2)); // Get the distance between the centers.

                        // Implement the WCA potential.
                        if( r < pow(2,1/6)*sigma ){
                            Fa[n1][0] -= (24*eps*dx/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                            Fa[n1][1] -= (24*eps*dy/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                            // Include the torque here.
                        }
                    }
                }
            }
        }

        // Get comparisons of the active-boundary forces.
        for(int seg1=0; seg1<3; seg1++){
            sigma = pow(2,-1/6)*(rSeg[seg1]);

            // Get the coordinates of the relative segments.
            x1 = xs[n1][0] + polarisation*sDV[seg1]*cos(xs[n1][2]);
            x2 = xs[n1][0] + polarisation*sDV[seg1]*cos(xs[n1][2]);
            y1 = xs[n1][1] + polarisation*sDV[seg1]*sin(xs[n1][2]);
            y2 = yT;
            y3 = yB;

            dx = x2-x1; // x-distance between centers.

            if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
            if(abs(y3-y1) > abs(y2-y1)){dy = y2-y1;}else{dy = y3-y1;} // Get wrapped y

            // actually do something????

        }

        // Get comparisons of active-passive forces (back-force).

    }

    // Calulate the forces on the passive particles.
    for(int n1=0; n1<np; n1++){

        // Reset force on current particle.
        Fp[n1][0] = 0;
        Fp[n1][1] = 0;

        // Get comparisons of passive-active forces.
        for(int n2=0; n2<ns; n2++){
            for(int seg1=0; seg1<3; seg1++){

                sigma = pow(2,-1/6)*(rP+rSeg[seg1]); // Get the interaction parameter

                // Get the coordinates of the comparison particles.
                x1 = xp[n1][0];
                x2 = xs[n2][0] + polarisation*sDV[seg1]*cos(xs[n2][2]);
                y1 = xp[n1][1];
                y2 = xs[n2][1] + polarisation*sDV[seg1]*sin(xs[n2][2]);

                dx = x2-x1; // x-distance between centers.
                dy = y2-y1; // y-distance between centers.

                sigma = pow(2,-1/6)*(rP+rA); // Get the interaction parameter.

                if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
                if(abs(dy) > (yT-yB)/2){dy = (dy/abs(dy))*(abs(dx) - (yT-yB)) ;} // Get wrapped y

                r = sqrt(pow(dx,2) + pow(dy,2)); // Distance between the centers.

                // Implement the WCA potential.
                if( r < sigma ){
                    Fp[n1][0] -= (24*eps*dx/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12));
                    Fp[n1][1] -= (24*eps*dy/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12));
                }

            }
        }

        // Get comparisons of the passive-boundary forces.
        // Get the coordinates of the relative segments.
        x1 = xp[n1][0];
        x2 = xp[n1][0];
        y1 = xp[n1][1];
        y2 = yT;
        y3 = yB;

        dx = x2-x1; // x-distance between centers.

        if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
        if(abs(y3-y1) > abs(y2-y1)){dy = y2-y1;}else{dy = y3-y1;} // Get wrapped y

        // actually do something????

        // Get comparisons of passive-passive forces, currently omitted.
        // This can be added in future.

    }

    return;
}

/* /////////////// */

// Stochastic Runge-Kutta scheme, 2nd order.
void update_pos_RKII(int ns, double xs[][3], double rA, double rSeg[], double sDV[], double polarisation, int np, double xp[][2], double rP, double vS, double sigR, double sigT, double xL, double xR, double yB, double yT, double dT, double DT, mt19937& gen){
    
    // Active forces
    double    F1a[ns][3];
    double    F2a[ns][3];
    double     Fa[ns][3]; // Forces store.

    // Passive forces
    double    F1p[np][4];
    double    F2p[np][2];
    double     Fp[np][2]; // Forces store.

    // Force constants
    double   f0 = 0.001; // Force coefficient.

    // Position stores
    double xtemp_a[ns][3]; // Temporary position store.
    double xtemp_p[np][2]; // Temporary position store.

    WCA_force_2(ns,Fa,xs,rA,rSeg,sDV,polarisation,np,Fp,xp,rP,xR,xL,yB,yT); // Current particle potentials.

    // Get F1, active.
    for (int n=0; n<ns; n++){ // Loop over particle index.

        // Figure out how to implement fx and fy here.
        F1a[n][0] = vS*dT*cos(xs[n][2]) - f0*Fa[n][0]*pow(dT,2);
        F1a[n][1] = vS*dT*sin(xs[n][2]) - f0*Fa[n][1]*pow(dT,2);
        // Scattering event here.
        F1a[n][2] = nrand(0,sigR,gen);

        // Get intermediate position.
        xtemp_a[n][0] = xs[n][0] + F1a[n][0];
        xtemp_a[n][1] = xs[n][1] + F1a[n][1];
        xtemp_a[n][2] = xs[n][2] + F1a[n][2];

        // Boundary conditions.
        xtemp_a[n][0] = xtemp_a[n][0] - (xR-xL)*((xtemp_a[n][0])/abs(xtemp_a[n][0]))*floor(abs(xtemp_a[n][0])/xR);
        xtemp_a[n][1] = xtemp_a[n][1] - (yT-yB)*((xtemp_a[n][1])/abs(xtemp_a[n][1]))*floor(abs(xtemp_a[n][1])/yT);

    }

    // Get F1, passive.
    for (int n=0; n<np; n++){ // Loop over particle index.

        F1p[n][0] =  -f0*Fp[n][0]*pow(dT,2);
        F1p[n][1] =  -f0*Fp[n][1]*pow(dT,2);
        F1p[n][2] = nrand(0,sigT,gen);
        F1p[n][3] = nrand(0,sigT,gen);

        // Get intermediate position.
        xtemp_p[n][0] = xp[n][0] + F1p[n][0] + F1p[n][2];
        xtemp_p[n][1] = xp[n][1] + F1p[n][1] + F1p[n][3];

        // Boundary conditions.
        xtemp_p[n][0] = xtemp_p[n][0] - (xR-xL)*((xtemp_p[n][0])/abs(xtemp_p[n][0]))*floor(abs(xtemp_p[n][0])/xR);
        xtemp_p[n][1] = xtemp_p[n][1] - (yT-yB)*((xtemp_p[n][1])/abs(xtemp_p[n][1]))*floor(abs(xtemp_p[n][1])/yT);

    }

    WCA_force_2(ns,Fa,xtemp_a,rA,rSeg,sDV,polarisation,np,Fp,xtemp_p,rP,xR,xL,yB,yT); // Intermediate particle potentials.

    // Get F2, active.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        // Figure out how to implement fx and fy here.
        F2a[n][0] = vS*dT*cos(xs[n][2]+F1a[n][2]) - f0*Fa[n][0]*pow(dT,2);
        F2a[n][1] = vS*dT*sin(xs[n][2]+F1a[n][2]) - f0*Fa[n][1]*pow(dT,2);
        // Scattering event here.
        F2a[n][2] = nrand(0,sigR,gen);
    }

    // Get F2, passive.
    for (int n=0; n<np; n++){ // Loop over particle index.
        F2p[n][0] = - f0*Fp[n][0]*pow(dT,2);
        F2p[n][1] = - f0*Fp[n][1]*pow(dT,2);
    }

    // Get x(t+\delta{t}), active.
    for (int n=0; n<ns; n++){ // Loop over particle index.
        // Get RKII-displacement.
        xs[n][0] += 0.5*(F1a[n][0] + F2a[n][0]);
        xs[n][1] += 0.5*(F1a[n][1] + F2a[n][1]);
        // Scattering event here.
        xs[n][2] += nrand(0,sigR,gen);
        
        // Boundary conditions.
        xs[n][0] = xs[n][0] - (xR-xL)*(xs[n][0]/abs(xs[n][0]))*floor(abs(xs[n][0])/xR);
        xs[n][1] = xs[n][1] - (yT-yB)*(xs[n][1]/abs(xs[n][1]))*floor(abs(xs[n][1])/yT);

    }

    // Get x(t+\delta{t}), passive.
    for (int n=0; n<np; n++){ // Loop over particle index.
        // Get RKII-displacement.
        xp[n][0] += 0.5*(F1p[n][0] + F2p[n][0]) + F1p[n][2];
        xp[n][1] += 0.5*(F1p[n][1] + F2p[n][1]) + F1p[n][3];
        
        // Boundary conditions.
        xp[n][0] = xp[n][0] - xR*((xp[n][0])/abs(xp[n][0]))*floor(abs(xp[n][0])/xR);
        xp[n][1] = xp[n][1] - yT*((xp[n][1])/abs(xp[n][1]))*floor(abs(xp[n][1])/yT);

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