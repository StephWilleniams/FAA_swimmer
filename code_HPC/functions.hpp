#ifndef functions
#define functions

#include "header.hpp"
#include "functions_rng.hpp"

// WCA potential, 3-part active particle, includes passive.
void WCA_force_2(int na, double Fa[][3], double xa[][4], double rA, double rSeg[], double sDV[], int np, double Fp[][2], double xp[][2], double rP, double xR, double xL, double yB, double yT){

    double sigma; double eps = 0.005;
    double r; double x1; double x2; double y1; double y2; double y3; double dx; double dy;

    // Calulate the forces on the active particles.
    //#pragma omp // parallel for
    for(int n1=0; n1<na; n1++){

        // Reset force on current particle.
        Fa[n1][0] = 0;
        Fa[n1][1] = 0;
        Fa[n1][2] = 0;

        // Get comparisons of active-active forces.
        for(int n2=0; n2<na; n2++){
            if(n1 != n2){

                // loop over the 3 segments in each swimmer
                for(int seg1=0; seg1<3; seg1++){
                    for(int seg2=0; seg2<3; seg2++){

                        sigma = pow(2,-1/6)*(rSeg[seg1]+rSeg[seg2]); // Get the interaction parameter

                        // Get the coordinates of the relative segments.
                        x1 = xa[n1][0] + sDV[seg1]*cos(xa[n1][2]);
                        x2 = xa[n2][0] + sDV[seg2]*cos(xa[n2][2]);
                        y1 = xa[n1][1] + sDV[seg1]*sin(xa[n1][2]);
                        y2 = xa[n2][1] + sDV[seg2]*sin(xa[n2][2]);
                    
                        dx = x2-x1; // x-distance between centers.
                        dy = y2-y1; // y-distance between centers.

                        if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
                        if(abs(dy) > (yT-yB)/2){dy = (dy/abs(dy))*(abs(dx) - (yT-yB)) ;} // Get wrapped y

                        r = sqrt(pow(dx,2) + pow(dy,2)); // Get the distance between the centers.

                        // Implement the WCA potential.
                        if( r < pow(2,1/6)*sigma ){
                            Fa[n1][0] -= (24*eps*dx/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                            Fa[n1][1] -= (24*eps*dy/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12));
                            Fa[n1][2] -= 4*eps*sDV[seg1]*(12*pow(sigma/r,11)-6*pow(sigma/r,6))*(dx*cos(xa[n1][2])-dy*sin(xa[n1][2]))/pow(r,3);
                        }
                    }
                }
            }
        }

        // Get comparisons of the active-boundary forces.
        for(int seg1=0; seg1<3; seg1++){

            sigma = pow(2,-1/6)*(rSeg[seg1]); // Get the interaction parameter.

            // Get the coordinates of the relative segments.
            x1 = xa[n1][0] + sDV[seg1]*cos(xa[n1][2]);
            x2 = xa[n1][0] + sDV[seg1]*cos(xa[n1][2]);
            y1 = xa[n1][1] + sDV[seg1]*sin(xa[n1][2]);
            y2 = yT;
            y3 = yB;

            dx = x2-x1; // x-distance between centers.

            if(abs(dx) > (xR-xL)/2){dx = (dx/abs(dx))*(abs(dx) - (xR-xL)) ;} // Get wrapped x
            if(abs(y3-y1) > abs(y2-y1)){dy = y2-y1;}else{dy = y3-y1;} // Get wrapped y
            r = abs(dy); // Get the displacement between the segment and boundary.
            Fa[n1][1] -= (24*eps*dy/pow(r,2))*(pow(sigma/r,6)-2*pow(sigma/r,12)); // Implement the force.
            Fa[n1][2] -= 4*eps*sDV[seg1]*(12*pow(sigma/r,11)-6*pow(sigma/r,6))*(-dy*sin(xa[n1][2]))/pow(r,3);
        }

        // Get comparisons of active-passive forces (back-force).

    }

    // Calulate the forces on the passive particles.
    for(int n1=0; n1<np; n1++){

        // Reset force on current particle.
        Fp[n1][0] = 0;
        Fp[n1][1] = 0;

        // Get comparisons of passive-active forces.
        for(int n2=0; n2<na; n2++){
            for(int seg1=0; seg1<3; seg1++){

                sigma = pow(2,-1/6)*(rP+rSeg[seg1]); // Get the interaction parameter

                // Get the coordinates of the comparison particles.
                x1 = xp[n1][0];
                x2 = xa[n2][0] + sDV[seg1]*cos(xa[n2][2]);
                y1 = xp[n1][1];
                y2 = xa[n2][1] + sDV[seg1]*sin(xa[n2][2]);

                dx = x2-x1; // x-distance between centers.
                dy = y2-y1; // y-distance between centers.

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
        sigma = rP; // Get the interaction parameter.
        // Get the coordinates of the relative segments.
        y1 = xp[n1][1];
        y2 = yT;
        y3 = yB;
        if(abs(y3-y1) > abs(y2-y1)){dy = y2-y1;}else{dy = y3-y1;} // Get wrapped y.
        r = abs(dy); // Get the distance to the boundary.
        Fp[n1][1] -= (24*eps*dy/pow(r,2))*(1*pow(sigma/r,6)-2*pow(sigma/r,12)); // Implement the force.

        // Get comparisons of passive-passive forces, currently omitted.
        // This can be added in future.

    }

    return;
}

/* /////////////// */

// Stochastic Runge-Kutta scheme, 2nd order.
void update_pos_RKII(int na, double xa[][4], double rA, double rSeg[], double sDV[], int np, double xp[][2], double rP, double vS, double sigR, double sigT, double xL, double xR, double yB, double yT, double f0, double fR, double fricPar, double fricPerp, double kickFreq, double kickStr, double dT, double DT, mt19937& gen){
    
    // Active forces
    double    F1a[na][4];
    double    F2a[na][3];
    double     Fa[na][3]; // Forces store.

    // Passive forces
    double    F1p[np][4];
    double    F2p[np][2];
    double     Fp[np][2]; // Forces store.

    // Position stores
    double xtemp_a[na][4]; // Temporary position store.
    double xtemp_p[np][2]; // Temporary position store.

    // Friction array stores
    double ct; double st; double fx; double fy; double det; double storeA; double storeB; double storeC; double storeD; // Friction rotation arrays.

    WCA_force_2(na,Fa,xa,rA,rSeg,sDV,np,Fp,xp,rP,xR,xL,yB,yT); // Current particle potentials.

    // Get F1, active.
    for (int n=0; n<na; n++){ // Loop over particle index.

        // Get friction in frame of current active particle.
        ct = cos(xa[n][2]); st = sin(xa[n][2]);
        // Get the components of the friction tensor [a,b;c,d], respectively..
        storeA = fricPar*pow(ct,2) + fricPerp*(1-pow(ct,2));
        storeB = (fricPar - fricPerp)*ct*st;storeC = storeB;
        storeD = fricPar*pow(st,2) + fricPerp*(1-pow(st,2));
        det = 1/(f0*(storeA*storeD - storeB*storeC)); // Get determinant of friction tensor.
        // Apply the mobility (inverse friction) to the components of the calculated force.
        fx = det*(storeD*Fa[n][0] - storeB*Fa[n][1])*pow(dT,2);
        fy = det*(-storeC*Fa[n][0] + storeA*Fa[n][1])*pow(dT,2);

        // Implement noise, kicks, fx, fy, and ftheta here.
        F1a[n][0] = vS*dT*cos(xa[n][2]) - fx;
        F1a[n][1] = vS*dT*sin(xa[n][2]) - fy;
        F1a[n][3] = nrand(0,sigR,gen); // Generate and store noise.
        F1a[n][2] = F1a[n][3] + (1/(fR*f0))*Fa[n][2] - kickStr*cos(xa[n][2])/abs(cos(xa[n][2]))*xa[n][3];

        // Get intermediate position.
        xtemp_a[n][0] = xa[n][0] + F1a[n][0];
        xtemp_a[n][1] = xa[n][1] + F1a[n][1];
        xtemp_a[n][2] = xa[n][2] + F1a[n][2];

        // Boundary conditions.
        xtemp_a[n][0] = xtemp_a[n][0] - (xR-xL)*floor((xtemp_a[n][0]+xR)/(xR-xL));
        xtemp_a[n][1] = xtemp_a[n][1] - (yT-yB)*floor((xtemp_a[n][1]+yT)/(yT-yB));
        xtemp_a[n][2] = xtemp_a[n][2] - (20*M_PI)*floor((xtemp_a[n][2]+10*M_PI)/(20*M_PI));
    }
    // Get F1, passive.
    for (int n=0; n<np; n++){ // Loop over particle index.

        F1p[n][0] =  -(1/f0)*Fp[n][0]*pow(dT,2);
        F1p[n][1] =  -(1/f0)*Fp[n][1]*pow(dT,2);
        F1p[n][2] = nrand(0,sigT,gen);
        F1p[n][3] = nrand(0,sigT,gen);

        // Get intermediate position.
        xtemp_p[n][0] = xp[n][0] + F1p[n][0] + F1p[n][2];
        xtemp_p[n][1] = xp[n][1] + F1p[n][1] + F1p[n][3];

        // Boundary conditions.
        xtemp_p[n][0] = xtemp_p[n][0] - (xR-xL)*floor((xtemp_p[n][0]+xR)/(xR-xL));
        xtemp_p[n][1] = xtemp_p[n][1] - (yT-yB)*floor((xtemp_p[n][1]+yT)/(yT-yB));
    }

    WCA_force_2(na,Fa,xtemp_a,rA,rSeg,sDV,np,Fp,xtemp_p,rP,xR,xL,yB,yT); // Intermediate particle potentials.

    // Get F2, active.
    for (int n=0; n<na; n++){ // Loop over particle index.

        // Get friction in frame of current active particle.
        ct = cos(xtemp_a[n][2]); st = sin(xtemp_a[n][2]);
        // Get the components of the friction tensor [a,b;c,d], respectively..
        storeA = fricPar*pow(ct,2) + fricPerp*(1-pow(ct,2));
        storeB = (fricPar - fricPerp)*ct*st;storeC = storeB;
        storeD = fricPar*pow(st,2) + fricPerp*(1-pow(st,2));
        det = 1/(f0*(storeA*storeD - storeB*storeC)); // Get determinant of friction tensor.
        // Apply the inverse to the components of the calculated force.
        fx = det*(storeD*Fa[n][0] - storeB*Fa[n][1])*pow(dT,2);
        fy = det*(-storeC*Fa[n][0] + storeA*Fa[n][1])*pow(dT,2);

        // Implement kicks, fx, fy, and ftheta here (noise is not needed).
        F2a[n][0] = vS*dT*cos(xa[n][2]+F1a[n][2]) - fx;
        F2a[n][1] = vS*dT*sin(xa[n][2]+F1a[n][2]) - fy;
        F2a[n][2] = (1/(fR*f0))*Fa[n][2] - kickStr*cos(xa[n][2])/abs(cos(xa[n][2]))*xa[n][3];
    }

    // Get F2, passive.
    for (int n=0; n<np; n++){ // Loop over particle index.
        F2p[n][0] = -(1/f0)*Fp[n][0]*pow(dT,2);
        F2p[n][1] = -(1/f0)*Fp[n][1]*pow(dT,2);
    }

    // Get x(t+\delta{t}), active.
    for (int n=0; n<na; n++){ // Loop over particle index.

        // Get RKII-displacement.
        xa[n][0] += 0.5*(F1a[n][0] + F2a[n][0]);
        xa[n][1] += 0.5*(F1a[n][1] + F2a[n][1]);
        if( (xa[n][1] > yT - 1.5*rA) || (xa[n][1] < yB + 1.5*rA) ){
            xa[n][3] = ((xa[n][1])/abs(xa[n][1]))*max(xa[n][3],ceil(urand(0,1,gen)-(exp(-kickFreq)))); // Determine if a kick is happening or if one is initiated.
        } else {
            xa[n][3] = 0;
        }
        xa[n][2] += 0.5*(F1a[n][2] + F2a[n][2]) + F1a[n][3];
        
        // Boundary conditions.
        xa[n][0] = xa[n][0] - (xR-xL)*floor((xa[n][0]+xR)/(xR-xL));
        xa[n][1] = xa[n][1] - (yT-yB)*floor((xa[n][1]+yT)/(yT-yB));
        xa[n][2] = xa[n][2] - (20*M_PI)*floor((xa[n][2]+10*M_PI)/(20*M_PI));
    }

    // Get x(t+\delta{t}), passive.
    for (int n=0; n<np; n++){ // Loop over particle index.
        // Get RKII-displacement.
        xp[n][0] += 0.5*(F1p[n][0] + F2p[n][0]) + F1p[n][2];
        xp[n][1] += 0.5*(F1p[n][1] + F2p[n][1]) + F1p[n][3];
        
        // Boundary conditions.
        xp[n][0] = xp[n][0] - (xR-xL)*floor((xp[n][0]+xR)/(xR-xL));
        xp[n][1] = xp[n][1] - (yT-yB)*floor((xp[n][1]+yT)/(yT-yB));
    }
}

#endif // 