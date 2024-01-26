#ifndef constants
#define constants

#include "header.hpp"

// Note: Captialized variables are dimensional (lower case ND).

// Time constants (s).
double T = 100; // Total runtime (s).
double DT = 0.0005; // Step size (s) (ND time). 
// Time constants ND.
double dT = DT/DT;
int nSteps = int(T/DT); // Number of calculated steps.

// Length constants (um).
double XR = 10;  // System right side.
double XL = -10; // System left side.
double YT = 10;  // System top side.
double YB = -10; // System bottom side.
double RA = 1; // Active particle size (this is the ND length)
double RP = 0.5; // Passive particle size.
// Length constants ND.
double xR = XR/RA;
double xL = XL/RA;
double yT = YT/RA;
double yB = YB/RA;
double rA = RA/RA;
double rP = RP/RA;

// Particle segement sizes (ND).
int pol = -1; // Active particle polarisation. +1 for pusher, -1 for puller.
double rSeg[3];
double  sDV[3];

// Velocity constants (um/s).
double VA = 100; // Swimmer speed.
// Velocity constants ND.
double vA = VA*DT/RA;

// Diffusion constants (xx^2/s).
double DR = 0.2; // Swimmer rotational diffusivity, rad^2/s.
double DTherm = 0.1; // Passive particle diffusivity, um^2/s.
// Diffusion constants ND.
double dR = DR*DT; 
double dTherm = DTherm*DT/pow(RA,2); 
double sigR = sqrt(dR);
double sigT = sqrt(4*dTherm*dT);

// Friction constants.
double a = 1.6; // Aspect ratio.
double f0 = 3; // Base friction.
double fR = M_PI*pow(a,2)/(3*(log(a) - 0.662 + (0.917/a)-(0.050/pow(a,2)))); // Rotational friction scaling.
double fricPar = 2*M_PI/(log(a) - 0.207 + (0.980/a)-(0.133/pow(a,2)));
double fricPerp = 4*M_PI/(log(a) + 0.839 + (0.185/a)+(0.233/pow(a,2))); ;

// ND constants.
const int na = 10; // Number of swimmers
const int np = 10; // Number of colloids

// arrays.
double xa[na][3]; // Store for swimmer center positions.
double xp[np][2]; // Store for colloid positions.

#endif // 