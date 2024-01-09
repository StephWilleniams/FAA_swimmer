#ifndef constants
#define constants

#include "header.hpp"

// Note: Captialized variables are dimensional (lower case ND).

// Time constants (s).
double T = 1; // Total runtime (s).
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
// Particle segement sizes (ND)
double rSeg[3] = {(4*sqrt(2)/(5*sqrt(2)+2))*rA,(4*sqrt(2)/(5*sqrt(2)+2))*rA/sqrt(2),(4*sqrt(2)/(5*sqrt(2)+2))*rA/2};
double  sDV[3] = {(rSeg[0] - rA), (rSeg[0]*2 - rA),(rSeg[0]*2 + rSeg[1] - rA)};

// Velocity constants (um/s).
double VA = 100; // Swimmer speed.
// Velocity constants ND.
double vA = VA*DT/RA;

// Diffusion constants (xx^2/s).
double DR = 0; // Swimmer rotational diffusivity, rad^2/s.
double DTherm = 0.1; // Passive particle diffusivity, um^2/s.
// Diffusion constants ND.
double dR = DR*DTherm; 
double dTherm = DTherm*DTherm/pow(RA,2); 
double sigR = sqrt(dR);
double sigT = sqrt(dT);

// ND constants
const int ns = 2; // Number of swimmers
const int np = 0; // Number of colloids
const int polarization = 1; // Active particle polarisation. +1 for pusher, -1 for puller.

// arrays.
double xs[ns][3]; // Store for swimmer center positions.
double xp[ns][2]; // Store for colloid positions.

#endif // 