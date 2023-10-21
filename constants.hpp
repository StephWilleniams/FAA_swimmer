#ifndef constants
#define constants

#include "header.hpp"

// Note: Captialized variables are dimensional.

// Time constants (s).
double T = 0.1; // Total runtime (s).
double DT = 0.0005; // Step size (s) (ND time). 
// Time constants ND.
double dT = DT/DT;
int nSteps = int(T/DT); // Number of calculated steps.

// Length constants (um).
double XR = 10;  // System right side.
double XL = -10; // System left side.
double YT = 10;  // System top side.
double YB = -10; // System bottom side.
double RS = 1; // Swimmer size (ND length)
// Length constants ND.
double xR = XR/RS;
double xL = XL/RS;
double yT = YT/RS;
double yB = YB/RS;
double rS = RS/RS;

// Velocity constants (um/s).
double VS = 100; // Swimmer speed.
// Velocity constants ND.
double vS = VS*DT/RS;

// Diffusion constants (xx^2/s).
double DR = 1; // Swimmer rotational diffusivity, rad^2/s.
// Diffusion constants ND.
double dR = DR*DT; 
double sigR = sqrt(dR);

// ND constants
const int ns = 2; // Number of swimmers

// arrays.
double x[ns][3]; // Store for particle positions.

#endif // 