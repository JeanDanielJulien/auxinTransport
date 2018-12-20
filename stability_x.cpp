#include "constants.h"
#include "stability_x.h"
#include <math.h>
#include <stdio.h>
#include <nlopt.h>

/* this compute the rate of the linear stability analysis, as a function of the wave number in the x-direction */

double stability_x(unsigned n, const double *kx, double *grad, void *parameters){

struct parameters_stability {double A; double B; double C;};

double Sx=(1./3.)*(2.*cos(sqrt(3.)*kx[0]/2.)+1);
double eigenvalue;
eigenvalue=-degradation_coeff_a-6.*((parameters_stability*)parameters)->A+6.*(1.-feedback)*((parameters_stability*)parameters)->B+6.*(2.*feedback-1.)*((parameters_stability*)parameters)->C+
            (6.*((parameters_stability*)parameters)->A-6.*(1.-feedback)*((parameters_stability*)parameters)->B+12.*(1.-feedback)*((parameters_stability*)parameters)->C)*Sx-
            6.*Sx*Sx*((parameters_stability*)parameters)->C;

return(eigenvalue);

}
