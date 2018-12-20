#include "constants.h"
#include "stability_y.h"
#include <math.h>
#include <stdio.h>
#include <nlopt.h>

/* this compute the rate of the linear stability analysis, as a function of the wave number in the y-direction */

double stability_y(unsigned n, const double *ky, double *grad, void *parameters){

struct parameters_stability {double A; double B; double C;};

double Sy=(1./3.)*(2.*cos(ky[0]/2.)+cos(ky[0]));
double eigenvalue;
eigenvalue=-degradation_coeff_a-6.*((parameters_stability*)parameters)->A+6.*(1.-feedback)*((parameters_stability*)parameters)->B+6.*(2.*feedback-1.)*((parameters_stability*)parameters)->C+
            (6.*((parameters_stability*)parameters)->A-6.*(1.-feedback)*((parameters_stability*)parameters)->B+12.*(1.-feedback)*((parameters_stability*)parameters)->C)*Sy-
            6.*Sy*Sy*((parameters_stability*)parameters)->C;

return(eigenvalue);

}
