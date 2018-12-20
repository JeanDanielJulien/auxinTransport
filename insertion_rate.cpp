#include "insertion_rate.h"
#include "constants.h"
#include <math.h>

// insertion rate of PIN1 proteins, as a function of the stress or the strain

double insertion_rate(double mechanics){

return pow(feedback_slope*mechanics,feedback_power);

}
