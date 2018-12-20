#include "stiffness.h"
#include "constants.h"
#include <math.h>

// compute the stiffness of the cell walls as a function of auxin concentration 

double stiffness(double a){

return min_stiffness+delta_stiffness*pow(stiffness_threshold,stiffness_power)/(pow(a,stiffness_power)+pow(stiffness_threshold,stiffness_power));

}
