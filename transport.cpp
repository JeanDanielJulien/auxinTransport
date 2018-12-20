#include "transport.h"
#include "constants.h"
#include <math.h>

// compute the saturated auxin transport as a function of auxin concentration

double transport(double a){

return pow(a,transport_power)/(pow(a,transport_power)+pow(transport_threshold,transport_power));

}
