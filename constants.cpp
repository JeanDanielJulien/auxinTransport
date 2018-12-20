#include <math.h>

int step_number=10001;
double stopping_tol=1.e-3;
int force_stop=0;

double pi=3.1415926;

double static_length=0.;

double production_coeff_a=1.;
double degradation_coeff_a=0.5e0;
double p_concentration=100.;
double diffusion=5.;

double transport_threshold=2.;
double transport_power=1.;

double feedback=0; //0 -> stress   1 -> strain

double feedback_power=3.;
double feedback_slope=30.;

double noise_amplitude=0.1;
int noise_amplitude_P=0.;
int noise_amplitude_pa=0.;

double min_stiffness=1.e6;
double delta_stiffness=4.e6;
double stiffness_power=2.;
double stiffness_threshold=2.;

double stress_threshold=3.e6;

double stress_x=1.4e5;
double stress_y=1.4e5;
double init_width=0.;
double init_height=0.;

int ablated_cell=-1;
