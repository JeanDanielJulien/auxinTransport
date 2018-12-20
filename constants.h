#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

extern int step_number;
extern double stopping_tol;
extern int force_stop;

extern double pi;

extern double static_length;

extern double production_coeff_a;
extern double degradation_coeff_a;
extern double p_concentration;
extern double diffusion;

extern double transport_threshold;
extern double transport_power;

extern double feedback;

extern double feedback_power;
extern double feedback_slope;

extern double noise_amplitude;
extern int noise_amplitude_P;
extern int noise_amplitude_pa;

extern double min_stiffness;
extern double delta_stiffness;
extern double stiffness_power;
extern double stiffness_threshold;

extern double stress_threshold;

extern double stress_x;
extern double stress_y;
extern double init_width;
extern double init_height;

extern int ablated_cell;

#pragma omp threadprivate(step_number,static_length,degradation_coeff_a,p_concentration,diffusion,feedback_power,feedback_slope, transport_power, transport_threshold, stiffness_power, stiffness_threshold, min_stiffness, min_stiffness2, delta_stiffness, delta_stiffness2, feedback, noise_amplitude, force_stop,stress_x,stress_y,ablated_cell)

#endif // CONSTANTS_H_INCLUDED
