#include "constants.h"
#include "auxin_dynamics.h"
#include <gsl/gsl_errno.h>
#include <math.h>
#include "stiffness.h"
#include "insertion_rate.h"
#include "transport.h"


/* this function computes the dynamics of auxin transport  */

int auxin_dynamics(double t, const double x[], double dx[], void *parameters){
    
    struct parameters_list2 {int cell_number; int *vertices_number_in_each_cell; int *cells_vertices; double *lengths; double *cell_areas; double *noise; double *strain; int *cells_neighbours; double *noise_P; double *noise_pa;};
    
    int i; double exchange_ik;
    
    double activation[ ((parameters_list2*)parameters)->vertices_number_in_each_cell[((parameters_list2*)parameters)->cell_number] ];
    double sum_activation, sum_lengths;
    for (i=0;i<(((parameters_list2*)parameters)->cell_number);i++){
        // compute the PIN1 insertion rate
        for (int j=((parameters_list2*)parameters)->vertices_number_in_each_cell[i];j<((parameters_list2*)parameters)->vertices_number_in_each_cell[i+1];j++){
            activation[j]=(((parameters_list2*)parameters)->strain)[j];
            activation[j]*=feedback+(1.-feedback)*stiffness(x[i])/stress_threshold;
            activation[j]=insertion_rate(activation[j]);
        }
        
        // check for ablated cell
        for (int j=((parameters_list2*)parameters)->vertices_number_in_each_cell[i];j<((parameters_list2*)parameters)->vertices_number_in_each_cell[i+1];j++){
            activation[j] *= (1.-(double)(i==ablated_cell)) * (1.-(double)(((parameters_list2*)parameters)->cells_neighbours[j]==ablated_cell)) ;
        }
        // normalize for PIN1 competition between the different walls
        sum_activation=0.; sum_lengths=0.;
        for (int j=((parameters_list2*)parameters)->vertices_number_in_each_cell[i];j<((parameters_list2*)parameters)->vertices_number_in_each_cell[i+1];j++){
            sum_activation+=(((parameters_list2*)parameters)->lengths)[j]*activation[j];
            sum_lengths+=(((parameters_list2*)parameters)->lengths)[j];
        }
        for (int j=((parameters_list2*)parameters)->vertices_number_in_each_cell[i];j<((parameters_list2*)parameters)->vertices_number_in_each_cell[i+1];j++){
            activation[j]/=1.+sum_activation/sum_lengths;
            activation[j]*=p_concentration*(((parameters_list2*)parameters)->noise_P)[i];
        }
    }
    // production/degradation
    for (i=0;i<(((parameters_list2*)parameters)->cell_number);i++){
        dx[i]=((((parameters_list2*)parameters)->noise)[i]*(((parameters_list2*)parameters)->noise_pa)[i]*production_coeff_a-degradation_coeff_a*x[i]);
    }
    
    // diffusion and transport
    for (i=0;i<(((parameters_list2*)parameters)->cell_number);i++){
        
        for (int j=((parameters_list2*)parameters)->vertices_number_in_each_cell[i];j<((parameters_list2*)parameters)->vertices_number_in_each_cell[i+1];j++){
            
            exchange_ik=activation[j]*transport(x[i]);
            exchange_ik+=diffusion*x[i];
            exchange_ik*=(1.-(double)(i==ablated_cell))*(1.-(double)(((parameters_list2*)parameters)->cells_neighbours[j]==ablated_cell));;
            
            dx[i]-=(((parameters_list2*)parameters)->lengths)[j]*exchange_ik/(((parameters_list2*)parameters)->cell_areas)[i];
            dx[((parameters_list2*)parameters)->cells_neighbours[j]]+=(((parameters_list2*)parameters)->lengths)[j]*exchange_ik/(((parameters_list2*)parameters)->cell_areas)[((parameters_list2*)parameters)->cells_neighbours[j]];
            
        }
    }
    
    return GSL_SUCCESS;
}
