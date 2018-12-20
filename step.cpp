#include "compute_cell_areas.h"
#include "compute_areas.h"
#include "compute_lengths.h"
#include "step.h"
#include <stdio.h>
#include "constants.h"
#include <math.h>
#include "minimize_BFGS.h"
#include "auxin_dynamics.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "compute_energy_and_forces.h"
#include "save_data.h"
#include <algorithm>

/* this function performs a time step (mechanics, then auxin dynamics) */

void step(double *vertices, int vertices_number, int *vertices_number_in_each_cell, int *cells_vertices, int cell_number,
          double *lengths, double *target_lengths, double *areas, double *cell_areas, double *forces,
          double *auxin, double *noise, int step_index, int *period, int *cells_neighbours, double *width, double *height, double *noise_P, double *noise_pa){
    
    double old_areas[cell_number];
    for (int i=0; i<cell_number; i++){
        old_areas[i]=cell_areas[i];
        
    }
    
    /* minimization of mechanical energy  */
    minimize_BFGS(vertices_number, cell_number, vertices, vertices_number_in_each_cell, cells_vertices, areas, cell_areas, lengths, target_lengths,auxin, period, width, height);
    
    // update geometry
    compute_areas(vertices, vertices_number_in_each_cell, cells_vertices, cell_number, areas, period, width[0], height[0]);
    compute_cell_areas(vertices, vertices_number_in_each_cell, cells_vertices, cell_number, areas, cell_areas);
    compute_lengths(vertices, vertices_number_in_each_cell, cells_vertices, cell_number, lengths, period, width[0], height[0]);
    
    // update auxin concentration due to dilution
    for (int i=0; i<cell_number; i++){
        auxin[i]*=old_areas[i]/cell_areas[i];
    }
    
    // compute strain
    double strain[vertices_number_in_each_cell[cell_number]];
    for (int i=0; i<cell_number; i++){
        for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            strain[j]=(lengths[j]-target_lengths[j])/target_lengths[j];
        }
    }
    
    // parameters structure for auxin dynamics
    struct parameters_list2 {int cell_number; int *vertices_number_in_each_cell; int *cells_vertices; double *lengths; double *cell_areas;
        double *noise; double *strain; int *cells_neighbours; double *noise_P; double *noise_da;};
    struct parameters_list2 parameters2;
    parameters2.cell_number=cell_number;
    parameters2.vertices_number_in_each_cell=vertices_number_in_each_cell;
    parameters2.cells_vertices=cells_vertices;
    parameters2.lengths=lengths;
    parameters2.cell_areas=cell_areas;
    parameters2.noise=noise;
    parameters2.strain=strain;
    parameters2.cells_neighbours=cells_neighbours;
    parameters2.noise_P=noise_P;
    parameters2.noise_da=noise_pa;
    void *pparameters2=&parameters2;
    int system_size=cell_number;
    
    double x[system_size];
    for (int i=0;i<cell_number;i++){ x[i]=auxin[i];}
    int s; double t=0.0;
    
    // time-step in auxin dynamics
    gsl_odeiv2_system sys={auxin_dynamics, NULL, static_cast<size_t>(system_size), pparameters2};
    gsl_odeiv2_driver *d=gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk2, 1.e-3, 1.e-3, 0.);
    s=gsl_odeiv2_driver_apply(d, &t, 1.e-1, x);
    if (s){printf("problem in rk2\n");}
    for (int i=0;i<cell_number;i++){ auxin[i]=x[i];}
    gsl_odeiv2_driver_free(d);
    
    
    
}





