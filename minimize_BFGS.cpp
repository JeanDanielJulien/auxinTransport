#include "compute_energy_and_forces.h"
#include "minimize_BFGS.h"
#include "compute_areas.h"
#include "compute_cell_areas.h"
#include "compute_lengths.h"
#include <nlopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

/* this function minimizes the mechanical energy with the BFGS algorithm (nlopt library necessary) */

void minimize_BFGS(int vertices_number, int cell_number, double *vertices, int *vertices_number_in_each_cell, int *cells_vertices, double *areas, double *cell_areas,
                   double *lengths, double *target_lengths, double *auxin, int *period, double *width, double *height){
    
    
    compute_areas(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, areas,period, width[0], height[0]);
    compute_cell_areas(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, areas, cell_areas);
    compute_lengths(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, lengths, period, width[0], height[0]);
    
    
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LD_LBFGS,2*vertices_number+2);
    
    nlopt_set_stopval(opt, -HUGE_VAL);
    nlopt_set_maxeval(opt, -1.);
    nlopt_set_maxtime(opt, -1.);
    nlopt_set_xtol_rel(opt, -1.);
    nlopt_set_xtol_abs1(opt, 1.e-6);
    nlopt_set_ftol_rel(opt, -1.);
    nlopt_set_ftol_abs(opt, -1.);
    
    // structure for the parameters
    struct parameters_list {int vertices_number; int cell_number; int *vertices_number_in_each_cell; int *cells_vertices; double *areas; double *cell_areas;
        double *lengths; double *target_lengths; double *auxin; int *period;};
    
    struct parameters_list parameters;
    parameters.vertices_number=vertices_number;
    parameters.cell_number=cell_number;
    parameters.vertices_number_in_each_cell=vertices_number_in_each_cell;
    parameters.cells_vertices=cells_vertices;
    parameters.areas=areas;
    parameters.cell_areas=cell_areas;
    parameters.lengths=lengths;
    parameters.target_lengths=target_lengths;
    parameters.auxin=auxin;
    parameters.period=period;
    void* pparameters=&parameters;
    
    // define function to minimize
    int reso=nlopt_set_min_objective(opt, compute_energy_and_forces, pparameters);
    
    // define variables of the function: vertices positions, width and height of the tissue
    double vertices_BFGS[2*vertices_number+2];
    for (int i=0; i<2*vertices_number;i++){ vertices_BFGS[i]=vertices[i];}
    vertices_BFGS[2*vertices_number]=width[0]; vertices_BFGS[2*vertices_number+1]=height[0];
    
    // actual minimization
    double min_energy;
    nlopt_optimize(opt, vertices_BFGS, &min_energy);
    printf("Optimization returns %i, size of the tissue from [%.4f x %.4f] to [%.4f x %.4f]  \n",reso,width[0],height[0],vertices_BFGS[2*vertices_number],vertices_BFGS[2*vertices_number+1]);
    
    for (int i=0; i<2*vertices_number;i++){vertices[i]=vertices_BFGS[i];}
    width[0]=vertices_BFGS[2*vertices_number]; height[0]=vertices_BFGS[2*vertices_number+1];
    
    nlopt_destroy(opt);
    
    // update geometry
    compute_areas(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, areas,period, width[0], height[0]);
    compute_cell_areas(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, areas, cell_areas);
    compute_lengths(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, lengths, period, width[0], height[0]);
    
}
