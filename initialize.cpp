#include "compute_areas.h"
#include "compute_cell_areas.h"
#include "compute_lengths.h"
#include "initialize.h"
#include "constants.h"
#include <stdio.h>

/* this function initializes the various lists used in the program */

void initialize(double *vertices, int vertices_number, int *vertices_number_in_each_cell, int *cells_vertices, int cell_number,
                double *lengths, double *target_lengths, double *areas, double *cell_areas, double *forces,
                double *auxin, int *period, double width, double height){
    
    int i;
    compute_areas(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, areas,period, width, height);
    compute_cell_areas(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, areas, cell_areas);
    compute_lengths(vertices, vertices_number_in_each_cell,cells_vertices, cell_number, lengths, period, width, height);
    for (i=0;i<vertices_number_in_each_cell[cell_number];i++){ target_lengths[i]=lengths[i];}
    for (i=0;i<cell_number;i++){ auxin[i]=production_coeff_a/degradation_coeff_a;
    }
    
}
