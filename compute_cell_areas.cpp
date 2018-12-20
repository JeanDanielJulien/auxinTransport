#include <math.h>
#include "constants.h"
#include "compute_cell_areas.h"
#include <stdio.h>
#include <stdlib.h>

/* this function computes the areas of each cell, knowing the areas of each triangles between two vertices and the barycenter */

void compute_cell_areas(double *vertices, int *vertices_number_in_each_cell, int *cells_vertices,
                    int cell_number, double *areas, double *cell_areas){

/* usual loop over the cells, then over the vertices of each cell */
for (int i=0;i<cell_number;i++){cell_areas[i]=0;
    for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
        cell_areas[i]+=areas[j];
        }
}

}
