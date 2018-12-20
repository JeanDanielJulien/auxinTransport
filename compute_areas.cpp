#include <math.h>
#include "constants.h"
#include "compute_areas.h"
#include <stdio.h>
#include <stdlib.h>

/* this function computes the areas of the triangles made with two vertices of the same cell and the barycenter of this cell
   the area is given by half the cross product of the centered coordinates of the vertices */

void compute_areas(double *vertices, int *vertices_number_in_each_cell, int *cells_vertices,
                    int cell_number, double *areas, int *period, double width, double height){

int i,j,previous;

/* usual loop over the cells, then over the vertices of each cell */
for (i=0;i<cell_number;i++){
    
    previous=vertices_number_in_each_cell[i+1]-1;
    for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
        areas[j]=0.5*((vertices[2*cells_vertices[previous]]+width*(double)period[2*previous])*(vertices[2*cells_vertices[j]+1]+height*(double)period[2*j+1])-
                    (vertices[2*cells_vertices[previous]+1]+height*(double)period[2*previous+1])*(vertices[2*cells_vertices[j]]+width*(double)period[2*j]));
        previous=j;
    }
}

}
