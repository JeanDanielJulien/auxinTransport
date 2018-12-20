#include <math.h>
#include "constants.h"
#include "compute_lengths.h"
#include <stdio.h>

/* this functions computes the lengths of each wall */

void compute_lengths(double *vertices, int *vertices_number_in_each_cell, int *cells_vertices, int cell_number, double *lengths, int *period, double width, double height){
    
    
    int previous;
    
    for (int i=0;i<cell_number;i++){
        previous=vertices_number_in_each_cell[i+1]-1;
        for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            
            lengths[j]=sqrt(pow(vertices[2*cells_vertices[j]]+width*(double)period[2*j]-vertices[2*cells_vertices[previous]]-width*(double)period[2*previous],2.)+
                            pow(vertices[2*cells_vertices[j]+1]+height*(double)period[2*j+1]-vertices[2*cells_vertices[previous]+1]-height*(double)period[2*previous+1],2.));
            previous=j;
        }
    }
    
}
