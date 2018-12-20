#include "compute_energy_and_forces.h"
#include "constants.h"
#include "compute_areas.h"
#include "compute_cell_areas.h"
#include "compute_lengths.h"
#include <math.h>
#include <stdio.h>
#include <nlopt.h>
#include "stiffness.h"

/* this function computes the mechanical energy and its derivatives */

double compute_energy_and_forces(unsigned n, const double *vertices, double *forces, void *parameters){
    
    // define the parameters structure
    struct parameters_list {int vertices_number; int cell_number; int *vertices_number_in_each_cell; int *cells_vertices; double *areas; double *cell_areas;
        double *lengths; double *target_lengths; double *auxin; int *period;};
    
    int j,previous; double energy=0,ablation;
    
    // update geometry
    compute_lengths((double*)vertices, ((parameters_list*)parameters)->vertices_number_in_each_cell, ((parameters_list*)parameters)->cells_vertices, ((parameters_list*)parameters)->cell_number,
                    ((parameters_list*)parameters)->lengths, ((parameters_list*)parameters)->period, (double)vertices[2*(((parameters_list*)parameters)->vertices_number)], (double)vertices[2*(((parameters_list*)parameters)->vertices_number)+1]);
    compute_areas((double*)vertices, ((parameters_list*)parameters)->vertices_number_in_each_cell, ((parameters_list*)parameters)->cells_vertices, ((parameters_list*)parameters)->cell_number,
                  ((parameters_list*)parameters)->areas, ((parameters_list*)parameters)->period, (double)vertices[2*(((parameters_list*)parameters)->vertices_number)], (double)vertices[2*(((parameters_list*)parameters)->vertices_number)+1]);
    compute_cell_areas((double*)vertices, ((parameters_list*)parameters)->vertices_number_in_each_cell, ((parameters_list*)parameters)->cells_vertices,
                       ((parameters_list*)parameters)->cell_number, ((parameters_list*)parameters)->areas, ((parameters_list*)parameters)->cell_areas);
    
    // check that forces are asked and compute them
    if (forces){ for (int i=0;i<2*((parameters_list*)parameters)->vertices_number+2;i++){forces[i]=0.;}}
    
    
    for (int i=0;i<((parameters_list*)parameters)->cell_number;i++){
        ablation=1;
        if (i==ablated_cell){ablation=0.;}
        previous=((parameters_list*)parameters)->vertices_number_in_each_cell[i+1]-1;
        for (j=((parameters_list*)parameters)->vertices_number_in_each_cell[i];j<((parameters_list*)parameters)->vertices_number_in_each_cell[i+1];j++){
            
            /* derivatives of the lengths of the walls */
            if (forces){
                
                forces[2*(((parameters_list*)parameters)->cells_vertices)[j]]+=ablation*stiffness(((parameters_list*)parameters)->auxin[i])*
                ((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j])*
                (vertices[2*(((parameters_list*)parameters)->cells_vertices)[j]]+vertices[2*((parameters_list*)parameters)->vertices_number]*(double)(((parameters_list*)parameters)->period)[2*j]
                 -vertices[2*(((parameters_list*)parameters)->cells_vertices)[previous]]-vertices[2*((parameters_list*)parameters)->vertices_number]*(double)(((parameters_list*)parameters)->period)[2*previous])/
                ((((parameters_list*)parameters)->lengths)[j]*(((parameters_list*)parameters)->target_lengths)[j]);
                
                
                forces[2*(((parameters_list*)parameters)->cells_vertices)[previous]]+=ablation*stiffness(((parameters_list*)parameters)->auxin[i])*
                ((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j])*
                (vertices[2*(((parameters_list*)parameters)->cells_vertices)[previous]]+vertices[2*((parameters_list*)parameters)->vertices_number]*(double)(((parameters_list*)parameters)->period)[2*previous]
                 -vertices[2*(((parameters_list*)parameters)->cells_vertices)[j]]-vertices[2*((parameters_list*)parameters)->vertices_number]*(double)(((parameters_list*)parameters)->period)[2*j])/
                ((((parameters_list*)parameters)->lengths)[j]*(((parameters_list*)parameters)->target_lengths)[j]);
                
                
                forces[2*(((parameters_list*)parameters)->cells_vertices)[j]+1]+=ablation*stiffness(((parameters_list*)parameters)->auxin[i])*
                ((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j])*
                (vertices[2*(((parameters_list*)parameters)->cells_vertices)[j]+1]+vertices[2*((parameters_list*)parameters)->vertices_number+1]*(double)(((parameters_list*)parameters)->period)[2*j+1]
                 -vertices[2*(((parameters_list*)parameters)->cells_vertices)[previous]+1]-vertices[2*((parameters_list*)parameters)->vertices_number+1]*(double)(((parameters_list*)parameters)->period)[2*previous+1])/
                ((((parameters_list*)parameters)->lengths)[j]*(((parameters_list*)parameters)->target_lengths)[j]);
                
                
                forces[2*(((parameters_list*)parameters)->cells_vertices)[previous]+1]+=ablation*stiffness(((parameters_list*)parameters)->auxin[i])*
                ((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j])*
                (vertices[2*(((parameters_list*)parameters)->cells_vertices)[previous]+1]+vertices[2*((parameters_list*)parameters)->vertices_number+1]*(double)(((parameters_list*)parameters)->period)[2*previous+1]
                 -vertices[2*(((parameters_list*)parameters)->cells_vertices)[j]+1]-vertices[2*((parameters_list*)parameters)->vertices_number+1]*(double)(((parameters_list*)parameters)->period)[2*j+1])/
                ((((parameters_list*)parameters)->lengths)[j]*(((parameters_list*)parameters)->target_lengths)[j]);
                
                
                
                
                forces[2*((parameters_list*)parameters)->vertices_number]+=ablation*(double)((((parameters_list*)parameters)->period)[2*j]-(((parameters_list*)parameters)->period)[2*previous])*
                stiffness(((parameters_list*)parameters)->auxin[i])*
                ((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j])*
                (vertices[2*(((parameters_list*)parameters)->cells_vertices)[j]]+vertices[2*((parameters_list*)parameters)->vertices_number]*(double)(((parameters_list*)parameters)->period)[2*j]
                 -vertices[2*(((parameters_list*)parameters)->cells_vertices)[previous]]-vertices[2*((parameters_list*)parameters)->vertices_number]*(double)(((parameters_list*)parameters)->period)[2*previous])/
                ((((parameters_list*)parameters)->lengths)[j]*(((parameters_list*)parameters)->target_lengths)[j]);
                
                
                forces[2*((parameters_list*)parameters)->vertices_number+1]+=ablation*(double)((((parameters_list*)parameters)->period)[2*j+1]-(((parameters_list*)parameters)->period)[2*previous+1])*
                stiffness(((parameters_list*)parameters)->auxin[i])*
                ((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j])*
                (vertices[2*(((parameters_list*)parameters)->cells_vertices)[j]+1]+vertices[2*((parameters_list*)parameters)->vertices_number+1]*(double)(((parameters_list*)parameters)->period)[2*j+1]
                 -vertices[2*(((parameters_list*)parameters)->cells_vertices)[previous]+1]-vertices[2*((parameters_list*)parameters)->vertices_number+1]*(double)(((parameters_list*)parameters)->period)[2*previous+1])/
                ((((parameters_list*)parameters)->lengths)[j]*(((parameters_list*)parameters)->target_lengths)[j]);
            }
            /* walls energy */
            energy+=ablation*0.5*stiffness(((parameters_list*)parameters)->auxin[i])*pow((((parameters_list*)parameters)->lengths)[j]-(((parameters_list*)parameters)->target_lengths)[j],2.)/(((parameters_list*)parameters)->target_lengths)[j];
            
            previous=j;
            
        }
    }
    
    // tension applied to the tissue
    if (forces){
        forces[2*(((parameters_list*)parameters)->vertices_number)]-=stress_x*init_height;
        forces[2*(((parameters_list*)parameters)->vertices_number)+1]-=stress_y*init_width;
    }
    
    energy-=stress_x*vertices[2*(((parameters_list*)parameters)->vertices_number)]*init_height
    +stress_y*vertices[2*(((parameters_list*)parameters)->vertices_number)+1]*init_width;
    
    return(energy);
    
}