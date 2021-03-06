#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "save_data.h"
#include <math.h>
#include "constants.h"
#include "stiffness.h"
#include "insertion_rate.h"

/* this function saves all the variables for further analysis */

void save_data(FILE *f, int vertices_number, double *vertices, int *vertices_number_in_each_cell, int cell_number, int *cells_vertices, int n,
               double *auxin, int *period, double width, double height, double *cell_areas, double *lengths, double *target_lengths, int *cells_neighbours){
    
    int i,j;
    
    fprintf(f,"cell_number= %i;\n",cell_number);
    fprintf(f,"vertices_number= %i;\n",vertices_number);
    
    fprintf(f,"vertices=[");
    for (i=0;i<vertices_number;i++){
        fprintf(f," %.4f %.4f ",vertices[2*i],vertices[2*i+1]);}
    fprintf(f,"];\n");
    
    fprintf(f,"vertices_number_per_cell=[");
    for (i=0;i<cell_number+1;i++){fprintf(f,"%i ",vertices_number_in_each_cell[i]);}
    fprintf(f," ];\n");
    
    fprintf(f,"cells_vertices=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f," %i ",cells_vertices[j]);
        }
        //    fprintf(f,"\n");
    }
    fprintf(f,"];\n");
    
    fprintf(f,"period=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f,"%i %i ",period[2*j],period[2*j+1]);}}
    fprintf(f,"];\n");
    
    fprintf(f,"width=%f;\nheight=%f;\n",width,height);
    
    fprintf(f,"auxin=[");
    for (i=0;i<cell_number;i++){fprintf(f,"%.6f ",auxin[i]);}
    fprintf(f," ];\n");
    
    fprintf(f,"lengths=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f,"%.6f ",lengths[j]);
        }
    }
    fprintf(f," ];\n");
    
    
    fprintf(f,"target_lengths=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f,"%.4f ",target_lengths[j]);
        }
    }
    fprintf(f," ];\n");
    
    double strain[vertices_number_in_each_cell[cell_number]];
    double stress[vertices_number_in_each_cell[cell_number]];
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            strain[j]=(lengths[j]-target_lengths[j])/target_lengths[j];
            stress[j]=strain[j]*stiffness(auxin[i])/stress_threshold;
        }
    }
    
    fprintf(f,"strain=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f,"%.6f ",strain[j]);
        }
    }
    fprintf(f," ];\n");
    
    
    fprintf(f,"stress=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f,"%.6f ",stress[j]);
        }
    }
    fprintf(f," ];\n");
    
    double pin1_density[ vertices_number_in_each_cell[cell_number] ];
    double sum_activation, sum_lengths;
    for (i=0;i<cell_number;i++){
        // compute the PIN1 insertion rate
        for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            pin1_density[j]=insertion_rate( feedback*strain[j] + (1.-feedback)*stress[j] );
        
            pin1_density[j] *= (double)(i!=ablated_cell) * (double)(cells_neighbours[j]!=ablated_cell) ;
        }
        // normalize for PIN1 competition between the different walls
        sum_activation=0.; sum_lengths=0.;
        for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            sum_activation+=lengths[j]*pin1_density[j];
            sum_lengths+=lengths[j];
        }
        for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            pin1_density[j]/=1.+sum_activation/sum_lengths;
            pin1_density[j]*=p_concentration;
        }
    }
    
    fprintf(f,"pin1_density=[");
    for (i=0;i<cell_number;i++){
        for (j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            fprintf(f,"%.6f ",pin1_density[j]);
        }
    }
    fprintf(f," ];\n");
    
    
    double total_stress[vertices_number_in_each_cell[cell_number]];
    
    for (i=0;i<cell_number;i++){
        for (int j=vertices_number_in_each_cell[i];j<vertices_number_in_each_cell[i+1];j++){
            
            total_stress[j] = stress[j] + stress[cells_neighbours[j]] ;
            
        }
    }
    
    
    
    fprintf(f,"ablated_cell=%i;\n",ablated_cell);
    
}
