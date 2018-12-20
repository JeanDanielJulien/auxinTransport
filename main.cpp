#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "load_initial_conditions.h"
#include "initialize.h"
#include "step.h"
#include "constants.h"
#include <math.h>
#include <fstream>
#include "compute_areas.h"
#include "compute_cell_areas.h"
#include "minimize_BFGS.h"
#include "compute_cell_areas.h"
#include <dirent.h>
#include "save_data.h"
#include <nlopt.h>
#include "stability_x.h"
#include "stability_y.h"
#include <unistd.h>

using namespace std;

int main(){
    
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        printf("Current working directory: %s\n", cwd);
    }
    
    srand((unsigned int)time(NULL));
    
    // path to inital conditions
    char initial_file_name[256];
    //sprintf(initial_file_name,"/Users/jjulien/Desktop/codes/auxinTransportStressStrain/initialConditions/1600_noise0_1.txt");
    sprintf(initial_file_name,"/Path/to/initial/conditions/1600_noise0_1.txt");
    
    // path to output files and basis of the output files names
    char file_name_basis[256];
    //sprintf(file_name_basis,"/Users/jjulien/Desktop/codes/auxinTransportStressStrain/output/name_of_the_simulation");
    sprintf(file_name_basis,"/Path/to/output/folder/name_of_the_simulation");
    
    
    
    
    /*  allocate memory and initialize tissue topology */
    
    int cell_number, vertices_number, i, k;
    
    // read the number of cells and vertices to allocate memory
    std::ifstream file(initial_file_name,ios::in);
    file >> cell_number;
    file >> vertices_number;
    file.close();
    
    double * vertices; vertices=(double*)calloc(2*vertices_number,sizeof(double));
    double * width; width=(double*)calloc(1,sizeof(double));
    double * height; height=(double*)calloc(1,sizeof(double));
    int * vertices_number_in_each_cell; vertices_number_in_each_cell=(int*)calloc(cell_number+1,sizeof(int));
    vertices_number_in_each_cell[0]=0;
    for (i=1;i<=cell_number;i++){vertices_number_in_each_cell[i]=vertices_number_in_each_cell[i-1]+6;}
    int * cells_vertices; cells_vertices=(int*)calloc(vertices_number_in_each_cell[cell_number],sizeof(int));
    int * period; period=(int*)calloc(2*vertices_number_in_each_cell[cell_number],sizeof(int));
    int * cells_neighbours; cells_neighbours=(int*)calloc(vertices_number_in_each_cell[cell_number],sizeof(int));
    
    // initialize tissue topology
    load_initial_conditions(initial_file_name, cell_number, vertices_number, vertices, vertices_number_in_each_cell, cells_vertices, period, cells_neighbours, width, height);
    
    printf("Initial condition loaded from %s\n",initial_file_name);
    
    /*  allocate memory for tissue geometry, mechanics, and auxin concentration */
    
    double * lengths; lengths=(double*)calloc(vertices_number_in_each_cell[cell_number],sizeof(double));
    double * target_lengths; target_lengths=(double*)calloc(vertices_number_in_each_cell[cell_number],sizeof(double));
    double * areas; areas=(double*)calloc(vertices_number_in_each_cell[cell_number],sizeof(double));
    double * cell_areas; cell_areas=(double*)calloc(cell_number, sizeof(double)); /* areas of cells */
    double * forces; forces=(double*)calloc(2*vertices_number, sizeof(double)); /* forces on vertices */
    double * auxin; auxin=(double*)calloc(cell_number, sizeof(double)); /* auxin concentration */

    // some noise in the position of the vertices if needed
    for (i=0;i<2*vertices_number;i++){vertices[i]=(1.+0.*(-1.+2.*(double)rand()/(double)RAND_MAX))*vertices[i];}
    
    /* from the vertices, initialize tissue geometry */
    initialize(vertices,vertices_number, vertices_number_in_each_cell, cells_vertices,cell_number,
               lengths,target_lengths,areas,cell_areas,forces,auxin,period, width[0], height[0]);
    
    
    
    printf("First energy minimization\n");
    // Uniform dilation to speed up the first minimization
    double dilatation=1.+(sqrt(3.)/2.)*((stress_x+stress_y)/2.)/(min_stiffness+delta_stiffness*pow(stiffness_threshold,stiffness_power)/(pow(stiffness_threshold,stiffness_power)+pow(production_coeff_a/degradation_coeff_a,stiffness_power)));
    
    init_width=width[0]; init_height=height[0];
    
    for (i=0;i<2*vertices_number;i++){vertices[i]*=dilatation;}
    width[0]*=dilatation; height[0]*=dilatation;
    
    minimize_BFGS(vertices_number, cell_number, vertices, vertices_number_in_each_cell, cells_vertices, areas, cell_areas, lengths, target_lengths, auxin, period, width, height);
    
    // compute the initial deformation of the tissue
    static_length=0.;
    for (i=0; i<vertices_number_in_each_cell[cell_number];i++){ static_length+=lengths[i];}
    static_length/=vertices_number_in_each_cell[cell_number];
    
    
    k=0; double noise[cell_number]; double noise_P[cell_number]; double noise_pa[cell_number];
    
    printf("Saving simulation parameters\n");
    FILE *f=NULL;
    char file_name[256];
    sprintf(file_name,"%s_parameters.m",file_name_basis);
    f=fopen(file_name,"w+");
    fprintf(f,"degradation_coeff_a=%f;\np_concentration=%f;\ntransport_threshold=%f;\ntransport_power=%f;\nfeedback_power=%f;\nfeedback_slope=%f;\nmin_stiffness=%f;\ndelta_stiffness=%f;\nstiffness_power=%f;\nstiffness_threshold=%f;\nstatic_length=%f;\nfeedback=%f;\nstress_x=%f;\nstress_y=%f;\nnoise_amplitude=%f;\ndiffusion=%f;\n",
            degradation_coeff_a,p_concentration,transport_threshold,transport_power,feedback_power,feedback_slope,min_stiffness,delta_stiffness,stiffness_power,stiffness_threshold,static_length,feedback,stress_x,stress_y,noise_amplitude,diffusion);
    
    // Compute the variables from the linear stability analysis
    double M[3],stable_feedback;
    stable_feedback=static_length/1.-1.;
    stable_feedback*=feedback_slope;
    if (feedback==0){
        stable_feedback *= min_stiffness+delta_stiffness*pow(stiffness_threshold,stiffness_power)/(pow(stiffness_threshold,stiffness_power)+pow(production_coeff_a/degradation_coeff_a,stiffness_power));
        stable_feedback /= stress_threshold;
    }
    stable_feedback=pow(stable_feedback,feedback_power);
    stable_feedback=stable_feedback/(1.+stable_feedback);
    M[0]=(2./(3.*sqrt(3.)))*(1/static_length)*p_concentration*stable_feedback;
    M[0]*=transport_power*pow(transport_threshold,transport_power)*pow(production_coeff_a/degradation_coeff_a,transport_power-1.);
    M[0]/=pow(pow(transport_threshold,transport_power)+pow(production_coeff_a/degradation_coeff_a,transport_power),2.);
    M[0]+=diffusion*(2./(3.*sqrt(3.)))*(1/static_length);
    M[1]=(2./(3.*sqrt(3.)))*(1/static_length)*p_concentration*stable_feedback;
    M[1]*=pow(production_coeff_a/degradation_coeff_a,transport_power)/(pow(transport_threshold,transport_power)+pow(production_coeff_a/degradation_coeff_a,transport_power));
    M[1]*=feedback_power;
    M[1]*=stiffness_power*delta_stiffness*pow(stiffness_threshold,stiffness_power)*pow(production_coeff_a/degradation_coeff_a,stiffness_power-1.);
    M[1]/=pow(pow(stiffness_threshold,stiffness_power)+pow(production_coeff_a/degradation_coeff_a,stiffness_power),2.);
    M[1]/=min_stiffness+delta_stiffness*pow(stiffness_threshold,stiffness_power)/(pow(stiffness_threshold,stiffness_power)+pow(production_coeff_a/degradation_coeff_a,stiffness_power));
    M[2]=M[1]*stable_feedback;
    M[2]/=2.;
    fprintf(f,"A=%f;\nB=%f;\nC=%f;\n",M[0],M[1],M[2]);
    
    struct parameters_stability {double A; double B; double C;};
    double eigenvalue,kx[1]={0.},ky[1]={0.};
    struct parameters_stability parameters;
    parameters.A=M[0]; parameters.B=M[1]; parameters.C=M[2];
    void* pparameters=&parameters;
    
    nlopt_opt opt;
    double lb[1]={0.}; double ub[1]={pi};
    
    // Find the most unstable wave number in the x-direction
    printf("Linear stability analysis, first minimization\n");
    opt=nlopt_create(NLOPT_LN_COBYLA,1);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_max_objective(opt, stability_x, pparameters);
    nlopt_optimize(opt, kx, &eigenvalue);
    if (eigenvalue>0){ fprintf(f,"unstable_x=%d;\n",1); fprintf(f,"L_x=%f;\n",2*pi/kx[0]);}
    else {  fprintf(f,"unstable_x=%d;\n",0); fprintf(f,"L_x=%f;\n",0.);}
    nlopt_destroy(opt);
    
    ub[0]=2.*pi/sqrt(3.);
    
    // Find the most unstable wave number in the y-direction
    printf("Linear stability analysis, second minimization\n");
    opt=nlopt_create(NLOPT_LN_COBYLA,1);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_max_objective(opt, stability_y, pparameters);
    nlopt_optimize(opt, ky, &eigenvalue);
    if (eigenvalue>0){ fprintf(f,"unstable_y=%d;\n",1); fprintf(f,"L_y=%f;\n",2*pi/ky[0]);}
    else {  fprintf(f,"unstable_y=%d;\n",0); fprintf(f,"L_y=%f;\n",0.);}
    nlopt_destroy(opt);
    
    fclose(f);
    
    printf("Writing done\n");
    
    // Initial noise, speeds up the few first steps
    for (i=0; i<cell_number;i++){noise[i]=1.+noise_amplitude*(2.*(double)rand()/(double)RAND_MAX-1.);}
    
    // Save the initial state of the tissue, before auxin dynamics
    FILE *finit=NULL; char file_name_init[256];
    sprintf(file_name_init,"%s_initial_state.m",file_name_basis);
    finit=fopen(file_name_init,"w+");
    save_data(finit, vertices_number, vertices, vertices_number_in_each_cell, cell_number, cells_vertices, k, auxin, period, width[0], height[0], cell_areas,lengths, target_lengths, cells_neighbours);
    fclose(finit);
    
    
    /* loop over the number of iteration */
    while (k<=step_number && force_stop==0){
        
        // update the noise
        for (i=0; i<cell_number;i++){noise_P[i]=1.+(double)(noise_amplitude_P)*(2.*(double)rand()/(double)RAND_MAX-1.);}
        for (i=0; i<cell_number;i++){noise_pa[i]=1.+(double)(noise_amplitude_pa)*(2.*(double)rand()/(double)RAND_MAX-1.);}
        
        // Perform one time step
        printf("Step %d\n",k);
        step(vertices, vertices_number, vertices_number_in_each_cell, cells_vertices, cell_number,
             lengths, target_lengths, areas, cell_areas, forces,
             auxin, noise, k, period, cells_neighbours, width, height,
             noise_P, noise_pa);
        k++;
        
        // Stop the initial noise
        if (k==10){
            for (i=0; i<cell_number;i++){noise[i]=1.;}
            noise_amplitude=0.;
        }
        
    }
    
    // Saves the variables
    FILE *f2=NULL; char file_name2[256];
    sprintf(file_name2,"%s_final_state.m", file_name_basis);
    
    f2=fopen(file_name2,"w+");
    save_data(f2, vertices_number, vertices, vertices_number_in_each_cell, cell_number, cells_vertices, k, auxin, period, width[0], height[0], cell_areas,lengths, target_lengths, cells_neighbours);
    fclose(f2);
    
    free(vertices); free(cells_vertices); free(width); free(height);
    free(lengths); free(areas); free(cell_areas); free(forces); free(auxin);
    free(target_lengths); free(period);
    free(cells_neighbours);
    
    printf("fin du programme\n");
    return 0;
}
