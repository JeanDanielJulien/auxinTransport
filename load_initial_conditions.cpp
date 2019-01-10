#include "load_initial_conditions.h"
#include <fstream>
#include <iostream>
#include <limits>

using namespace std;

// this function reads the file containing the initial conditions

void load_initial_conditions(const char* file_name, int cell_number, int vertices_number, int cumulated_vertices_number, double *vertices, int *vertices_number_in_each_cell, int *cells_vertices, int *period, int* cells_neighbours, double *width, double *height){
   
    std::ifstream file(file_name,ios::in);
    
    file >> cell_number;
    file >> vertices_number;
    file >> cumulated_vertices_number;
    
    for (int i=0;i<2*vertices_number;i++){
        file >> vertices[i];
    }
    
    for (int i=0;i<=cell_number;i++){
        file >> vertices_number_in_each_cell[i];
        }
    
    for (int i=0;i<vertices_number_in_each_cell[cell_number];i++){
        file >> cells_vertices[i];
    }
    
    for (int i=0;i<2*vertices_number_in_each_cell[cell_number];i++){
        file >> period[i];
    }
    
    for (int i=0;i<vertices_number_in_each_cell[cell_number];i++){
        file >> cells_neighbours[i];
    }
    
    file >> width[0];
    file >> height[0];
    
    file.close();
    
}
