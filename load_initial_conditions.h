#ifndef LOAD_INITIAL_CONDITIONS_H_INCLUDED
#define LOAD_INITIAL_CONDITIONS_H_INCLUDED

void load_initial_conditions(const char* file_name, int cell_number, int vertices_number, double *vertices, int *vertices_number_in_each_cell, int *cells_vertices, int *period, int *cells_neighbours, double *width, double *height);

#endif // LOAD_INITIAL_CONDITIONS_H_INCLUDED
