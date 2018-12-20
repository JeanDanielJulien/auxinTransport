#ifndef INITIALIZE_H_INCLUDED
#define INITIALIZE_H_INCLUDED

void initialize(double *vertices, int vertices_number, int *vertices_number_in_each_cell, int *cells_vertices, int cell_number,
                double *lengths, double *target_lengths, double *areas, double *cell_areas, double *forces,
                double *auxin, int *period, double width, double height);

#endif // INITIALIZE_H_INCLUDED
