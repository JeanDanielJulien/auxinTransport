#ifndef SAVE_DATA_H_INCLUDED
#define SAVE_DATA_H_INCLUDED

void save_data(FILE *f, int vertices_number, double *vertices, int *vertices_number_in_each_cell, int cell_number, int *cells_vertices, int n,
               double *auxin, int *period, double width, double height, double *cell_areas, double *lengths, double *target_lengths, int *cells_neighbours);

#endif // SAVE_DATA_H_INCLUDED
