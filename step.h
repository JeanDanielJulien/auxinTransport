#ifndef STEP_H_INCLUDED
#define STEP_H_INCLUDED

void step(double *vertices, int vertices_number, int *vertices_number_in_each_cell, int *cells_vertices, int cell_number,
          double *lengths, double *target_lengths, double *areas, double *cell_areas, double *forces,
          double *auxin, double *noise, int step_index, int *period, int *cells_neighbours, double *width, double *height, double *noise_P, double *noise_pa);

#endif // STEP_H_INCLUDED
