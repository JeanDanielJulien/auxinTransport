#ifndef MINIMIZE_BFGS_H_INCLUDED
#define MINIMIZE_BFGS_H_INCLUDED

void minimize_BFGS(int vertices_number, int cell_number, double *vertices, int *vertices_number_in_each_cell, int *cells_vertices, double *areas, double *cell_areas,
                   double *lengths, double *target_lengths, double *auxin, int *period, double *width, double *height);

#endif // MINIMIZE_BFGS_H_INCLUDED
