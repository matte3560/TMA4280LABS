#pragma once

#include "mpi_poisson.h"

#include <mpi.h>


poisson_result_t hybrid_poisson(int n);

void hybrid_grid(double *grid, double h, int n);
void hybrid_diag(double *diag, int m, int n);
void hybrid_gen_rhs(double **b, double *grid, double h, int m);
void hybrid_dst(double **b, int m, int n, bool inv);
void hybrid_solve_tu(double **b, double *diag, int m);
