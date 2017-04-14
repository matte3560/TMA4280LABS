#pragma once

#include "common_poisson.h"

double serial_poisson(int n);

void serial_grid(double *grid, double h, int n);
void serial_diag(double *diag, int m, int n);
void serial_gen_rhs(double **b, double *grid, double h, int m);
void serial_transpose(double **bt, double **b, int m);
void serial_dst(double **b, int m, int n, bool inv);
void serial_solve_tu(double **b, double *diag, int m);
double serial_u_max(double **b, int m);
