#pragma once

void serial_grid(double *grid, double h, int n);
void serial_diag(double *diag, int n, int m);
void serial_gen_rhs(double **b, double *grid, double h, int m);
void serial_transpose(double **bt, double **b, size_t m);
void serial_dst(double** b, int n, int m, bool inv);
