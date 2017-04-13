#pragma once

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define PI 3.14159265358979323846

double poisson_rhs(double x, double y);
double *mk_1D_array(size_t n, bool zero);
double **mk_2D_array(size_t n1, size_t n2, bool zero);
