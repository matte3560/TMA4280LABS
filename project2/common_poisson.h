#pragma once

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define PI 3.14159265358979323846

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
extern void fst_(double *v, int *n, double *w, int *nn);
extern void fstinv_(double *v, int *n, double *w, int *nn);

double poisson_rhs(double x, double y);
double *mk_1D_array(size_t n, bool zero);
double **mk_2D_array(size_t n1, size_t n2, bool zero);
void free_2D_array(double **array);
