#pragma once

#include <math.h>
#include <assert.h>

#define TESTPRINT(msg) printf("\n\n=====  %s  =====\n", #msg)
#define TESTPASS() puts("\n...PASSED")
#define APPROX(a,b) assert( fabs((a)-(b)) < 0.0001 )


void test_print_mat(double** mat, int rows, int cols);
void test_print_vec(double* vec, int size);

void test_cmp_vec(double* v1, double* v2, int size);
void test_cmp_mat(double** m1, double** m2, int rows, int cols);
