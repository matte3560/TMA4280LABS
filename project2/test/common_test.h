#pragma once

#include <math.h>

#define TESTPRINT(msg) printf("\n\n=====  %s  =====\n", #msg)
#define TESTPASS() puts("\n...PASSED")
#define APPROX(a,b) assert( fabs((a)-(b)) < 0.0001 )


void test_print_mat(double** mat, int size);
void test_print_vec(double* vec, int size);
