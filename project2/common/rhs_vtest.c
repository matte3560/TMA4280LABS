#include "common_poisson.h"

/*
 * Right hand side with known exact solution
 * u(x,y) = sin(PI*x)*sin(2*PI*y)
 */

double poisson_rhs(double x, double y)
{
    return 5*PI*PI*sin(PI*x)*sin(2*PI*y);
}
