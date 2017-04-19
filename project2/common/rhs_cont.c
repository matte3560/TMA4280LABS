/* Continuous right hand side function */

#include "common_poisson.h"

#define M_E (2.7182818284590452353602874713527)

double poisson_rhs(double x, double y)
{
    return pow(M_E, x)*sin(2*PI*x)*sin(2*PI*y);
}
