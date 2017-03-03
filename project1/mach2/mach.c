#include "mach.h"

#include <math.h>


double mach_term( double x, int32_t i )
{
	return pow( -1, i-1 ) * ( pow( x, 2*i-1 ) / ( 2*i-1 ) );
}

double mach_sum( double x, int32_t n )
{
	double sum = 0;

#pragma omp parallel for reduction(+:sum)
	for ( int32_t i = 1; i <= n; i++ )
	{
		sum += mach_term( x, i );
	}

	return sum;
}

double mach_pi( int32_t n )
{
	return 4 * ( 4*mach_sum( 1.0/5.0, n ) - mach_sum( 1.0/239.0, n ) );
}

