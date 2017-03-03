#pragma once

#include <stdint.h>

typedef struct {
	double value;
	double reduc_time;
} result_t;

double zeta_term( int32_t i );
result_t zeta_sum( int32_t n );
result_t zeta_pi( int32_t n );

