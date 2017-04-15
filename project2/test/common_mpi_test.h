#pragma once

#include <unistd.h>

#include "common_test.h"

#define MPI_TESTPRINT(msg) do{ if (mpi_rank == 0) TESTPRINT(msg); }while(false)
#define MPI_TESTPASS() do{ if (mpi_rank == 0) TESTPASS(); }while(false)

#define MPI_RANK0(code) do{ if (mpi_rank == 0) { code } }while(false)
#define MPI_SYNC(code) \
	do{ \
		MPI_Barrier(mpi_comm); \
		sleep(1); \
		do{ code }while(false); \
		sleep(1); \
		MPI_Barrier(mpi_comm); \
	}while(false)
