#pragma once

#include "common_poisson.h"

#include <mpi.h>


/* MPI related globals */
MPI_Comm mpi_comm;
int mpi_rank;
int mpi_size;

double mpi_poisson(int n);

/* Commented out -> not implemented yet */
int mpi_init(int *argc, char ***argv);
int mpi_finalize();
void mpi_grid(double *grid, double h, int n);
void mpi_diag(double *diag, int m, int n);
void mpi_gen_rhs(double **b, double *grid, double h, int m);
void mpi_transpose(double **bt, double **b, int m);
void mpi_dst(double **b, int m, int n, bool inv);
//void mpi_solve_tu(double **b, double *diag, int m);
//double mpi_u_max(double **b, int m);

/* Array related functions */
size_t mpi_idx_start(const size_t size);
size_t mpi_idx_end(const size_t size);
size_t mpi_idx_size(const size_t size);
size_t mpi_idx_part(const size_t size);
size_t mpi_padded_size(const size_t size);
void mpi_allgather_mat(double** mat, int m, int n);
void mpi_allgather_vec(double* vec, int size);
