#pragma once

#include "mtx_repmem.h"
#include "mtx_arithmetic.h"
#include "mtx_actions.h"
#include "mtx_logs.h"

matrix *mtx_exp(const matrix *mtx, double eps);


/**
 * @brief Solves linear system AX = B using Gaussian elimination with partial pivoting
 * @param A Square matrix (n x n)
 * @param B Right-hand side matrix (n x m)
 * @return Solution matrix X (n x m), NULL on failure
 */
matrix* mtx_solve_gauss(const matrix* A, const matrix* B);

/**
 * @brief Verifies solution by computing residual ||AX - B||
 * @param A Input matrix
 * @param X Solution matrix
 * @param B Right-hand side
 * @return Residual norm, -1.0 on error
 */
double mtx_verify_solution(const matrix* A, const matrix* X, const matrix* B);