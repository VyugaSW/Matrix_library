#pragma once

#include "mtx_repmem.h"
#include "mtx_arithmetic.h"
#include "mtx_actions.h"
#include "mtx_logs.h"

/**
 * @brief Computes matrix exponential e^A using Taylor series expansion
 * @param mtx Input matrix to compute exponential of
 * @param eps Precision threshold for stopping the series expansion
 * @return Pointer to newly allocated matrix containing result,
 *         NULL on error or if input is invalid
 * @note Uses Taylor series: e^A = I + A + A?/2! + A?/3! + ...
 * 
 * The function continues adding terms until the norm of the next term
 * is smaller than eps (||A^k/k!|| < eps)
 */
matrix *mtx_exp(const matrix *mtx, double eps);


/**
 * @brief Solves linear system AX = B using Gaussian elimination with partial pivoting
 * @param A Square matrix (n x n)
 * @param B Right-hand side matrix (n x m)
 * @return Solution matrix X (n x m), NULL on failure
 * 
 * @note det(A) != 0
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