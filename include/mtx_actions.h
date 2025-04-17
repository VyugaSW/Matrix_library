#pragma once

#include "mtx_repmem.h"
#include "mtx_arithmetic.h"

/* ================== Matrix Transformations ================== */

/**
 * @brief Transposes a square matrix in-place
 * @param mtx Matrix to transpose (must be square)
 * @return 0 on success, 1 if NULL pointer, -1 if non-square matrix
 */
int mtx_transpose(matrix *mtx);

/* ================== Row/Column Operations ================== */

/**
 * @brief Swaps two rows in a matrix
 * @param mtx Matrix to modify
 * @param row1 First row index (0-based)
 * @param row2 Second row index (0-based)
 * @return 0 on success, 1 if NULL pointer, -1 if invalid row indices
 */
int mtx_swap_rows(matrix *mtx, size_t row1, size_t row2);

/**
 * @brief Swaps two columns in a matrix
 * @param mtx Matrix to modify
 * @param col1 First column index (0-based)
 * @param col2 Second column index (0-based)
 * @return 0 on success, 1 if NULL pointer, -1 if invalid column indices
 */
int mtx_swap_cols(matrix *mtx, size_t col1, size_t col2);

/**
 * @brief Multiplies a row by a scalar factor
 * @param mtx Matrix to modify
 * @param row Row index to multiply (0-based)
 * @param factor Multiplication factor
 * @return 0 on success, 1 if NULL pointer, -1 if invalid row index
 */
int mtx_row_mult(matrix *mtx, size_t row, double factor);

/**
 * @brief Divides a row by a scalar value
 * @param mtx Matrix to modify
 * @param row Row index to divide (0-based)
 * @param divisor Divisor (operation skipped if |divisor| < 1e-20)
 * @return 0 on success, 1 if NULL pointer, -1 if invalid row index or division by zero
 */
int mtx_row_div(matrix *mtx, size_t row, double divisor);

/**
 * @brief Adds a multiple of one row to another
 * @param mtx Matrix to modify
 * @param target_row Row to be modified (0-based)
 * @param source_row Row to multiply and add (0-based)
 * @param factor Multiplication factor for source row
 * @return 0 on success, 1 if NULL pointer, -1 if invalid row indices
 */
int mtx_row_add(matrix *mtx, size_t target_row, size_t source_row, double factor);

/* ================== Matrix Norms ================== */

/**
 * @brief Computes the infinity norm (maximum row sum) of a matrix
 * @param mtx Input matrix
 * @return Matrix norm, -1.0 if NULL pointer or invalid matrix
 */
double mtx_norm(const matrix *mtx);
