#pragma once

#include <stdio.h>

/**
 * @brief Matrix structure (row-major storage)
 */
struct matrix; 
typedef struct matrix matrix;


/* ================== Core Allocation ================== */

/**
 * @brief Allocates uninitialized matrix
 * @param w Number of columns
 * @param h Number of rows
 * @return Pointer to allocated matrix, NULL on failure
 */
matrix* mtx_alloc(size_t w, size_t h);

int mtx_assign(matrix *mtx1, const matrix *mtx2);

/**
 * @brief Creates a deep copy of matrix
 * @param mtx Source matrix to copy
 * @return New matrix copy, NULL on failure
 */
matrix* mtx_copy(const matrix* mtx);

/**
 * @brief Releases matrix memory
 * @param m Matrix to deallocate (safe with NULL)
 */
void mtx_free(matrix* m);

/* ================== Accessors ================== */

/**
 * @brief Gets mutable pointer to element
 * @param mtx Target matrix
 * @param i Row index (0-based)
 * @param j Column index (0-based)
 * @return Pointer to element, NULL on invalid indices
 */
double* mtx_ptr(matrix* mtx, size_t i, size_t j);

/**
 * @brief Gets const pointer to element
 * @param mtx Target matrix
 * @param i Row index (0-based)
 * @param j Column index (0-based)
 * @return Const pointer to element, NULL on invalid indices
 */
const double* mtx_cptr(const matrix* mtx, size_t i, size_t j);

/* ================== Initialization ================== */

/**
 * @brief Sets all elements to zero
 * @param mtx Matrix to initialize
 */
void mtx_set_zero(matrix* mtx);

/**
 * @brief Converts matrix to identity (must be square)
 * @param mtx Matrix to modify
 */
void mtx_set_id(matrix* mtx);

/**
 * @brief Creates zero-initialized matrix
 * @param w Number of columns
 * @param h Number of rows
 * @return New matrix, NULL on failure
 */
matrix* mtx_alloc_zero(size_t w, size_t h);

/**
 * @brief Creates identity matrix
 * @param w Matrix dimension (must be square)
 * @return New identity matrix, NULL on failure
 */
matrix* mtx_alloc_id(size_t w, size_t h);

/* ================== I/O Operations ================== */

/**
 * @brief Reads matrix from stdin
 * @param mtx Matrix to populate
 * @return 0 on success, -1 on error
 */
int mtx_input(matrix* mtx);

/**
 * @brief Prints matrix to stdout
 * @param mtx Matrix to print
 * @param precision Decimal places to show
 */
void mtx_print(const matrix* mtx, int precision); 

size_t mtx_get_width(const matrix *mtx);

size_t mtx_get_height(const matrix *mtx);

