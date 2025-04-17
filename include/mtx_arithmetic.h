#pragma once

#include "mtx_repmem.h"

/* ================== In-place Operations ================== */

/**
 * @brief Performs matrix addition (mtx1 += mtx2)
 * @param mtx1 Target matrix to modify
 * @param mtx2 Source matrix to add
 * @return 0 on success, 1 if any pointer is NULL, -1 if size mismatch
 */
int mtx_add(matrix *mtx1, const matrix *mtx2);

/**
 * @brief Performs matrix subtraction (mtx1 -= mtx2)
 * @param mtx1 Target matrix to modify
 * @param mtx2 Source matrix to subtract
 * @return 0 on success, 1 if any pointer is NULL, -1 if size mismatch
 */
int mtx_sub(matrix *mtx1, const matrix *mtx2);

/**
 * @brief Scales matrix by scalar (mtx *= d)
 * @param mtx Matrix to scale
 * @param d Scaling factor
 */
void mtx_smul(matrix *mtx, double d);

/**
 * @brief Divides matrix by scalar (mtx /= d)
 * @param mtx Matrix to divide
 * @param d Divisor (operation skipped if |d| < 1e-10)
 */
void mtx_sdiv(matrix *mtx, double d);

/* ================== Result-storing Operations ================== */

/**
 * @brief Computes matrix addition (mtx = mtx1 + mtx2)
 * @param mtx Destination matrix
 * @param mtx1 First operand
 * @param mtx2 Second operand
 * @return 0 on success, 1 if any pointer is NULL, -1 if size mismatch
 */
int mtx_add2(matrix *mtx, const matrix *mtx1, const matrix *mtx2);

/**
 * @brief Computes matrix subtraction (mtx = mtx1 - mtx2)
 * @param mtx Destination matrix
 * @param mtx1 First operand
 * @param mtx2 Second operand
 * @return 0 on success, 1 if any pointer is NULL, -1 if size mismatch
 */
int mtx_sub2(matrix *mtx, const matrix *mtx1, const matrix *mtx2);

/**
 * @brief Computes scalar multiplication (mtx = mtx1 * d)
 * @param mtx Destination matrix
 * @param mtx1 Source matrix
 * @param d Scaling factor
 * @return 0 on success, 1 if any pointer is NULL, -1 if size mismatch
 */
int mtx_smul2(matrix *mtx, const matrix *mtx1, double d);

/**
 * @brief Computes scalar division (mtx = mtx1 / d)
 * @param mtx Destination matrix
 * @param mtx1 Source matrix
 * @param d Divisor (operation skipped if |d| < 1e-10)
 * @return 0 on success, 1 if any pointer is NULL, -1 if size mismatch
 */
int mtx_sdiv2(matrix *mtx, const matrix *mtx1, double d);

int mtx_mul (matrix *m1, const matrix *m2);

int mtx_mul2 (matrix *m, const matrix *m1, const matrix *m2);