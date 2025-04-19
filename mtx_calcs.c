#include "mtx_repmem.h"
#include "mtx_arithmetic.h"
#include "mtx_actions.h"
#include "mtx_logs.h"
#include <math.h>
#include <string.h>

struct matrix
{
    double *data; // data + w * i + j
    size_t w, h;    
};

matrix *mtx_exp(const matrix *mtx, double eps) {
    if(eps <= MTX_MIN_DIVISOR) {
        MTX_LOG_ERROR("Epsilon cant equal to zero");
        return NULL;
    }

    if(!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid matrix in exp calculation");
        return NULL;
    }

    matrix *res = mtx_alloc_id(mtx->w, mtx->h);
    if(!res) return NULL;

    matrix *term = mtx_copy(mtx);
    if(!term) {
        mtx_free(res);
        return NULL;
    }

    matrix *mtx_pow = mtx_copy(mtx);
    if(!mtx_pow) {  
        mtx_free(res);
        mtx_free(term);
        return NULL;
    }

    double factorial = 1.0;
    int k = 1;

    while(mtx_norm(term) >= eps) {
        if(mtx_add(res, term) != 0) break;

        factorial *= ++k;
        if(mtx_mul(mtx_pow, mtx) != 0) break;
        
        matrix* new_term = mtx_copy(mtx_pow);
        if(!new_term) break;
        
        if(mtx_sdiv(new_term, factorial) != 0) {
            mtx_free(new_term);
            break;
        }
        
        if(mtx_assign(term, new_term) != 0) {
            mtx_free(new_term);
            break;
        }
        mtx_free(new_term);
    }

    mtx_free(term);
    mtx_free(mtx_pow);
    return res;
}

matrix* mtx_solve_gauss(const matrix* A, const matrix* B){
    // Check inputs
    if (!A || !B || !A->data || !B->data) {
        MTX_LOG_ERROR("Null matrix in solver");
        return NULL;
    }
    if (A->w != A->h) {
        MTX_LOG_ERROR("Matrix A must be square");
        return NULL;
    }
    if (A->h != B->h) {
        MTX_LOG_ERROR("Dimension mismatch between A and B");
        return NULL;
    }

    const size_t n = A->h;
    const size_t m = B->w;

    // Create augmented matrix [A|B]
    matrix* aug = mtx_alloc(n + m, n);
    if (!aug) {
        MTX_LOG_ERROR("Failed to allocate augmented matrix");
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        memcpy(aug->data + (n+m)*i, A->data + n*i, sizeof(double)*n);
        memcpy(aug->data + (n+m)*i + n, B->data + m*i, sizeof(double)*m);
    }

    // Gaussian elimination with partial pivoting
    for (size_t k = 0; k < n; ++k) {
        // Partial pivoting: find row with maximum element in current column
        size_t max_row = k;
        double max_val = fabs(*mtx_cptr(aug, k, k));
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(*mtx_cptr(aug, i, k));
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        // Swap rows if necessary
        if (max_row != k) {
            if (mtx_swap_rows(aug, k, max_row) != 0) {
                mtx_free(aug);
                MTX_LOG_ERROR("Row swap failed");
                return NULL;
            }
        }

        // Check for zero pivot (matrix is singular)
        if (fabs(*mtx_cptr(aug, k, k)) < MTX_MIN_DIVISOR) {
            mtx_free(aug);
            MTX_LOG_ERROR("Matrix is singular (zero pivot)");
            return NULL;
        }

        // Eliminate below and above the current row
        for (size_t j = 0; j < n; j++) {
            if (j != k) {
                double factor = *mtx_cptr(aug, j, k) / *mtx_cptr(aug, k, k);
                if (mtx_row_add(aug, j, k, -factor) != 0) {
                    mtx_free(aug);
                    MTX_LOG_ERROR("Eliminating rows operation failed");
                    return NULL;
                }
            }
        }

        // Normalize the current row
        double pivot = *mtx_cptr(aug, k, k);
        if (mtx_row_div(aug, k, pivot) != 0) {
            mtx_free(aug);
            MTX_LOG_ERROR("Row normalization failed");
            return NULL;
        }
        
    }

    matrix* X = mtx_alloc(m, n);
    if (!X) {
        mtx_free(aug);
        MTX_LOG_ERROR("Failed to allocate solution matrix");
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < B->w; j++) {
            *mtx_ptr(X, i, j) = *mtx_cptr(aug, i, n + j);
        }
    }

    mtx_free(aug);
    return X;
}

double mtx_verify_solution(const matrix* A, const matrix* X, const matrix* B) {
    if (!A || !X || !B || A->w != X->h || X->w != B->w || A->h != B->h) {
        MTX_LOG_ERROR("Invalid dimensions in solution verification");
        return -1.0;
    }

    matrix* AX = mtx_alloc(A->h, X->w);
    if (!AX) {
        MTX_LOG_ERROR("Failed to allocate temp matrix");
        return -1.0;
    }

    if (mtx_mul2(AX, A, X) != 0) {
        mtx_free(AX);
        MTX_LOG_ERROR("Matrix multiplication failed");
        return -1.0;
    }

    if (mtx_sub(AX, B) != 0) {
        mtx_free(AX);
        MTX_LOG_ERROR("Matrix subtraction failed");
        return -1.0;
    }

    double residual = mtx_norm(AX);
    mtx_free(AX);
    return residual;
}