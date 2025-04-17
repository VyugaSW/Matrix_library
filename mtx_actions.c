#include "mtx_actions.h"
#include <math.h>
#include "mtx_logs.h"

struct matrix
{
    double *data; // data + w * i + j
    size_t w, h;    
};

int mtx_transpose(matrix *mtx) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in transpose");
        return 1;
    }
    
    if (mtx->w != mtx->h) {
        MTX_LOG_ERROR("Non-square matrix in transpose");
        return -1;
    }
    
    for (size_t i = 0; i < mtx->h; i++) {
        for (size_t j = i + 1; j < mtx->w; j++) {
            double tmp = mtx->data[i * mtx->w + j];
            mtx->data[i * mtx->w + j] = mtx->data[j * mtx->w + i];
            mtx->data[j * mtx->w + i] = tmp;
        }
    }

    MTX_LOG("Matrix transposed");
    return 0;
}

int mtx_swap_rows(matrix *mtx, size_t row1, size_t row2) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in row swap");
        return 1;
    }
    if (row1 >= mtx->h || row2 >= mtx->h) {
        MTX_LOG_ERROR("Invalid row indices in swap");
        return -1;
    }
    
    for (size_t j = 0; j < mtx->w; j++) {
        double tmp = mtx->data[row1 * mtx->w + j];
        mtx->data[row1 * mtx->w + j] = mtx->data[row2 * mtx->w + j];
        mtx->data[row2 * mtx->w + j] = tmp;
    }

    MTX_LOG("Rows were swapped");
    return 0;
}

int mtx_swap_cols(matrix *mtx, size_t col1, size_t col2) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in column swap");
        return 1;
    }
    if (col1 >= mtx->w || col2 >= mtx->w) {
        MTX_LOG_ERROR("Invalid column indices in swap");
        return -1;
    }
    
    for (size_t i = 0; i < mtx->h; i++) {
        double tmp = mtx->data[i * mtx->w + col1];
        mtx->data[i * mtx->w + col1] = mtx->data[i * mtx->w + col2];
        mtx->data[i * mtx->w + col2] = tmp;
    }

    MTX_LOG("Columns were swapped");
    return 0;
}

int mtx_row_mult(matrix *mtx, size_t row, double factor) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in row multiply");
        return 1;
    }
    if (row >= mtx->h) {
        MTX_LOG_ERROR("Invalid row index in multiply");
        return -1;
    }
    
    for (size_t j = 0; j < mtx->w; j++) {
        mtx->data[row * mtx->w + j] *= factor;
    }

    MTX_LOG("Row was multiplied");
    return 0;
}

int mtx_row_div(matrix *mtx, size_t row, double divisor) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in row division");
        return 1;
    }
    if (row >= mtx->h) {
        MTX_LOG_ERROR("Invalid row index in division");
        return -1;
    }
    if (fabs(divisor) < MTX_MIN_DIVISOR) {
        MTX_LOG_ERROR("Division by zero in row division");
        return -1;
    }
    
    for (size_t j = 0; j < mtx->w; j++) {
        mtx->data[row * mtx->w + j] /= divisor;
    }

    MTX_LOG("Row was divided");
    return 0;
}

int mtx_row_add(matrix *mtx, size_t target_row, size_t source_row, double factor) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in row addition");
        return 1;
    }
    if (target_row >= mtx->h || source_row >= mtx->h) {
        MTX_LOG_ERROR("Invalid row indices in addition");
        return -1;
    }
    
    for (size_t j = 0; j < mtx->w; j++) {
        mtx->data[target_row * mtx->w + j] += 
            mtx->data[source_row * mtx->w + j] * factor;
    }

    MTX_LOG("Rows were added");
    return 0;
}

double mtx_norm(const matrix *mtx) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in norm calculation");
        return -1.0;
    }
    
    double max_norm = 0.0;
    for (size_t i = 0; i < mtx->h; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < mtx->w; j++) {
            row_sum += fabs(mtx->data[i * mtx->w + j]);
        }
        if (row_sum > max_norm) {
            max_norm = row_sum;
        }
    }

    MTX_LOG("Matrix norm calculated");
    return max_norm;
}