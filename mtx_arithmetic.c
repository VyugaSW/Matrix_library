#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mtx_arithmetic.h"
#include "mtx_logs.h"

struct matrix
{
    double *data; // data + w * i + j
    size_t w, h;    
};

int mtx_add(matrix *mtx1, const matrix *mtx2) {
    if (!mtx1 || !mtx2 || !mtx1->data || !mtx2->data) {
        MTX_LOG_ERROR("Null matrix pointer in addition");
        return 1;
    }
    if (mtx1->w != mtx2->w || mtx1->h != mtx2->h) {
        MTX_LOG_ERROR("Matrix size mismatch in addition");
        return -1;
    }

    for (size_t i = 0; i < mtx1->h * mtx1->w; i++) {
        mtx1->data[i] += mtx2->data[i];
    }
    MTX_LOG("Matrix addition completed");
    return 0;
}

int mtx_sub(matrix *mtx1, const matrix *mtx2) {
    if (!mtx1 || !mtx2 || !mtx1->data || !mtx2->data) {
        MTX_LOG_ERROR("Null matrix pointer in subtraction");
        return 1;
    }
    if (mtx1->w != mtx2->w || mtx1->h != mtx2->h) {
        MTX_LOG_ERROR("Matrix size mismatch in subtraction");
        return -1;
    }

    for (size_t i = 0; i < mtx1->h * mtx1->w; i++) {
        mtx1->data[i] -= mtx2->data[i];
    }
    MTX_LOG("Matrix subtraction completed");
    return 0;
}

void mtx_smul(matrix *mtx, double d) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in scalar multiplication");
        return;
    }
    
    for (size_t i = 0; i < mtx->h * mtx->w; i++) {
        mtx->data[i] *= d;
    }
    MTX_LOG("Matrix scalar multiplication completed");
}

void mtx_sdiv(matrix *mtx, double d) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Null matrix in scalar division");
        return;
    }
    if (fabs(d) < 1e-20) {
        MTX_LOG_ERROR("Attempt to divide by zero in scalar division");
        return;
    }
    
    mtx_smul(mtx, 1.0/d);
    MTX_LOG("Matrix scalar division completed");
}

int mtx_add2(matrix *mtx, const matrix *mtx1, const matrix *mtx2) {
    if (!mtx || !mtx1 || !mtx2 || !mtx->data || !mtx1->data || !mtx2->data) {
        MTX_LOG_ERROR("Null matrix pointer in add2 operation");
        return 1;
    }
    if (mtx1->w != mtx2->w || mtx1->h != mtx2->h) {
        MTX_LOG_ERROR("Source matrices size mismatch in add2");
        return -1;
    }
    if (mtx->w != mtx1->w || mtx->h != mtx1->h) {
        MTX_LOG_ERROR("Destination matrix size mismatch in add2");
        return -1;
    }

    for (size_t i = 0; i < mtx->h * mtx->w; i++) {
        mtx->data[i] = mtx1->data[i] + mtx2->data[i];
    }
    MTX_LOG("Matrix add2 operation completed");
    return 0;
}

int mtx_sub2(matrix *mtx, const matrix *mtx1, const matrix *mtx2) {
    if (!mtx || !mtx1 || !mtx2 || !mtx->data || !mtx1->data || !mtx2->data) {
        MTX_LOG_ERROR("Null matrix pointer in sub2 operation");
        return 1;
    }
    if (mtx1->w != mtx2->w || mtx1->h != mtx2->h) {
        MTX_LOG_ERROR("Source matrices size mismatch in sub2");
        return -1;
    }
    if (mtx->w != mtx1->w || mtx->h != mtx1->h) {
        MTX_LOG_ERROR("Destination matrix size mismatch in sub2");
        return -1;
    }

    for (size_t i = 0; i < mtx->h * mtx->w; i++) {
        mtx->data[i] = mtx1->data[i] - mtx2->data[i];
    }
    MTX_LOG("Matrix sub2 operation completed");
    return 0;
}

int mtx_smul2(matrix *mtx, const matrix *mtx1, double d) {
    if (!mtx || !mtx1 || !mtx->data || !mtx1->data) {
        MTX_LOG_ERROR("Null matrix pointer in smul2 operation");
        return 1;
    }
    if (mtx->w != mtx1->w || mtx->h != mtx1->h) {
        MTX_LOG_ERROR("Matrix size mismatch in smul2");
        return -1;
    }

    for (size_t i = 0; i < mtx->h * mtx->w; i++) {
        mtx->data[i] = mtx1->data[i] * d;
    }
    MTX_LOG("Matrix smul2 operation completed");
    return 0;
}

int mtx_sdiv2(matrix *mtx, const matrix *mtx1, double d) {
    if (!mtx || !mtx1 || !mtx->data || !mtx1->data) {
        MTX_LOG_ERROR("Null matrix pointer in sdiv2 operation");
        return 1;
    }
    if (fabs(d) < 1e-20) {
        MTX_LOG_ERROR("Attempt to divide by zero in sdiv2");
        return 1;
    }
    if (mtx->w != mtx1->w || mtx->h != mtx1->h) {
        MTX_LOG_ERROR("Matrix size mismatch in sdiv2");
        return -1;
    }

    mtx_smul2(mtx, mtx1, 1.0/d);

    MTX_LOG("Matrix sdiv2 operation completed");
    return 0;
}

int mtx_mul (matrix *mtx1, const matrix *mtx2){
    if(!mtx1 || !mtx1->data || !mtx2 || !mtx2->data){
        MTX_LOG_ERROR("Null matrix pointer in mul operation");
        return 1;
    }

    if(mtx1->w != mtx2->h){
        MTX_LOG_ERROR("Matrixes different sizes for mul");
        return -1;
    }

    matrix *temp = mtx_alloc(mtx1->w, mtx2->h);
    if(!temp){
        MTX_LOG_ERROR("Allocation in mul is failure");
        return -1;
    }

    for (size_t i = 0; i < mtx1->h; i++) {
        for (size_t j = 0; j < mtx2->w; j++) {
            double sum = 0;
            for (size_t k = 0; k < mtx1->w; k++) {
                sum += (*mtx_cptr(mtx1,i,k)) * (*mtx_cptr(mtx2,k,j));
            }
            *mtx_ptr(temp,i,j) = sum;
        }
    }

    mtx_free(mtx1);
    mtx1 = mtx_copy(temp);
    mtx_free(temp);

    return 0;
}

int mtx_mul2 (matrix *mtx, const matrix *mtx1, const matrix *mtx2){
    if(!mtx1 || !mtx1->data || !mtx2 || !mtx2->data || !mtx || !mtx->data){
        MTX_LOG_ERROR("Null matrix pointer in mul2 operation");
        return 1;
    }

    if(mtx1->w != mtx2->h || mtx->w != mtx1->w || mtx->h != mtx2->h){
        MTX_LOG_ERROR("Incompatible matrix sizes for mul2");
        return -1;
    }


    matrix *temp = NULL;
    matrix *result = NULL;

    if(mtx == mtx1 || mtx == mtx2) {
        temp = mtx_alloc(mtx2->w, mtx1->h);
        if(!temp) {
            MTX_LOG_ERROR("Allocation for temp matrix failed");
            return -1;
        }
        result = temp;
    } 
    else {
        if(mtx->w != mtx2->w || mtx->h != mtx1->h) {
            MTX_LOG_ERROR("Result matrix has wrong dimensions");
            return -1;
        }
        result = mtx;
    }

    for(size_t i = 0; i < mtx1->h; i++) {
        for(size_t j = 0; j < mtx2->w; j++) {
            double sum = 0;
            for(size_t k = 0; k < mtx1->w; k++) {
                sum += (*mtx_cptr(mtx1, i, k)) * (*mtx_cptr(mtx2, k, j));
            }
            *mtx_ptr((matrix*)result, i, j) = sum;
        }
    }

    if(temp){
        free(mtx->data);
        mtx_assign(mtx, temp);
        mtx_free(temp);
    }

    return 0;
}