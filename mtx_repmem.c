#include "mtx_repmem.h"
#include <stdlib.h>
#include <string.h>
#include "mtx_logs.h"

struct matrix
{
    double *data; // data + w * i + j
    size_t w, h;    
};


matrix* mtx_alloc(size_t w, size_t h) {
    if(w == 0 || h == 0) {
        MTX_LOG_ERROR("Attempt to allocate matrix with zero dimensions");
        return NULL;
    }

    matrix *mtx = malloc(sizeof(matrix));
    if (!mtx) {
        MTX_LOG_ERROR("Failed to allocate matrix struct");
        return NULL; 
    }

    mtx->data = (double*)malloc(w * h * sizeof(double));
    if (!mtx->data) {
        MTX_LOG_ERROR("Failed to allocate matrix data");
        free(mtx); 
        return NULL; 
    }

    mtx->w = w;
    mtx->h = h;
    MTX_LOG("Allocated matrix.");

    return mtx;
}

int mtx_assign(matrix *mtx1, const matrix *mtx2) {
    if (!mtx1 || !mtx2 || !mtx1->data || !mtx2->data) {
        MTX_LOG_ERROR("Invalid matrix pointers in assignment");
        return 1;
    }
    
    if (mtx1->w != mtx2->w || mtx1->h != mtx2->h) {
        MTX_LOG_ERROR("Matrix size mismatch in assignment");
        return -1;
    }
    
    memcpy(mtx1->data, mtx2->data, mtx1->w * mtx1->h * sizeof(double));
    MTX_LOG("Matrix assignment completed");
    
    return 0; 
}

matrix* mtx_copy(const matrix *mtx) {
    if(!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid source matrix for copy");
        return NULL;
    }

    matrix *new_mtx = mtx_alloc(mtx->w, mtx->h);
    if(!new_mtx) {
        MTX_LOG_ERROR("Failed to allocate matrix during copy");
        return NULL;
    }

    mtx_assign(new_mtx, mtx);
    MTX_LOG("Matrix copy completed");

    return new_mtx;
}

void mtx_free(matrix *mtx) {
    if (!mtx) {
        MTX_LOG_ERROR("Attempt to free NULL matrix");
        return;
    }

    free(mtx->data);
    free(mtx);
    MTX_LOG("Freed matrix");
}

double* mtx_ptr(matrix* mtx, size_t i, size_t j) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid matrix access attempt");
        return NULL;
    }
    return mtx->data + mtx->w * i + j;
}

const double* mtx_cptr(const matrix* mtx, size_t i, size_t j) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid matrix access attempt (const)");
        return NULL;
    }
    return mtx->data + mtx->w * i + j;
}

void mtx_set_zero(matrix *mtx) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid matrix in set_zero");
        return;
    } 
    
    memset(mtx->data, 0, mtx->w * mtx->h * sizeof(double));
    MTX_LOG("Matrix set to zero");
}

void mtx_set_id(matrix *mtx) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid matrix in set_id");
        return;
    }
    
    mtx_set_zero(mtx);
    
    if (mtx->w != mtx->h) {
        MTX_LOG_ERROR("Attempt to set identity for non-square matrix");
        return;
    } 
    
    for (size_t i = 0; i < mtx->w; i++) {
        *mtx_ptr(mtx, i, i) = 1.0;  
    }

    MTX_LOG("Matrix set to identity");
}

matrix* mtx_alloc_zero(size_t w, size_t h) {
    matrix *mtx = mtx_alloc(w, h);
    if(!mtx) {
        MTX_LOG_ERROR("Failed to allocate zero matrix");
        return NULL;
    }
    
    mtx_set_zero(mtx);
    return mtx;
}

matrix* mtx_alloc_id(size_t w, size_t h) {
    matrix *mtx = mtx_alloc(w, h);
    if(!mtx) {
        MTX_LOG_ERROR("Failed to allocate identity matrix");
        return NULL;
    }
    
    mtx_set_id(mtx);
    return mtx;
}

int mtx_input(matrix *mtx) {
    if (!mtx || !mtx->data) {
        MTX_LOG_ERROR("Invalid matrix in input");
        return -1;
    } 

    printf("Enter the matrix %zux%zu:\n", mtx->h, mtx->w);
    for (size_t i = 0; i < mtx->h; i++) {
        printf("Row %zu: ", i + 1);
        for (size_t j = 0; j < mtx->w; j++) {
            if (scanf("%lf", mtx->data + i * mtx->w + j) != 1) {
                MTX_LOG_ERROR("Matrix input error");
                return -1;
            }
        }
    }
    MTX_LOG("Matrix input completed");

    return 0;
}

void mtx_print(const matrix *m, int precision) {
    if (!m || !m->data) {
        MTX_LOG_ERROR("Attempt to print invalid matrix");
        printf("Wrong matrix!\n");
        return;
    }

    printf("Matrix %zux%zu:\n", m->h, m->w);
    for (size_t i = 0; i < m->h; i++) {
        for (size_t j = 0; j < m->w; j++) {
            printf("%.*f ", precision, *(m->data + i * m->w + j));
        }
        printf("\n");
    }
    MTX_LOG("Matrix printed");
}

size_t mtx_get_width(const matrix *mtx){
    if(!mtx){
        MTX_LOG_ERROR("Matrix is Null. Width cannot be gotten");
        return 0;
    }
    return mtx->w;
}

size_t mtx_get_height(const matrix *mtx){
    if(!mtx){
        MTX_LOG_ERROR("Matrix is Null. Height cannot be gotten");
        return 0;
    }
    
    return mtx->h;
}