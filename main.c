#include <stdio.h>

#include "mtx_repmem.h"
#include "mtx_actions.h"
#include "mtx_arithmetic.h"
#include "mtx_calcs.h"

int main(){
    matrix *mtx = mtx_alloc_id(3,3);
    *mtx_ptr(mtx,1,1) = 2;
    *mtx_ptr(mtx,2,2) = -1;
    matrix *res = mtx_exp(mtx, 1e-18);
    mtx_print(res, 5);
    mtx_free(mtx);
    mtx_free(res);

    // { 2x_1 + 4x_2 + 6x_3  = 8
    // |  x_1 + 8x_2 + 12x_3 = 4
    // { 4x_1 + 6x_2 + 10x_3 = 2
    // x_1 =  4
    // x_2 = 21
    // x_3 = -14

    matrix *A = mtx_alloc(3,3);
    matrix *B = mtx_alloc(1,3);

    mtx_input(A);
    mtx_input(B);

    matrix *X = mtx_solve_gauss(A,B);
    mtx_print(X,2);

    mtx_free(A);
    mtx_free(B);
    mtx_free(X);
    return 0;
}