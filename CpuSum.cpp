#include "CpuSum.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "ArrayUtils.h"
#include <stdio.h>
#include "math.h"
#include <cmath>

void cpu_abs_sub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C,
                            int size, bool *divergence) {
    *divergence = false;
    for (int ind=0; ind < size; ind++){
        if(A[ind]!=A[ind]||B[ind]!=B[ind]) {
            *divergence=true;
        }
        C[ind] = abs(A[ind] - B[ind]);
    }
}

void cpu_abs_relSub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C, int size, bool *divergence) {

    *divergence = false;
    for (int ind=0; ind < size; ind++){
        if(A[ind]!=A[ind]||B[ind]!=B[ind]) {
            *divergence=true;
        }
        C[ind] = abs(A[ind] - B[ind]) / A[ind];

    }
}

void cpu_sqsub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C, int size) {
    for (int ind=0; ind < size; ind++){
        C[ind] = (A[ind] - B[ind]) * (A[ind] - B[ind]);
    }
}

void cpu_sqsubi(int *A, int *B, FLOAT_TYPE *C, int size) {
    for (int ind=0; ind < size; ind++){
        C[ind] = (FLOAT_TYPE) (A[ind] - B[ind]) * (FLOAT_TYPE) (A[ind] - B[ind]);
    }
}

void cpu_sum(FLOAT_TYPE *A, FLOAT_TYPE *B, int size) {
    B[0] = 0;
    for (int ind=0; ind < size; ind++){
        B[0] += A[ind]
    }
}

void cpu_sum_int(int *A, int *B, int size) {
    B[0] = 0;
    for (int ind=0; ind < size; ind++){
        B[0] += A[ind]
    }
}

void cpu_max256(FLOAT_TYPE *A, FLOAT_TYPE *B, int size) {
    FLOAT_TYPE max = 0;
    for (int ind = 0; ind <size; ind ++){
        if (max < A[ind]) max = A[ind];
    }
    B[0] = max;
}

FLOAT_TYPE cpu_sum_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size) {
    cpu_sum(C, D, size);
    return D[0];
}

int cpu_sum_int_h(int *C, int *D, int size) {
    cpu_sum(C, D, size);
    return D[0];
}

FLOAT_TYPE cpu_max_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size) {
    cpu_max256(C, D, size);
    return D[0];
}

void cpu_cond_copy_mask2D(FLOAT_TYPE *A, FLOAT_TYPE *B, int *mask,
                                     int value, int size) {
    for (int ind=0; ind < size; ind++) {
        A[ind] = ((mask[ind] & BND_ID_ALL) == BOUND_ID(value)) ? B[ind] : 0.0;
    }
}

void cpu_cond_copy_mask3D(FLOAT_TYPE *A, FLOAT_TYPE *B,
                                     int* bcBoundId_d, int value, int size) {
    for (int ind=0; ind < size; ind++) {
        A[ind] = bcBoundId_d[ind] == value ? B[ind] : 0.0;
    }
}



