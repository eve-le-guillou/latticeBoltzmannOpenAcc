/**
 * Functions for sum on GPU
 * @file GpuSum.h
 * @author Adam Koleszar (adam.koleszar@gmail.com) - functions used for the 2D solver
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - functions used for the 3D solver, to calculate new types of residuals
 */
#ifndef CPU_SUM_H
#define CPU_SUM_H

#include "FloatType.h"
void cpu_abs_sub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C, int size,bool *divergence);
/**
 * @brief Compute square of the difference of two vectors: \f$(A-B)^2\f$
 *
 * @param[in] A,B input vector
 * @param[out] C result vector
 * @param[in] size vector size
 */
void cpu_sqsub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C, int size);

/**
 * @brief Compute square of the difference of two vectors: \f$(A-B)^2\f$
 *
 * @param[in] A,B input vector
 * @param[out] C result vector
 * @param[in] size vector size
 */
void cpu_sqsubi(int *A, int *B, FLOAT_TYPE *C, int size);

/**
 * @brief Sum of a vector
 *
 * @param[in] A input vector
 * @param[out] B sum of vector (in the first element)
 * @param[in] size vector size
 */
void cpu_sum(FLOAT_TYPE *A, FLOAT_TYPE *B, int size);

void cpu_sum_int(FLOAT_TYPE *A, FLOAT_TYPE *B, int size);
/**
 * @brief Conditional vector copy for BC bitmask
 *
 * @param[in] A input vector
 * @param[out] B output vector
 * @param[in] mask array of conditional values
 * @param[in] value value to compare to
 * @param[in] size vector size
 */
void cpu_cond_copy_mask2D(FLOAT_TYPE *A, FLOAT_TYPE *B, int *mask, int value, int size);
void cpu_cond_copy_mask3D(FLOAT_TYPE *A, FLOAT_TYPE *B, int *bcBoundId_d, int value, int size);

/**
 * @brief Host function for GPU vector sum (calls #cpu_sum256)
 *
 * @param C input vector
 * @param D temporary vector
 * @param size vector size
 * @return sum of the vector
 */
FLOAT_TYPE cpu_sum_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size);
FLOAT_TYPE cpu_max_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size);

void cpu_abs_relSub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C,
                               int size, bool *divergence);

int cpu_sum_int_h(int *C, int *D, int size);
#endif
