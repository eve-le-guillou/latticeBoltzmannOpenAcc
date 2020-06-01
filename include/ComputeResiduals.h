/**
 * Functions for residual computations and result comparing
 * @file ComputeResiduals.h
 * @author Adam Koleszar (adam.koleszar@gmail.com) - 2D version
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - development of the existing functions into 3D version
 */

#ifndef ComputeResiduals_H
#define ComputeResiduals_H

#include "FloatType.h"




/**
 * @brief Compute residual on GPU
 * @details compute the norm of the difference of f and fColl
 *  
 * @param f_d               distribution function
 * @param fColl_d           distribution function collision step
 * @param temp9a_d,temp9b_d temporary array for summation (size: 9xMxN)
 * @param m                 number of columns
 * @param n                 number of rows
 * @return norm of the difference
 */
__host__ FLOAT_TYPE computeResidual2D(FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d,
                                    FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d,
                                    int m, int n);

__host__ FLOAT_TYPE computeResidual3D(FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d,
                                    FLOAT_TYPE *temp19a_d, FLOAT_TYPE *temp19b_d,
                                    int m, int n, int h);
__host__ FLOAT_TYPE computeNewResidual3D(FLOAT_TYPE *fn, FLOAT_TYPE *fnprev,
                                    FLOAT_TYPE *f1, FLOAT_TYPE *temp19a_d, FLOAT_TYPE *temp19b_d,
                                    int m, int n, int h);
__global__ void gpu_NewResidual(FLOAT_TYPE *fn, FLOAT_TYPE *fnprev, FLOAT_TYPE *f1, FLOAT_TYPE *res, int size);
/**
 * @brief Compute drag or lift on GPU
 * 
 * @param bcMask_d         BC bitmask
 * @param dl_d             drag or lift
 * @param tempA_d, tempB_d temporary array for summation (size: MxN)
 * @param m                number of columns
 * @param n                number of rows
 * @param boundaryId       boundary to calculate drag/lift on
 * @return drag/lift
 */
__host__ FLOAT_TYPE computeDragLift2D(int *bcMask_d, FLOAT_TYPE *dl_d,
                                    FLOAT_TYPE *tempA_d, FLOAT_TYPE *tempB_d,
                                    int m, int n, int boundaryId);
									
__host__ FLOAT_TYPE computeDragLift3D(int *bcBoundId_d, FLOAT_TYPE *dl_d,
                                    FLOAT_TYPE *tempA_d, FLOAT_TYPE *tempB_d,
                                    int m, int n, int h, int boundaryId);
									

#endif
