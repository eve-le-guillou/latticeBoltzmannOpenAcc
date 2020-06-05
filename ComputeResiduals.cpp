#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "ComputeResiduals.h"
#include "CpuFunctions.h"
#include "FilesReading.h"
#include "ShellFunctions.h"
#include "ArrayUtils.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "CpuSum.h"

FLOAT_TYPE computeResidual2D(FLOAT_TYPE *f_d, FLOAT_TYPE *fTemp_d, FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d, int m, int n) {
	cpu_sqsub(f_d, fTemp_d, temp9a_d, 9 * m * n);
	return sqrt(cpu_sum_h(temp9a_d, temp9b_d, 9 * m * n));
}

FLOAT_TYPE computeResidual3D(FLOAT_TYPE *f_d, FLOAT_TYPE *fTemp_d, FLOAT_TYPE *temp19a_d, FLOAT_TYPE *temp19b_d, int m, int n, int h) {
	cpu_sqsub(f_d, fTemp_d, temp19a_d, 19 * m * n * h);
	return sqrt(cpu_sum_h(temp19a_d, temp19b_d, 19 * m * n * h));
}

FLOAT_TYPE computeNewResidual3D(FLOAT_TYPE *fn, FLOAT_TYPE *fnprev, FLOAT_TYPE *f1, FLOAT_TYPE *temp19a_d, FLOAT_TYPE *temp19b_d,
        int m, int n, int h) {
	cpu_NewResidual(fn, fnprev, f1, temp19a_d, 19 * m * n * h);
	return cpu_max_h(temp19a_d, temp19b_d, 19 * m * n * h);
}


FLOAT_TYPE computeDragLift2D(int *bcMask_d, FLOAT_TYPE *dl_d, FLOAT_TYPE *tempA_d, FLOAT_TYPE *tempB_d, int m, int n, int boundaryId) {
	cpu_cond_copy_mask2D(tempA_d, dl_d, bcMask_d, boundaryId,
			m * n);
	return cpu_sum_h(tempA_d, tempB_d, m * n);
}

FLOAT_TYPE computeDragLift3D(int *bcBoundId_d, FLOAT_TYPE *dl_d, FLOAT_TYPE *tempA_d, FLOAT_TYPE *tempB_d, int m, int n, int h, int boundaryId) {
	cpu_cond_copy_mask3D(tempA_d, dl_d, bcBoundId_d,
			boundaryId, m * n);
	return cpu_sum_h(tempA_d, tempB_d, m * n);
}

void cpu_NewResidual(FLOAT_TYPE *fn, FLOAT_TYPE *fnprev, FLOAT_TYPE *f1, FLOAT_TYPE *res, int size) {
	for (int ind=0; ind < size; ind++) {
		res[ind] = abs(abs(fn[ind] - fnprev[ind]) / f1[ind]);
	}
}
