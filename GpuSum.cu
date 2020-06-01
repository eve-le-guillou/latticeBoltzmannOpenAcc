#include "GpuSum.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "ArrayUtils.h"
#include "Check.h"
#include <stdio.h>
#include "math.h"
#include <cmath>

__global__ void gpu_abs_sub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C,
		int size, bool *divergence) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int ind = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
			+ threadIdx.x;
	*divergence = false;
	if (ind < size) {
		if(A[ind]!=A[ind]||B[ind]!=B[ind]) {
			*divergence=true;
		}
		C[ind] = abs(A[ind] - B[ind]);

	}
}

__global__ void gpu_abs_relSub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C,
		int size, bool *divergence) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int ind = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
			+ threadIdx.x;
	*divergence = false;
	if (ind < size) {
		if(A[ind]!=A[ind]||B[ind]!=B[ind]) {
			*divergence=true;
		}
		C[ind] = abs(A[ind] - B[ind]) / A[ind];

	}
}

__global__ void gpu_sqsub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C,
		int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int ind = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
			+ threadIdx.x;

	if (ind < size) {
		C[ind] = (A[ind] - B[ind]) * (A[ind] - B[ind]);
	}
}

__global__ void gpu_sqsubi(int *A, int *B, FLOAT_TYPE *C, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int ind = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
			+ threadIdx.x;
	if (ind < size) {
		C[ind] = (FLOAT_TYPE) (A[ind] - B[ind])
				* (FLOAT_TYPE) (A[ind] - B[ind]);
	}
}

__global__ void gpu_sum(FLOAT_TYPE *A, FLOAT_TYPE *B, int size) {
	int tid = threadIdx.x;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int gid = bid * blockDim.x * 2 + tid;

	extern __shared__ float temp[];
	temp[tid] = (gid < size) ? A[gid] : 0.0;
	temp[tid] += (gid + blockDim.x < size) ? A[gid + blockDim.x] : 0.0;
	__syncthreads();
	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s) {
			temp[tid] = temp[tid] + temp[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0) {
		B[bid] = temp[0];
	}
}

__global__ void gpu_sum256(FLOAT_TYPE *A, FLOAT_TYPE *B, int size) {
	int tid = threadIdx.x;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int gid = bid * blockDim.x * 2 + tid;

	extern __shared__ float temp[];
	temp[tid] = (gid < size) ? A[gid] : 0.0;
	temp[tid] += (gid + blockDim.x < size) ? A[gid + blockDim.x] : 0.0;
	__syncthreads();
	if (tid < 128)
		temp[tid] += temp[tid + 128];
	__syncthreads();
	if (tid < 64)
		temp[tid] += temp[tid + 64];
	__syncthreads();

	if (tid < 32) {
		temp[tid] += temp[tid + 32];
		__syncthreads();
		temp[tid] += temp[tid + 16];
		__syncthreads();
		temp[tid] += temp[tid + 8];
		__syncthreads();
		temp[tid] += temp[tid + 4];
		__syncthreads();
		temp[tid] += temp[tid + 2];
		__syncthreads();
		temp[tid] += temp[tid + 1];
	}

	if (tid == 0) {
		B[bid] = temp[0];
	}
}

__global__ void gpu_sum256_int(int *A, int *B, int size) {
	int tid = threadIdx.x;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int gid = bid * blockDim.x * 2 + tid;

	extern __shared__ int temp1[];
	temp1[tid] = (gid < size) ? A[gid] : 0;
	temp1[tid] += (gid + blockDim.x < size) ? A[gid + blockDim.x] : 0;
	__syncthreads();
	if (tid < 128)
		temp1[tid] += temp1[tid + 128];
	__syncthreads();
	if (tid < 64)
		temp1[tid] += temp1[tid + 64];
	__syncthreads();

	if (tid < 32) {
		temp1[tid] += temp1[tid + 32];
		__syncthreads();
		temp1[tid] += temp1[tid + 16];
		__syncthreads();
		temp1[tid] += temp1[tid + 8];
		__syncthreads();
		temp1[tid] += temp1[tid + 4];
		__syncthreads();
		temp1[tid] += temp1[tid + 2];
		__syncthreads();
		temp1[tid] += temp1[tid + 1];
	}

	if (tid == 0) {
		B[bid] = temp1[0];
	}
}

__global__ void gpu_max256(FLOAT_TYPE *A, FLOAT_TYPE *B, int size) {
	int tid = threadIdx.x;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int gid = bid * blockDim.x * 2 + tid;

	extern __shared__ float temp[];
	temp[tid] = (gid < size) ? A[gid] : 0.0;

	temp[tid] += (gid + blockDim.x < size) ? A[gid + blockDim.x] : 0.0;
	__syncthreads();
	if (tid < 128)
		temp[tid] =
				(temp[tid] >= temp[tid + 128]) ? temp[tid] : temp[tid + 128];
	__syncthreads();
	if (tid < 64)
		temp[tid] = (temp[tid] >= temp[tid + 64]) ? temp[tid] : temp[tid + 64];
	__syncthreads();

	if (tid < 32) {
		temp[tid] = (temp[tid] >= temp[tid + 32]) ? temp[tid] : temp[tid + 32];
		__syncthreads();
		temp[tid] = (temp[tid] >= temp[tid + 16]) ? temp[tid] : temp[tid + 16];
		__syncthreads();
		temp[tid] = (temp[tid] >= temp[tid + 8]) ? temp[tid] : temp[tid + 8];
		__syncthreads();
		temp[tid] = (temp[tid] >= temp[tid + 4]) ? temp[tid] : temp[tid + 4];
		__syncthreads();
		temp[tid] = (temp[tid] >= temp[tid + 2]) ? temp[tid] : temp[tid + 2];
		__syncthreads();
		temp[tid] = (temp[tid] >= temp[tid + 1]) ? temp[tid] : temp[tid + 1];
	}

	if (tid == 0) {
		B[bid] = temp[0];
	}
}

__host__ FLOAT_TYPE gpu_sum_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size) {
	dim3 grid_dim;

	int remaining = size;
	int shared_size = 256 * sizeof(FLOAT_TYPE);
	int req_blocks = 0;

	while (remaining > 1) {
		req_blocks = (remaining - 1) / 256 / 2 + 1;
		grid_dim.x = static_cast<int>(ceil(sqrt(req_blocks)));
		grid_dim.y = (req_blocks - 1) / grid_dim.x + 1;
		gpu_sum256<<<grid_dim, 256, shared_size>>>(C, D, remaining);

		//swap
		FLOAT_TYPE *temp = C;
		C = D;
		D = temp;

		remaining = req_blocks;
	}

	FLOAT_TYPE result = 0.0;
	CHECK(cudaMemcpy(&result, C, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost));
	return result;
}

__host__ int gpu_sum_int_h(int *C, int *D, int size) {
	dim3 grid_dim;

	int remaining = size;
	int shared_size = 256 * sizeof(int);
	int req_blocks = 0;

	while (remaining > 1) {
		req_blocks = (remaining - 1) / 256 / 2 + 1;
		grid_dim.x = static_cast<int>(ceil(sqrt(req_blocks)));
		grid_dim.y = (req_blocks - 1) / grid_dim.x + 1;
		gpu_sum256_int<<<grid_dim, 256, shared_size>>>(C, D, remaining);

		//swap
		int *temp = C;
		C = D;
		D = temp;

		remaining = req_blocks;
	}

	int result = 0;
	CHECK(cudaMemcpy(&result, C, sizeof(int), cudaMemcpyDeviceToHost));
	return result;
}

__host__ FLOAT_TYPE gpu_max_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size) {
	dim3 grid_dim;

	int remaining = size;
	int shared_size = 256 * sizeof(FLOAT_TYPE);
	int req_blocks = 0;
	while (remaining > 1) {
		req_blocks = (remaining - 1) / 256 / 2 + 1;
		grid_dim.x = static_cast<int>(ceil(sqrt(req_blocks)));
		grid_dim.y = (req_blocks - 1) / grid_dim.x + 1;
		gpu_max256<<<grid_dim, 256, shared_size>>>(C, D, remaining);

		//swap
		FLOAT_TYPE *temp = C;
		C = D;
		D = temp;

		remaining = req_blocks;
	}

	FLOAT_TYPE result = 0.0;
	CHECK(cudaMemcpy(&result, C, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost));
	return result;
}

__global__ void gpu_cond_copy_mask2D(FLOAT_TYPE *A, FLOAT_TYPE *B, int *mask,
		int value, int size) {
	int ind = ind = blockIdx.x * blockDim.x * blockDim.y * blockDim.z
			+ threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x
			+ threadIdx.x;
	if (ind < size) {
		A[ind] = ((mask[ind] & BND_ID_ALL) == BOUND_ID(value)) ? B[ind] : 0.0;
	}
}

__global__ void gpu_cond_copy_mask3D(FLOAT_TYPE *A, FLOAT_TYPE *B,
		int* bcBoundId_d, int value, int size) {
	int ind = ind = blockIdx.x * blockDim.x * blockDim.y * blockDim.z
			+ threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x
			+ threadIdx.x;
	if (ind < size) {
		A[ind] = bcBoundId_d[ind] == value ? B[ind] : 0.0;
	}
}



