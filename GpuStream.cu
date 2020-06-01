#include <stdio.h>

#include "GpuFunctions.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "GpuConstants.h"

__global__ void gpuStreaming2D(int* fluid_d, int* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d*length_d;
	FLOAT_TYPE *f, *mf;
	int n = length_d;
	if (ind < ms && fluid_d[ind] == 1)
	{
		f_d[ind] = fColl_d[ind];	//Update fNewStep = fColl
		f = f_d + ms;				// f is f_d memory positions but f starts in f_d 1st level==1st lattice direction
		mf = fColl_d + ms;
		f[ind]      = (stream_d[ind]      == 1) ? mf[ind-1]        : mf[ind];		// stream_d == 1 means that
		f[ind+ms]   = (stream_d[ind+ms]   == 1) ? mf[ind+ms-n]     : mf[ind+ms]; 	// the streaming is allowed
		f[ind+2*ms] = (stream_d[ind+2*ms] == 1) ? mf[ind+2*ms+1]   : mf[ind+2*ms];  // "the regular case"
		f[ind+3*ms] = (stream_d[ind+3*ms] == 1) ? mf[ind+3*ms+n]   : mf[ind+3*ms];  // stream_d != 1 means
		f[ind+4*ms] = (stream_d[ind+4*ms] == 1) ? mf[ind+4*ms-n-1] : mf[ind+4*ms]; 	// wall or node outside dom.
		f[ind+5*ms] = (stream_d[ind+5*ms] == 1) ? mf[ind+5*ms-n+1] : mf[ind+5*ms];
		f[ind+6*ms] = (stream_d[ind+6*ms] == 1) ? mf[ind+6*ms+n+1] : mf[ind+6*ms];
		f[ind+7*ms] = (stream_d[ind+7*ms] == 1) ? mf[ind+7*ms+n-1] : mf[ind+7*ms];
	}
}

__global__ void gpuStreaming2DCG(int* fluid_d, int* stream_d, FLOAT_TYPE* r_f_d, FLOAT_TYPE* r_fColl_d, FLOAT_TYPE* b_f_d, FLOAT_TYPE* b_fColl_d, int *cg_dir_d)
{
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d*length_d;
	FLOAT_TYPE *r_f, *r_mf, *b_f, *b_mf;
	int n = length_d;
	if (ind < ms)
	{
		int ori = cg_dir_d[ind];
		r_f_d[ind] = r_fColl_d[ind];	//Update fNewStep = fColl
		r_f = r_f_d + ms;				// f is r_f_d memory positions but f starts in r_f_d 1st level==1st lattice direction
		r_mf = r_fColl_d + ms;
		b_f_d[ind] = b_fColl_d[ind];	//Update fNewStep = fColl
		b_f = b_f_d + ms;				// f is r_f_d memory positions but f starts in r_f_d 1st level==1st lattice direction
		b_mf = b_fColl_d + ms;

		switch(ori){
		case 0:
			r_f[ind]      = r_mf[ind-1];
			r_f[ind+ms]   = r_mf[ind+ms-n];
			r_f[ind+2*ms] = r_mf[ind+2*ms+1];
			r_f[ind+3*ms] = r_mf[ind+3*ms+n];
			r_f[ind+4*ms] = r_mf[ind+4*ms-n-1];
			r_f[ind+5*ms] = r_mf[ind+5*ms-n+1];
			r_f[ind+6*ms] = r_mf[ind+6*ms+n+1];
			r_f[ind+7*ms] = r_mf[ind+7*ms+n-1];

			b_f[ind]      = b_mf[ind-1];
			b_f[ind+ms]   = b_mf[ind+ms-n];
			b_f[ind+2*ms] = b_mf[ind+2*ms+1];
			b_f[ind+3*ms] = b_mf[ind+3*ms+n];
			b_f[ind+4*ms] = b_mf[ind+4*ms-n-1];
			b_f[ind+5*ms] = b_mf[ind+5*ms-n+1];
			b_f[ind+6*ms] = b_mf[ind+6*ms+n+1];
			b_f[ind+7*ms] = b_mf[ind+7*ms+n-1];
			break;
		case 1: //NORTH
			r_f[ind]      = r_mf[ind-1];
			r_f[ind+ms]   = r_mf[ind+ms-n];
			r_f[ind+2*ms] = r_mf[ind+2*ms+1];
			r_f[ind+4*ms] = r_mf[ind+4*ms-n-1];
			r_f[ind+5*ms] = r_mf[ind+5*ms-n+1];

			b_f[ind]      = b_mf[ind-1];
			b_f[ind+ms]   = b_mf[ind+ms-n];
			b_f[ind+2*ms] = b_mf[ind+2*ms+1];
			b_f[ind+4*ms] = b_mf[ind+4*ms-n-1];
			b_f[ind+5*ms] = b_mf[ind+5*ms-n+1];
			break;
		case 2: //SOUTH
			r_f[ind]      = r_mf[ind-1];
			r_f[ind+2*ms] = r_mf[ind+2*ms+1];
			r_f[ind+3*ms] = r_mf[ind+3*ms+n];
			r_f[ind+6*ms] = r_mf[ind+6*ms+n+1];
			r_f[ind+7*ms] = r_mf[ind+7*ms+n-1];

			b_f[ind]      = b_mf[ind-1];
			b_f[ind+2*ms] = b_mf[ind+2*ms+1];
			b_f[ind+3*ms] = b_mf[ind+3*ms+n];
			b_f[ind+6*ms] = b_mf[ind+6*ms+n+1];
			b_f[ind+7*ms] = b_mf[ind+7*ms+n-1];
			break;
		case 3: //EAST
			r_f[ind]      = r_mf[ind-1];
			r_f[ind+ms]   = r_mf[ind+ms-n];
			r_f[ind+3*ms] = r_mf[ind+3*ms+n];
			r_f[ind+4*ms] = r_mf[ind+4*ms-n-1];
			r_f[ind+7*ms] = r_mf[ind+7*ms+n-1];

			b_f[ind]      = b_mf[ind-1];
			b_f[ind+ms]   = b_mf[ind+ms-n];
			b_f[ind+3*ms] = b_mf[ind+3*ms+n];
			b_f[ind+4*ms] = b_mf[ind+4*ms-n-1];
			b_f[ind+7*ms] = b_mf[ind+7*ms+n-1];
			break;
		case 4: //WEST
			r_f[ind+ms]   = r_mf[ind+ms-n];
			r_f[ind+2*ms] = r_mf[ind+2*ms+1];
			r_f[ind+3*ms] = r_mf[ind+3*ms+n];
			r_f[ind+5*ms] = r_mf[ind+5*ms-n+1];
			r_f[ind+6*ms] = r_mf[ind+6*ms+n+1];

			b_f[ind+ms]   = b_mf[ind+ms-n];
			b_f[ind+2*ms] = b_mf[ind+2*ms+1];
			b_f[ind+3*ms] = b_mf[ind+3*ms+n];
			b_f[ind+5*ms] = b_mf[ind+5*ms-n+1];
			b_f[ind+6*ms] = b_mf[ind+6*ms+n+1];
			break;
		case 5: //NE
			r_f[ind]      = r_mf[ind-1];
			r_f[ind+ms]   = r_mf[ind+ms-n];
			r_f[ind+4*ms] = r_mf[ind+4*ms-n-1];

			b_f[ind]      = b_mf[ind-1];
			b_f[ind+ms]   = b_mf[ind+ms-n];
			b_f[ind+4*ms] = b_mf[ind+4*ms-n-1];
			break;
		case 6: //NW
			r_f[ind+ms]   = r_mf[ind+ms-n];
			r_f[ind+2*ms] = r_mf[ind+2*ms+1];
			r_f[ind+5*ms] = r_mf[ind+5*ms-n+1];

			b_f[ind+ms]   = b_mf[ind+ms-n];
			b_f[ind+2*ms] = b_mf[ind+2*ms+1];
			b_f[ind+5*ms] = b_mf[ind+5*ms-n+1];
			break;
		case 7: //SE
			r_f[ind]      = r_mf[ind-1];
			r_f[ind+3*ms] = r_mf[ind+3*ms+n];
			r_f[ind+7*ms] = r_mf[ind+7*ms+n-1];

			b_f[ind]      = b_mf[ind-1];
			b_f[ind+3*ms] = b_mf[ind+3*ms+n];
			b_f[ind+7*ms] = b_mf[ind+7*ms+n-1];
			break;
		case 8: //SW
			r_f[ind+2*ms] = r_mf[ind+2*ms+1];
			r_f[ind+3*ms] = r_mf[ind+3*ms+n];
			r_f[ind+6*ms] = r_mf[ind+6*ms+n+1];

			b_f[ind+2*ms] = b_mf[ind+2*ms+1];
			b_f[ind+3*ms] = b_mf[ind+3*ms+n];
			b_f[ind+6*ms] = b_mf[ind+6*ms+n+1];
			break;
		default:
			break;
		}
	}
}

__global__ void gpuStreaming3D(int* fluid_d, bool* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
	int blockId = blockIdx.x
			+ blockIdx.y * gridDim.x;
	int ind =  blockId * (blockDim.x * blockDim.y)
																+ (threadIdx.y * blockDim.x)
																+ threadIdx.x;

	int ms = depth_d*length_d*height_d;
	FLOAT_TYPE *f, *mf;
	if (ind < ms && fluid_d[ind] == 1)
	{
		f_d[ind] = fColl_d[ind];	//Update fNewStep = fColl
		f = f_d + ms;				// f is f_d memory position but f starts in f_d 1st level==1st lattice direction
		mf = fColl_d + ms;
		f[ind+0  *ms]	=	(stream_d[ind+0	 *ms]	==	1)	?	mf[ind+0  *ms +	c3D_d[1	]]:	f[ind+0  *ms];
		f[ind+1	 *ms]	=	(stream_d[ind+1	 *ms]	==	1)	?	mf[ind+1  *ms +	c3D_d[2	]]: f[ind+1  *ms];
		f[ind+2	 *ms]	=	(stream_d[ind+2	 *ms]	==	1)	?	mf[ind+2  *ms +	c3D_d[3	]]:	f[ind+2  *ms];
		f[ind+3	 *ms]	=	(stream_d[ind+3	 *ms]	==	1)	?	mf[ind+3  *ms +	c3D_d[4	]]:	f[ind+3  *ms];
		f[ind+4	 *ms]	=	(stream_d[ind+4	 *ms]	==	1)	?	mf[ind+4  *ms +	c3D_d[5	]]:	f[ind+4  *ms];
		f[ind+5	 *ms]	=	(stream_d[ind+5	 *ms]	==	1)	?	mf[ind+5  *ms +	c3D_d[6	]]:	f[ind+5  *ms];
		f[ind+6	 *ms]	=	(stream_d[ind+6	 *ms]	==	1)	?	mf[ind+6  *ms +	c3D_d[7	]]:	f[ind+6  *ms];
		f[ind+7	 *ms]	=	(stream_d[ind+7	 *ms]	==	1)	?	mf[ind+7  *ms +	c3D_d[8	]]:	f[ind+7  *ms];
		f[ind+8	 *ms]	=	(stream_d[ind+8	 *ms]	==	1)	?	mf[ind+8  *ms +	c3D_d[9	]]:	f[ind+8  *ms];
		f[ind+9	 *ms]	=	(stream_d[ind+9	 *ms]	==	1)	?	mf[ind+9  *ms +	c3D_d[10]]:	f[ind+9  *ms];
		f[ind+10 *ms]	=	(stream_d[ind+10 *ms]	==	1)	?	mf[ind+10 *ms +	c3D_d[11]]:	f[ind+10 *ms];
		f[ind+11 *ms]	=	(stream_d[ind+11 *ms]	==	1)	?	mf[ind+11 *ms +	c3D_d[12]]:	f[ind+11 *ms];
		f[ind+12 *ms]	=	(stream_d[ind+12 *ms]	==	1)	?	mf[ind+12 *ms +	c3D_d[13]]:	f[ind+12 *ms];
		f[ind+13 *ms]	=	(stream_d[ind+13 *ms]	==	1)	?	mf[ind+13 *ms +	c3D_d[14]]:	f[ind+13 *ms];
		f[ind+14 *ms]	=	(stream_d[ind+14 *ms]	==	1)	?	mf[ind+14 *ms +	c3D_d[15]]:	f[ind+14 *ms];
		f[ind+15 *ms]	=	(stream_d[ind+15 *ms]	==	1)	?	mf[ind+15 *ms +	c3D_d[16]]:	f[ind+15 *ms];
		f[ind+16 *ms]	=	(stream_d[ind+16 *ms]	==	1)	?	mf[ind+16 *ms +	c3D_d[17]]:	f[ind+16 *ms];
		f[ind+17 *ms]	=	(stream_d[ind+17 *ms]	==	1)	?	mf[ind+17 *ms +	c3D_d[18]]:	f[ind+17 *ms];
	}
}
