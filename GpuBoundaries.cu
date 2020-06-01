/**
 * @author2D Adam Koleszar (adam.koleszar@gmail.com)
 * @author3D Maciej Kubat (m.j.kubat@cranfield.ac.uk) upgraded 2016 : In charge of BCmask recognition
 * @author3D Alfonso Aguilar Pontes (a.aguilar-pontes@cranfield.ac.uk) upgraded 2016 :  In charge of Boundary Conditions Equations and models
 *
 */
#include "GpuFunctions.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "GpuConstants.h"

__global__ void gpuBcInlet2D(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE* f_d,
		FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int size) {
	int bci = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d * length_d;

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC_FLUID) {
			if (bcMask_d[bci] & BC_INLT_B) {
				if ((bcMask_d[bci] & BC_INLT_N)
						&& !(bcMask_d[bci] & BC_CORNER)) {
					///@todo code: simplify code even though it will change accuracy
					f_d[ind + 4 * ms] = f_d[ind + 2 * ms];

					f_d[ind + 8 * ms] = f_d[ind + 6 * ms]
					                        + (f_d[ind] + f_d[ind + 1 * ms] + f_d[ind + 3 * ms]
					                                                              + 2
					                                                              * (f_d[ind + 2 * ms]
					                                                                     + f_d[ind + 6 * ms]
					                                                                           + f_d[ind + 5 * ms]))
					                                                                           * u0_d[ind] / 6.0;

					f_d[ind + 7 * ms] = f_d[ind + 5 * ms]
					                        - (f_d[ind] + f_d[ind + 1 * ms] + f_d[ind + 3 * ms]
					                                                              + 2
					                                                              * (f_d[ind + 2 * ms]
					                                                                     + f_d[ind + 6 * ms]
					                                                                           + f_d[ind + 5 * ms]))
					                                                                           * u0_d[ind] / 6.0;
					//printf("BCF%d 4:%f, 7:%f, 8:%f\n", ind, f_d[ind+4*ms], f_d[ind+7*ms], f_d[ind+8*ms]);
				}
				if ((bcMask_d[bci] & BC_INLT_W)
						&& !(bcMask_d[bci] & BC_CORNER)) {
					f_d[ind + 1 * ms] = f_d[ind + 3 * ms]
					                        + 2. * ((f_d[ind] + f_d[ind + 2 * ms] + f_d[ind + 4 * ms] + 2. * (f_d[ind + 3 * ms] + f_d[ind + 6 * ms]+ f_d[ind  + 7 * ms]))
					                        		/ (1.0 - u0_d[ind]))
					                        		* u0_d[ind] / 3;

					f_d[ind + 5 * ms] = f_d[ind + 7 * ms]
					                        + ((f_d[ind] + f_d[ind + 2 * ms] + f_d[ind + 4 * ms]
					                                                               + 2.
					                                                               * (f_d[ind + 3 * ms]
					                                                                      + f_d[ind + 6 * ms]
					                                                                            + f_d[ind + 7 * ms]))
					                        		/ (1.0 - u0_d[ind])) * u0_d[ind] / 6;

					f_d[ind + 8 * ms] = f_d[ind + 6 * ms]
					                        + ((f_d[ind] + f_d[ind + 2 * ms] + f_d[ind + 4 * ms]
					                                                               + 2.
					                                                               * (f_d[ind + 3 * ms]
					                                                                      + f_d[ind + 6 * ms]
					                                                                            + f_d[ind + 7 * ms]))
					                        		/ (1.0 - u0_d[ind])) * u0_d[ind] / 6;
				}
				///@todo code: compute inlet on other sides
			}
		}

	}
}

__global__ void gpuBcInlet3D(int *bcIdx_d, unsigned long long *bcMask_d,
		FLOAT_TYPE* f_d, FLOAT_TYPE* u1_d, FLOAT_TYPE* v1_d, FLOAT_TYPE* w1_d, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int bci = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																									+ threadIdx.x;

	int ms = depth_d * length_d * height_d;
	FLOAT_TYPE uW, vW, wW, dW;
	FLOAT_TYPE uE, vE, wE, dE;
	FLOAT_TYPE uS, vS, wS, dS;
	FLOAT_TYPE uN, vN, wN, dN;
	FLOAT_TYPE uB, vB, wB, dB;
	FLOAT_TYPE uT, vT, wT, dT;
	FLOAT_TYPE Nxy, Nxz, Nyx, Nyz, Nzx, Nzy;
	if (bci < size) {
		int ind = bcIdx_d[bci];
		//    printf("bcMask[ind] |= BC3D_MASK((unsigned long long)bcType[bci], dir); %#016lX\n", (bcMask_d[bci] & BC3D_FLUID) );
		if (bcMask_d[bci] & BC3D_FLUID) {
			if (bcMask_d[bci] & BC3D_INLT_B) {
				/*WEST*/if ((bcMask_d[bci] & BC3D_INLT_2)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 2))
								== (bcMask_d[bci] & BC3D_INLT_2))
								&& !(bcMask_d[bci] & BC3D_CORNER)) {
					//					printf("i am the west inlet %d %f \n",ind, -w1_d[ind]);//Check comment

					uW = u1_d[ind];
					vW = v1_d[ind];
					wW = w1_d[ind];
					dW = 1.0 / (1.0 + uW)
																																													* (f_d[ind + 0 * ms] + f_d[ind + 3 * ms]
																																													                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																													                                                     + f_d[ind + 6 * ms] + f_d[ind + 15 * ms]
																																													                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																													                                                                                                          + f_d[ind + 18 * ms]
																																													                                                                                                                + 2.0
																																													                                                                                                                * (f_d[ind + 2 * ms]
																																													                                                                                                                       + f_d[ind + 8 * ms]
																																													                                                                                                                             + f_d[ind + 10 * ms]
																																													                                                                                                                                   + f_d[ind + 12 * ms]
																																													                                                                                                                                         + f_d[ind + 14 * ms]));

					Nxy = 0.5
							* (f_d[ind + 3 * ms] + f_d[ind + 15 * ms]
							                           + f_d[ind + 17 * ms]
							                                 - (f_d[ind + 4 * ms] + f_d[ind + 16 * ms]
							                                                            + f_d[ind + 18 * ms]))
							                                                            - 1.0 / 3.0 * dW * vW;

					Nxz = 0.5
							* (f_d[ind + 5 * ms] + f_d[ind + 15 * ms]
							                           + f_d[ind + 16 * ms]
							                                 - (f_d[ind + 6 * ms] + f_d[ind + 17 * ms]
							                                                            + f_d[ind + 18 * ms]))
							                                                            - 1.0 / 3.0 * dW * wW;

					f_d[ind + 1 * ms] = f_d[ind + 2 * ms] + 1.0 / 3.0 * dW * uW;
					f_d[ind + 7 * ms] = f_d[ind + 10 * ms]
					                        + 1.0 / 6.0 * dW * (uW + vW) - Nxy;
					f_d[ind + 9 * ms] = f_d[ind + 8 * ms]
					                        + 1.0 / 6.0 * dW * (uW - vW) + Nxy;
					f_d[ind + 11 * ms] = f_d[ind + 14 * ms]
					                         + 1.0 / 6.0 * dW * (uW + wW) - Nxz;
					f_d[ind + 13 * ms] = f_d[ind + 12 * ms]
					                         + 1.0 / 6.0 * dW * (uW - wW) + Nxz;
				}
				/*EAST*/if ((bcMask_d[bci] & BC3D_INLT_1)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 1))
								== BC3D_INLT_1)
								&& !(bcMask_d[bci] & BC3D_CORNER)) {
					//					printf("i am the east inlet %d %f \n",ind, w1_d[ind]);//Check comment
					uE = u1_d[ind];
					vE = v1_d[ind];
					wE = w1_d[ind];
					dE = 1.0 / (1.0 - uE)
																																													* (f_d[ind + 0 * ms] + f_d[ind + 3 * ms]
																																													                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																													                                                     + f_d[ind + 6 * ms] + f_d[ind + 15 * ms]
																																													                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																													                                                                                                          + f_d[ind + 18 * ms]
																																													                                                                                                                + 2.0
																																													                                                                                                                * (f_d[ind + 1 * ms]
																																													                                                                                                                       + f_d[ind + 7 * ms]
																																													                                                                                                                             + f_d[ind + 9 * ms]
																																													                                                                                                                                   + f_d[ind + 11 * ms]
																																													                                                                                                                                         + f_d[ind + 13 * ms]));

					Nxy = 0.5
							* (f_d[ind + 3 * ms] + f_d[ind + 15 * ms]
							                           + f_d[ind + 17 * ms]
							                                 - (f_d[ind + 4 * ms] + f_d[ind + 16 * ms]
							                                                            + f_d[ind + 18 * ms]))
							                                                            - 1.0 / 3.0 * dE * vE;

					Nxz = 0.5
							* (f_d[ind + 5 * ms] + f_d[ind + 15 * ms]
							                           + f_d[ind + 16 * ms]
							                                 - (f_d[ind + 6 * ms] + f_d[ind + 17 * ms]
							                                                            + f_d[ind + 18 * ms]))
							                                                            - 1.0 / 3.0 * dE * wE;

					f_d[ind + 2 * ms] = f_d[ind + 1 * ms] - 1.0 / 3.0 * dE * uE;
					f_d[ind + 8 * ms] = f_d[ind + 9 * ms]
					                        - 1.0 / 6.0 * dE * (uE - vE) - Nxy;
					f_d[ind + 10 * ms] = f_d[ind + 7 * ms]
					                         - 1.0 / 6.0 * dE * (uE + vE) + Nxy;
					f_d[ind + 12 * ms] = f_d[ind + 13 * ms]
					                         - 1.0 / 6.0 * dE * (uE - wE) - Nxz;
					f_d[ind + 14 * ms] = f_d[ind + 11 * ms]
					                         - 1.0 / 6.0 * dE * (uE + wE) + Nxz;
				}
				/*SOUTH*/if ((bcMask_d[bci] & BC3D_INLT_4)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 4))
								== (bcMask_d[bci] & BC3D_INLT_4))
								&& !(bcMask_d[bci] & BC3D_CORNER)) {
					//					printf("i am the south inlet %d %f \n",ind, u1_d[ind]);//Check comment
					uS = u1_d[ind];
					vS = v1_d[ind];
					wS = w1_d[ind];
					dS = 1.0 / (1.0 - vS)
																																													* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																													                           + f_d[ind + 2 * ms] + f_d[ind + 5 * ms]
																																													                                                     + f_d[ind + 6 * ms] + f_d[ind + 11 * ms]
																																													                                                                               + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																													                                                                                                          + f_d[ind + 14 * ms]
																																													                                                                                                                + 2.0
																																													                                                                                                                * (f_d[ind + 4 * ms]
																																													                                                                                                                       + f_d[ind + 9 * ms]
																																													                                                                                                                             + f_d[ind + 10 * ms]
																																													                                                                                                                                   + f_d[ind + 16 * ms]
																																													                                                                                                                                         + f_d[ind + 18 * ms]));

					Nyx = 0.5
							* (f_d[ind + 1 * ms] + f_d[ind + 11 * ms]
							                           + f_d[ind + 13 * ms]
							                                 - (f_d[ind + 2 * ms] + f_d[ind + 12 * ms]
							                                                            + f_d[ind + 14 * ms]))
							                                                            - 1.0 / 3.0 * dS * uS;

					Nyz = 0.5
							* (f_d[ind + 5 * ms] + f_d[ind + 11 * ms]
							                           + f_d[ind + 12 * ms]
							                                 - (f_d[ind + 6 * ms] + f_d[ind + 13 * ms]
							                                                            + f_d[ind + 14 * ms]))
							                                                            - 1.0 / 3.0 * dS * wS;

					f_d[ind + 3 * ms] = f_d[ind + 4 * ms] + 1.0 / 3.0 * dS * vS;
					f_d[ind + 7 * ms] = f_d[ind + 10 * ms]
					                        + 1.0 / 6.0 * dS * (vS + uS) - Nyx;
					f_d[ind + 8 * ms] = f_d[ind + 9 * ms]
					                        + 1.0 / 6.0 * dS * (vS - uS) + Nyx;
					f_d[ind + 15 * ms] = f_d[ind + 18 * ms]
					                         + 1.0 / 6.0 * dS * (vS + wS) - Nyz;
					f_d[ind + 17 * ms] = f_d[ind + 16 * ms]
					                         + 1.0 / 6.0 * dS * (vS - wS) + Nyz;
				}
				/*NORTH*/if ((bcMask_d[bci] & BC3D_INLT_3)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 3))
								== (bcMask_d[bci] & BC3D_INLT_3))
								&& !(bcMask_d[bci] & BC3D_CORNER)) {
					//					printf("i am the north inlet %d %f \n",ind, u1_d[ind]);//Check comment
					uN = u1_d[ind];
					vN = v1_d[ind];
					wN = w1_d[ind];
					dN = 1.0 / (1.0 + vN)
																																													* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																													                           + f_d[ind + 2 * ms] + f_d[ind + 5 * ms]
																																													                                                     + f_d[ind + 6 * ms] + f_d[ind + 11 * ms]
																																													                                                                               + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																													                                                                                                          + f_d[ind + 14 * ms]
																																													                                                                                                                + 2.0
																																													                                                                                                                * (f_d[ind + 3 * ms]
																																													                                                                                                                       + f_d[ind + 7 * ms]
																																													                                                                                                                             + f_d[ind + 8 * ms]
																																													                                                                                                                                   + f_d[ind + 15 * ms]
																																													                                                                                                                                         + f_d[ind + 17 * ms]));

					Nyx = 0.5
							* (f_d[ind + 1 * ms] + f_d[ind + 11 * ms]
							                           + f_d[ind + 13 * ms]
							                                 - (f_d[ind + 2 * ms] + f_d[ind + 12 * ms]
							                                                            + f_d[ind + 14 * ms]))
							                                                            - 1.0 / 3.0 * dN * uN;

					Nyz = 0.5
							* (f_d[ind + 5 * ms] + f_d[ind + 11 * ms]
							                           + f_d[ind + 12 * ms]
							                                 - (f_d[ind + 6 * ms] + f_d[ind + 13 * ms]
							                                                            + f_d[ind + 14 * ms]))
							                                                            - 1.0 / 3.0 * dN * wN;

					f_d[ind + 4 * ms] = f_d[ind + 3 * ms] - 1.0 / 3.0 * dN * vN;
					f_d[ind + 9 * ms] = f_d[ind + 8 * ms]
					                        - 1.0 / 6.0 * dN * (vN - uN) - Nyx;
					f_d[ind + 10 * ms] = f_d[ind + 7 * ms]
					                         - 1.0 / 6.0 * dN * (vN + uN) + Nyx;
					f_d[ind + 16 * ms] = f_d[ind + 17 * ms]
					                         - 1.0 / 6.0 * dN * (vN - wN) - Nyz;
					f_d[ind + 18 * ms] = f_d[ind + 15 * ms]
					                         - 1.0 / 6.0 * dN * (vN + wN) + Nyz;
				}
				/*BOTTOM*/if ((bcMask_d[bci] & BC3D_INLT_6)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 6))
								== (bcMask_d[bci] & BC3D_INLT_6))
								&& !(bcMask_d[bci] & BC3D_CORNER)) {
					//					printf("i am the bottom inlet %d %f \n",ind, -u1_d[ind]);//Check comment

					uB = u1_d[ind];
					vB = v1_d[ind];
					wB = w1_d[ind];
					dB = 1.0 / (1.0 - wB)
																																													* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																													                           + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																													                                                     + f_d[ind + 4 * ms] + f_d[ind + 7 * ms]
																																													                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																													                                                                                                         + f_d[ind + 10 * ms]
																																													                                                                                                               + 2.0
																																													                                                                                                               * (f_d[ind + 6 * ms]
																																													                                                                                                                      + f_d[ind + 13 * ms]
																																													                                                                                                                            + f_d[ind + 14 * ms]
																																													                                                                                                                                  + f_d[ind + 17 * ms]
																																													                                                                                                                                        + f_d[ind + 18 * ms]));

					Nzx = 0.5
							* (f_d[ind + 1 * ms] + f_d[ind + 7 * ms]
							                           + f_d[ind + 9 * ms]
							                                 - (f_d[ind + 2 * ms] + f_d[ind + 8 * ms]
							                                                            + f_d[ind + 10 * ms]))
							                                                            - 1.0 / 3.0 * dB * uB;

					Nzy = 0.5
							* (f_d[ind + 3 * ms] + f_d[ind + 7 * ms]
							                           + f_d[ind + 8 * ms]
							                                 - (f_d[ind + 4 * ms] + f_d[ind + 9 * ms]
							                                                            + f_d[ind + 10 * ms]))
							                                                            - 1.0 / 3.0 * dB * vB;

					f_d[ind + 5 * ms] = f_d[ind + 6 * ms] + 1.0 / 3.0 * dB * wB;
					f_d[ind + 11 * ms] = f_d[ind + 14 * ms]
					                         + 1.0 / 6.0 * dB * (wB + uB) - Nzx;
					f_d[ind + 12 * ms] = f_d[ind + 13 * ms]
					                         + 1.0 / 6.0 * dB * (wB - uB) + Nzx;
					f_d[ind + 15 * ms] = f_d[ind + 18 * ms]
					                         + 1.0 / 6.0 * dB * (wB + vB) - Nzy;
					f_d[ind + 16 * ms] = f_d[ind + 17 * ms]
					                         + 1.0 / 6.0 * dB * (wB - vB) + Nzy;

				}
				/*TOP*/if ((bcMask_d[bci] & BC3D_INLT_5)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 5))
								== (bcMask_d[bci] & BC3D_INLT_5))
								&& !(bcMask_d[bci] & BC3D_CORNER)) {

					//printf("i am the top inlet %d %f \n",ind, u1_d[ind]);//Check comment
					/// Zou-He MODEL ///
					uT = u1_d[ind];
					vT = v1_d[ind];
					wT = w1_d[ind];
					dT = (1.0 / (1.0 + wT))
																																													* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																													                           + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																													                                                     + f_d[ind + 4 * ms] + f_d[ind + 7 * ms]
																																													                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																													                                                                                                         + f_d[ind + 10 * ms]
																																													                                                                                                               + 2.0
																																													                                                                                                               * (f_d[ind + 5 * ms]
																																													                                                                                                                      + f_d[ind + 11 * ms]
																																													                                                                                                                            + f_d[ind + 12 * ms]
																																													                                                                                                                                  + f_d[ind + 15 * ms]
																																													                                                                                                                                        + f_d[ind + 16 * ms]));

					Nzx = 0.5
							* (f_d[ind + 1 * ms] + f_d[ind + 7 * ms]
							                           + f_d[ind + 9 * ms]
							                                 - (f_d[ind + 2 * ms] + f_d[ind + 8 * ms]
							                                                            + f_d[ind + 10 * ms]))
							                                                            - 1.0 / 3.0 * dT * uT;

					Nzy = 0.5
							* (f_d[ind + 3 * ms] + f_d[ind + 7 * ms]
							                           + f_d[ind + 8 * ms]
							                                 - (f_d[ind + 4 * ms] + f_d[ind + 9 * ms]
							                                                            + f_d[ind + 10 * ms]))
							                                                            - 1.0 / 3.0 * dT * vT;

					f_d[ind + 6 * ms] = f_d[ind + 5 * ms] - 1.0 / 3.0 * dT * wT;
					f_d[ind + 13 * ms] = f_d[ind + 12 * ms]
					                         - 1.0 / 6.0 * dT * (wT - uT) - Nzx;
					f_d[ind + 14 * ms] = f_d[ind + 11 * ms]
					                         - 1.0 / 6.0 * dT * (wT + uT) + Nzx;
					f_d[ind + 17 * ms] = f_d[ind + 16 * ms]
					                         - 1.0 / 6.0 * dT * (wT - vT) - Nzy;
					f_d[ind + 18 * ms] = f_d[ind + 15 * ms]
					                         - 1.0 / 6.0 * dT * (wT + vT) + Nzy;
					//		            if(ind==26739) printf("i am here TOP lid u= %.14f f1: %.14f f2: %.14f f7: %.14f f8: %.14f "
					//		                        		"f9: %.14f f10: %.14f f11: %.14f f12: %.14f f13: %.14f f14: %.14f  \n",
					//		                        		u1_d[ind], f_d[ind+ 1*ms],f_d[ind+ 2*ms],f_d[ind+ 7*ms], f_d[ind+ 8*ms], f_d[ind+ 9*ms],
					//		                        		f_d[ind+10*ms],f_d[ind+11*ms],f_d[ind+12*ms],f_d[ind+13*ms],f_d[ind+14*ms]);
					//		            //				printf("f6: %f  f13:  %f  f14: %f  f17: %f  f18: %f  Nzx: %f  Nzy: %f\n",f_d[ind +  6*ms],f_d[ind +  13*ms],f_d[ind +  14*ms],f_d[ind +  17*ms],f_d[ind +  18*ms], Nzx, Nzy);

					//	/// He-Lou MODEL ///
					//				uT = u1_d[ind];
					//				vT = v1_d[ind];
					//				wT = w1_d[ind];
					////printf("ind: %d, u: %f, v: %f, w: %f\n",ind, uT, vT, wT);
					//				Nzx = 0.5*( f_d[ind +  1*ms] + f_d[ind +  7*ms] + f_d[ind +  9*ms] -
					//						( f_d[ind +  2*ms] + f_d[ind +  8*ms] + f_d[ind + 10*ms] ) ) - 1.0/3.0 * uT;
					//
					//				Nzy = 0.5*( f_d[ind +  3*ms] + f_d[ind +  7*ms] + f_d[ind +  8*ms] -
					//						( f_d[ind +  4*ms] + f_d[ind +  9*ms] + f_d[ind + 10*ms] ) ) - 1.0/3.0 * vT;
					//
					//				f_d[ind +  6*ms] = f_d[ind +  5*ms] - 1.0/3.0*dT*wT;
					//				f_d[ind + 13*ms] = f_d[ind + 12*ms] - 1.0/6.0*dT*(wT - uT) - Nzx;
					//				f_d[ind + 14*ms] = f_d[ind + 11*ms] - 1.0/6.0*dT*(wT + uT) + Nzx;
					//				f_d[ind + 17*ms] = f_d[ind + 16*ms] - 1.0/6.0*dT*(wT - vT) - Nzy;
					//				f_d[ind + 18*ms] = f_d[ind + 15*ms] - 1.0/6.0*dT*(wT + vT) + Nzy;
					//				printf("f6: %f  f13:  %f  f14: %f  f17: %f  f18: %f  Nzx: %f  Nzy: %f\n",f_d[ind +  6],f_d[ind +  6],f_d[ind +  6],f_d[ind +  6],f_d[ind +  6], Nzx, Nzy);
				}
			}
		}
	}
}
__global__ void gpuBcWall2D(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE *f_d,
		FLOAT_TYPE *fColl_d, FLOAT_TYPE *q_d, int size) {
	int bci = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d * length_d;
	int dir;

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC_FLUID) {
			for (dir = 1; dir < 9; ++dir) {
				if ((bcMask_d[bci] & BC_MASK(BC_ALL, dir))
						== (bcMask_d[bci] & BC_MASK(BC_WALL, dir))
						&& (bcMask_d[bci] & BC_MASK(BC_WALL, dir))) {
					// printf("%d: %X-%X-%X\n", ind, bcMask_d[bci], bcMask_d[bci] & BC_MASK(BC_ALL, dir), bcMask_d[bci] & BC_MASK(BC_WALL, dir));
					switch (boundaryType_d) {
					case CURVED: //curved
						if (q_d[bci * 8 + dir - 1] < 0.5) // if the distance from the boundary is less than 0.5?
						{
							f_d[ind + opp2D_d[dir]] = 2 * q_d[bci * 8 + dir]
							                                  * fColl_d[ind + dir * ms]
							                                            + (1 - 2 * q_d[bci * 8 + dir - 1])
							                                            * fColl_d[ind + dir * ms
							                                                      + c2D_d[dir]];

							// printf("if %d: dir: %d, Q:%f, F:%f; ", ind, dir, q_d[bci*8+dir], f_d[ind+opp_d[dir]]);
						} else {
							f_d[ind + opp2D_d[dir]] = fColl_d[ind + dir * ms]
							                                  / 2 / q_d[bci * 8 + dir - 1]
							                                            + (2 * q_d[bci * 8 + dir - 1] - 1)
							                                            / (2 * q_d[bci * 8 + dir - 1])
							                                            * fColl_d[ind + opp2D_d[dir]];

							// printf("else %d: dir: %d, Q:%f,  F: %f; ", ind, dir, q_d[bci*8+dir], f_d[ind+opp_d[dir]]);
						}
						break;
					case STRAIGHT: //half-way bounce back
						f_d[ind + opp2D_d[dir]] = f_d[ind + dir * ms];
						break;
					}
				}
			}
		}
	}
}

__global__ void gpuBcSimpleWall3D(int *bcIdx_d, unsigned long long *bcMask_d,
		FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d, FLOAT_TYPE *q_d, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int bci = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																									+ threadIdx.x;
	int ms = depth_d * length_d * height_d;
	int dir;

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC3D_FLUID) {
			for (dir = 1; dir < 19; ++dir) {
				if ((bcMask_d[bci] & BC3D_MASK(BC3D_ALL, dir))
						== (bcMask_d[bci] & BC3D_MASK(BC3D_WALL, dir))
						&& (bcMask_d[bci] & BC3D_MASK(BC3D_WALL, dir))) {
					// printf("%d: %X-%X-%X\n", ind, bcMask_d[bci], bcMask_d[bci] & BC3D_MASK(BC3D_ALL, dir), bcMask_d[bci] & BC3D_MASK(BC3D_WALL, dir));
					switch (boundaryType_d) { //TODO review curved 3D
					case CURVED: //curved
						if (q_d[bci * 18 + dir - 1] < 0.5) // if the distance from the boundary is less than 0.5?
						{
							f_d[ind + opp3D_d[dir]] = 2 * q_d[bci * 18 + dir]
							                                  * fColl_d[ind + dir * ms]
							                                            + (1 - 2 * q_d[bci * 18 + dir - 1])
							                                            * fColl_d[ind + dir * ms
							                                                      + c3D_d[dir]];

							// printf("if %d: dir: %d, Q:%f, F:%f; ", ind, dir, q_d[bci*8+dir], f_d[ind+opp_d[dir]]);
						} else {
							f_d[ind + opp3D_d[dir]] = fColl_d[ind + dir * ms]
							                                  / 2 / q_d[bci * 18 + dir - 1]
							                                            + (2 * q_d[bci * 18 + dir - 1] - 1)
							                                            / (2 * q_d[bci * 18 + dir - 1])
							                                            * fColl_d[ind + opp3D_d[dir]];

							// printf("else %d: dir: %d, Q:%f,  F: %f; ", ind, dir, q_d[bci*8+dir], f_d[ind+opp_d[dir]]);
						}
						break;
					case STRAIGHT: //half-way bounce back
						f_d[ind + opp3D_d[dir]] = f_d[ind + dir * ms]; //TODO imposed two times, the directions are paired
						break;
					}
				}
			}
		}
	}
}

__global__ void gpuBcComplexWall3D(int *bcIdx_d, unsigned long long *bcMask_d,
		FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d, FLOAT_TYPE *q_d, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int bci = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																									+ threadIdx.x;
	int ms = depth_d * length_d * height_d;
	int dir;

	int NoOfDirs;
	int x[3] = { 0, 0, 0 }; //normal vector
	int y[3] = { 0, 0, 0 }; //normal vector
	int z[3] = { 0, 0, 0 }; //normal vector
	FLOAT_TYPE uW, vW, wW, dW;
	FLOAT_TYPE uE, vE, wE, dE;
	FLOAT_TYPE uS, vS, wS, dS;
	FLOAT_TYPE uN, vN, wN, dN;
	FLOAT_TYPE uB, vB, wB, dB;
	FLOAT_TYPE uT, vT, wT, dT;
	FLOAT_TYPE Nxy, Nxz, Nyx, Nyz, Nzx, Nzy;
	FLOAT_TYPE Nx, Ny, Nz, C18, C27, C36, C45;
	if (bci < size) {
		int ind = bcIdx_d[bci];

		NoOfDirs = 0;
		for (dir = 1; dir < 19; ++dir) {
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_ALL, dir))
					== (bcMask_d[bci]
					             & BC3D_MASK((unsigned long long)BC3D_WALL, dir))
					             && (bcMask_d[bci]
					                          & BC3D_MASK((unsigned long long)BC3D_WALL, dir))) {
				NoOfDirs++;
			}
		}
		for (int j = 0; j < 3; j++) {
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_WALL, 1))
					&& x[0] == 0 && x[1] == 0) {
				x[j] = -1;
				continue;
			}
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_WALL, 2))
					&& x[0] == 0 && x[1] == 0) {
				x[j] = 1;
				continue;
			}
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_WALL, 3))
					&& y[0] == 0 && y[1] == 0) {
				y[j] = -1;
				continue;
			}
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_WALL, 4))
					&& y[0] == 0 && y[1] == 0) {
				y[j] = 1;
				continue;
			}
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_WALL, 5))
					&& z[0] == 0 && z[1] == 0) {
				z[j] = -1;
				continue;
			}
			if ((bcMask_d[bci] & BC3D_MASK((unsigned long long)BC3D_WALL, 6))
					&& z[0] == 0 && z[1] == 0) {
				z[j] = 1;
				continue;
			}
		}

		if (NoOfDirs == 5) { //MIDDLE OF WALL, TODO
			if (bcMask_d[bci] & BC3D_FLUID) {
				if (bcMask_d[bci] & BC3D_WALL_B) {
					/*WEST*/if ((bcMask_d[bci] & BC3D_WALL_2)
							&& ((bcMask_d[bci]
							              & BC3D_MASK((unsigned long long)BC3D_ALL, 2))
									== (bcMask_d[bci] & BC3D_WALL_2))) {
						//				        	printf("i am here W\n");

						uW = 0.0;
						vW = 0.0;
						wW = 0.0;
						dW = 1.0 / (1.0 + uW)
																																														* (f_d[ind + 0 * ms] + f_d[ind + 3 * ms]
																																														                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																														                                                     + f_d[ind + 6 * ms] + f_d[ind + 15 * ms]
																																														                                                                               + f_d[ind + 16 * ms]
																																														                                                                                     + f_d[ind + 17 * ms]
																																														                                                                                           + f_d[ind + 18 * ms]
																																														                                                                                                 + 2.0
																																														                                                                                                 * (f_d[ind + 2 * ms]
																																														                                                                                                        + f_d[ind + 8 * ms]
																																														                                                                                                              + f_d[ind + 10 * ms]
																																														                                                                                                                    + f_d[ind + 12 * ms]
																																														                                                                                                                          + f_d[ind + 14 * ms]));

						Nxy = 0.5
								* (f_d[ind + 3 * ms] + f_d[ind + 15 * ms]
								                           + f_d[ind + 17 * ms]
								                                 - (f_d[ind + 4 * ms]
								                                        + f_d[ind + 16 * ms]
								                                              + f_d[ind + 18 * ms]))
								                                              - 1.0 / 3.0 * dW * vW;

						Nxz = 0.5
								* (f_d[ind + 5 * ms] + f_d[ind + 15 * ms]
								                           + f_d[ind + 16 * ms]
								                                 - (f_d[ind + 6 * ms]
								                                        + f_d[ind + 17 * ms]
								                                              + f_d[ind + 18 * ms]))
								                                              - 1.0 / 3.0 * dW * wW;

						f_d[ind + 1 * ms] = f_d[ind + 2 * ms]
						                        + 1.0 / 3.0 * dW * uW;
						f_d[ind + 7 * ms] = f_d[ind + 10 * ms]
						                        + 1.0 / 6.0 * dW * (uW + vW) - Nxy;
						f_d[ind + 9 * ms] = f_d[ind + 8 * ms]
						                        + 1.0 / 6.0 * dW * (uW - vW) + Nxy;
						f_d[ind + 11 * ms] = f_d[ind + 14 * ms]
						                         + 1.0 / 6.0 * dW * (uW + wW) - Nxz;
						f_d[ind + 13 * ms] = f_d[ind + 12 * ms]
						                         + 1.0 / 6.0 * dW * (uW - wW) + Nxz;
					}
					/*EAST*/if ((bcMask_d[bci] & BC3D_WALL_1)
							&& ((bcMask_d[bci]
							              & BC3D_MASK((unsigned long long)BC3D_ALL, 1))
									== (bcMask_d[bci] & BC3D_WALL_1))) {
						//					printf("i am here E\n");

						uE = 0.0;
						vE = 0.0;
						wE = 0.0;
						dE = 1.0 / (1.0 - uE)
																																														* (f_d[ind + 0 * ms] + f_d[ind + 3 * ms]
																																														                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																														                                                     + f_d[ind + 6 * ms] + f_d[ind + 15 * ms]
																																														                                                                               + f_d[ind + 16 * ms]
																																														                                                                                     + f_d[ind + 17 * ms]
																																														                                                                                           + f_d[ind + 18 * ms]
																																														                                                                                                 + 2.0
																																														                                                                                                 * (f_d[ind + 1 * ms]
																																														                                                                                                        + f_d[ind + 7 * ms]
																																														                                                                                                              + f_d[ind + 9 * ms]
																																														                                                                                                                    + f_d[ind + 11 * ms]
																																														                                                                                                                          + f_d[ind + 13 * ms]));

						Nxy = 0.5
								* (f_d[ind + 3 * ms] + f_d[ind + 15 * ms]
								                           + f_d[ind + 17 * ms]
								                                 - (f_d[ind + 4 * ms]
								                                        + f_d[ind + 16 * ms]
								                                              + f_d[ind + 18 * ms]))
								                                              - 1.0 / 3.0 * dE * vE;

						Nxz = 0.5
								* (f_d[ind + 5 * ms] + f_d[ind + 15 * ms]
								                           + f_d[ind + 16 * ms]
								                                 - (f_d[ind + 6 * ms]
								                                        + f_d[ind + 17 * ms]
								                                              + f_d[ind + 18 * ms]))
								                                              - 1.0 / 3.0 * dE * wE;

						f_d[ind + 2 * ms] = f_d[ind + 1 * ms]
						                        - 1.0 / 3.0 * dE * uE;
						f_d[ind + 8 * ms] = f_d[ind + 9 * ms]
						                        - 1.0 / 6.0 * dE * (uE - vE) - Nxy;
						f_d[ind + 10 * ms] = f_d[ind + 7 * ms]
						                         - 1.0 / 6.0 * dE * (uE + vE) + Nxy;
						f_d[ind + 12 * ms] = f_d[ind + 13 * ms]
						                         - 1.0 / 6.0 * dE * (uE - wE) - Nxz;
						f_d[ind + 14 * ms] = f_d[ind + 11 * ms]
						                         - 1.0 / 6.0 * dE * (uE + wE) + Nxz;
					}
					/*SOUTH*/if ((bcMask_d[bci] & BC3D_WALL_4)
							&& ((bcMask_d[bci]
							              & BC3D_MASK((unsigned long long)BC3D_ALL, 4))
									== (bcMask_d[bci] & BC3D_WALL_4))) {
						//					printf("i am here S\n");

						uS = 0.0;
						vS = 0.0;
						wS = 0.0;
						dS = 1.0 / (1.0 - vS)
																																														* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																														                           + f_d[ind + 2 * ms] + f_d[ind + 5 * ms]
																																														                                                     + f_d[ind + 6 * ms] + f_d[ind + 11 * ms]
																																														                                                                               + f_d[ind + 12 * ms]
																																														                                                                                     + f_d[ind + 13 * ms]
																																														                                                                                           + f_d[ind + 14 * ms]
																																														                                                                                                 + 2.0
																																														                                                                                                 * (f_d[ind + 4 * ms]
																																														                                                                                                        + f_d[ind + 9 * ms]
																																														                                                                                                              + f_d[ind + 10 * ms]
																																														                                                                                                                    + f_d[ind + 16 * ms]
																																														                                                                                                                          + f_d[ind + 18 * ms]));

						Nyx = 0.5
								* (f_d[ind + 1 * ms] + f_d[ind + 11 * ms]
								                           + f_d[ind + 13 * ms]
								                                 - (f_d[ind + 2 * ms]
								                                        + f_d[ind + 12 * ms]
								                                              + f_d[ind + 14 * ms]))
								                                              - 1.0 / 3.0 * dS * uS;

						Nyz = 0.5
								* (f_d[ind + 5 * ms] + f_d[ind + 11 * ms]
								                           + f_d[ind + 12 * ms]
								                                 - (f_d[ind + 6 * ms]
								                                        + f_d[ind + 13 * ms]
								                                              + f_d[ind + 14 * ms]))
								                                              - 1.0 / 3.0 * dS * wS;

						f_d[ind + 3 * ms] = f_d[ind + 4 * ms]
						                        + 1.0 / 3.0 * dS * vS;
						f_d[ind + 7 * ms] = f_d[ind + 10 * ms]
						                        + 1.0 / 6.0 * dS * (vS + uS) - Nyx;
						f_d[ind + 8 * ms] = f_d[ind + 9 * ms]
						                        + 1.0 / 6.0 * dS * (vS - uS) + Nyx;
						f_d[ind + 15 * ms] = f_d[ind + 18 * ms]
						                         + 1.0 / 6.0 * dS * (vS + wS) - Nyz;
						f_d[ind + 17 * ms] = f_d[ind + 16 * ms]
						                         + 1.0 / 6.0 * dS * (vS - wS) + Nyz;
					}
					/*NORTH*/if ((bcMask_d[bci] & BC3D_WALL_3)
							&& ((bcMask_d[bci]
							              & BC3D_MASK((unsigned long long)BC3D_ALL, 3))
									== (bcMask_d[bci] & BC3D_WALL_3))) {
						//					printf("i am here N\n");

						uN = 0.0;
						vN = 0.0;
						wN = 0.0;
						dN = 1.0 / (1.0 + vN)
																																														* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																														                           + f_d[ind + 2 * ms] + f_d[ind + 5 * ms]
																																														                                                     + f_d[ind + 6 * ms] + f_d[ind + 11 * ms]
																																														                                                                               + f_d[ind + 12 * ms]
																																														                                                                                     + f_d[ind + 13 * ms]
																																														                                                                                           + f_d[ind + 14 * ms]
																																														                                                                                                 + 2.0
																																														                                                                                                 * (f_d[ind + 3 * ms]
																																														                                                                                                        + f_d[ind + 7 * ms]
																																														                                                                                                              + f_d[ind + 8 * ms]
																																														                                                                                                                    + f_d[ind + 15 * ms]
																																														                                                                                                                          + f_d[ind + 17 * ms]));

						Nyx = 0.5
								* (f_d[ind + 1 * ms] + f_d[ind + 11 * ms]
								                           + f_d[ind + 13 * ms]
								                                 - (f_d[ind + 2 * ms]
								                                        + f_d[ind + 12 * ms]
								                                              + f_d[ind + 14 * ms]))
								                                              - 1.0 / 3.0 * dN * uN;

						Nyz = 0.5
								* (f_d[ind + 5 * ms] + f_d[ind + 11 * ms]
								                           + f_d[ind + 12 * ms]
								                                 - (f_d[ind + 6 * ms]
								                                        + f_d[ind + 13 * ms]
								                                              + f_d[ind + 14 * ms]))
								                                              - 1.0 / 3.0 * dN * wN;

						f_d[ind + 4 * ms] = f_d[ind + 3 * ms]
						                        - 1.0 / 3.0 * dN * vN;
						f_d[ind + 9 * ms] = f_d[ind + 8 * ms]
						                        - 1.0 / 6.0 * dN * (vN - uN) - Nyx;
						f_d[ind + 10 * ms] = f_d[ind + 7 * ms]
						                         - 1.0 / 6.0 * dN * (vN + uN) + Nyx;
						f_d[ind + 16 * ms] = f_d[ind + 17 * ms]
						                         - 1.0 / 6.0 * dN * (vN - wN) - Nyz;
						f_d[ind + 18 * ms] = f_d[ind + 15 * ms]
						                         - 1.0 / 6.0 * dN * (vN + wN) + Nyz;
					}
					/*BOTTOM*/if ((bcMask_d[bci] & BC3D_WALL_6)
							&& ((bcMask_d[bci]
							              & BC3D_MASK((unsigned long long)BC3D_ALL, 6))
									== (bcMask_d[bci] & BC3D_WALL_6))) {
						//					printf("i am here B\n");

						//        	printf("i am bottom inlet\n");

						uB = 0.0;
						vB = 0.0;
						wB = 0.0;
						dB = 1.0 / (1.0 - wB)
																																														* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																														                           + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																														                                                     + f_d[ind + 4 * ms] + f_d[ind + 7 * ms]
																																														                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																														                                                                                                         + f_d[ind + 10 * ms]
																																														                                                                                                               + 2.0
																																														                                                                                                               * (f_d[ind + 6 * ms]
																																														                                                                                                                      + f_d[ind + 13 * ms]
																																														                                                                                                                            + f_d[ind + 14 * ms]
																																														                                                                                                                                  + f_d[ind + 17 * ms]
																																														                                                                                                                                        + f_d[ind + 18 * ms]));

						Nzx = 0.5
								* (f_d[ind + 1 * ms] + f_d[ind + 7 * ms]
								                           + f_d[ind + 9 * ms]
								                                 - (f_d[ind + 2 * ms] + f_d[ind + 8 * ms]
								                                                            + f_d[ind + 10 * ms]))
								                                                            - 1.0 / 3.0 * dB * uB;

						Nzy = 0.5
								* (f_d[ind + 3 * ms] + f_d[ind + 7 * ms]
								                           + f_d[ind + 8 * ms]
								                                 - (f_d[ind + 4 * ms] + f_d[ind + 9 * ms]
								                                                            + f_d[ind + 10 * ms]))
								                                                            - 1.0 / 3.0 * dB * vB;

						f_d[ind + 5 * ms] = f_d[ind + 6 * ms]
						                        + 1.0 / 3.0 * dB * wB;
						f_d[ind + 11 * ms] = f_d[ind + 14 * ms]
						                         + 1.0 / 6.0 * dB * (wB + uB) - Nzx;
						f_d[ind + 12 * ms] = f_d[ind + 13 * ms]
						                         + 1.0 / 6.0 * dB * (wB - uB) + Nzx;
						f_d[ind + 15 * ms] = f_d[ind + 18 * ms]
						                         + 1.0 / 6.0 * dB * (wB + vB) - Nzy;
						f_d[ind + 16 * ms] = f_d[ind + 17 * ms]
						                         + 1.0 / 6.0 * dB * (wB - vB) + Nzy;

					}
					/*TOP*/if ((bcMask_d[bci] & BC3D_WALL_5)
							&& ((bcMask_d[bci]
							              & BC3D_MASK((unsigned long long)BC3D_ALL, 5))
									== (bcMask_d[bci] & BC3D_WALL_5))) {
						//		        	printf("i am here T\n");

						/// Zou-He MODEL ///
						uT = 0.0;
						vT = 0.0;
						wT = 0.0;
						dT = (1.0 / (1.0 + wT))
																																														* (f_d[ind + 0 * ms] + f_d[ind + 1 * ms]
																																														                           + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																														                                                     + f_d[ind + 4 * ms] + f_d[ind + 7 * ms]
																																														                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																														                                                                                                         + f_d[ind + 10 * ms]
																																														                                                                                                               + 2.0
																																														                                                                                                               * (f_d[ind + 5 * ms]
																																														                                                                                                                      + f_d[ind + 11 * ms]
																																														                                                                                                                            + f_d[ind + 12 * ms]
																																														                                                                                                                                  + f_d[ind + 15 * ms]
																																														                                                                                                                                        + f_d[ind + 16 * ms]));

						Nzx = 0.5
								* (f_d[ind + 1 * ms] + f_d[ind + 7 * ms]
								                           + f_d[ind + 9 * ms]
								                                 - (f_d[ind + 2 * ms] + f_d[ind + 8 * ms]
								                                                            + f_d[ind + 10 * ms]))
								                                                            - 1.0 / 3.0 * dT * uT;

						Nzy = 0.5
								* (f_d[ind + 3 * ms] + f_d[ind + 7 * ms]
								                           + f_d[ind + 8 * ms]
								                                 - (f_d[ind + 4 * ms] + f_d[ind + 9 * ms]
								                                                            + f_d[ind + 10 * ms]))
								                                                            - 1.0 / 3.0 * dT * vT;

						f_d[ind + 6 * ms] = f_d[ind + 5 * ms]
						                        - 1.0 / 3.0 * dT * wT;
						f_d[ind + 13 * ms] = f_d[ind + 12 * ms]
						                         - 1.0 / 6.0 * dT * (wT - uT) - Nzx;
						f_d[ind + 14 * ms] = f_d[ind + 11 * ms]
						                         - 1.0 / 6.0 * dT * (wT + uT) + Nzx;
						f_d[ind + 17 * ms] = f_d[ind + 16 * ms]
						                         - 1.0 / 6.0 * dT * (wT - vT) - Nzy;
						f_d[ind + 18 * ms] = f_d[ind + 15 * ms]
						                         - 1.0 / 6.0 * dT * (wT + vT) + Nzy;
						//				printf("f6: %f  f13:  %f  f14: %f  f17: %f  f18: %f  Nzx: %f  Nzy: %f\n",f_d[ind +  6*ms],f_d[ind +  13*ms],f_d[ind +  14*ms],f_d[ind +  17*ms],f_d[ind +  18*ms], Nzx, Nzy);
					}
				}
			}
		}
		if (NoOfDirs == 9) { //EDGES, TODO
			if (x[0] == 1 && z[1] == 1) { //BOTTOM WEST
				//					if(ind==30) printf("i am here BOTTOM WEST_1 f3: %f f4: %f\n",f_d[ind+ 3*ms],f_d[ind+ 4*ms]);
				//					if(ind==840) printf("i am here BOTTOM WEST f3: %f f4: %f\n",f_d[ind+ 3*ms],f_d[ind+ 4*ms]);
				//		        	printf("i am here BOTTOM WEST\n");
				Ny = -0.25 * (f_d[ind + 3 * ms] - f_d[ind + 4 * ms]);

				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms] + Ny;
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms] - Ny;
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms] + Ny;
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms] - Ny;
				//Buried links
				f_d[ind + 12 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 13 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (x[0] == -1 && z[1] == 1) { //BOTTOM EAST
				//		        	printf("i am here BOTTOM EAST\n");

				Ny = -0.25 * (f_d[ind + 3 * ms] - f_d[ind + 4 * ms]);

				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms];
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms] + Ny;
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms] - Ny;
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms] + Ny;
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms] - Ny;

				//Buried links
				f_d[ind + 11 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                                                    + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 14 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                                                    + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                                                    + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (y[0] == 1 && z[1] == 1) { //BOTTOM SOUTH
				Nx = -0.25 * (f_d[ind + 1 * ms] - f_d[ind + 2 * ms]);

				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms] + Nx;
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms] - Nx;
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms] + Nx;
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms] - Nx;
				//Buried links
				f_d[ind + 16 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 17 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (y[0] == -1 && z[1] == 1) { //BOTTOM NORTH
				Nx = -0.25 * (f_d[ind + 1 * ms] - f_d[ind + 2 * ms]);

				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms];
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms] + Nx;
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms] - Nx;
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms] + Nx;
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms] - Nx;
				//Buried links
				f_d[ind + 15 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 17 * ms]);
				f_d[ind + 18 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 17 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 17 * ms]);

			}
			if (x[0] == 1 && z[1] == -1) { //TOP WEST
				//		        	printf("i am here TOP WEST\n");

				Ny = -0.25 * (f_d[ind + 3 * ms] - f_d[ind + 4 * ms]);

				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms] + Ny;
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms] - Ny;
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms] + Ny;
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms] - Ny;
				//Buried links
				f_d[ind + 11 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                                                    + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 14 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                                                    + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                                                    + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (x[0] == -1 && z[1] == -1) { //TOP EAST
				//		        	printf("i am here TOP EAST\n");
				Ny = -0.25 * (f_d[ind + 3 * ms] - f_d[ind + 4 * ms]);

				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms];
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms] + Ny;
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms] - Ny;
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms] + Ny;
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms] - Ny;
				//Buried links
				f_d[ind + 12 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 13 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (y[0] == 1 && z[1] == -1) { //TOP SOUTH
				Nx = -0.25 * (f_d[ind + 1 * ms] - f_d[ind + 2 * ms]);

				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms] + Nx;
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms] - Nx;
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms] + Nx;
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms] - Nx;
				//Buried links
				f_d[ind + 15 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 17 * ms]);
				f_d[ind + 18 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 17 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 17 * ms]);
			}
			if (y[0] == -1 && z[1] == -1) { //TOP NORTH
				Nx = -0.25 * (f_d[ind + 1 * ms] - f_d[ind + 2 * ms]);

				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms];
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms] + Nx;
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms] - Nx;
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms] + Nx;
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms] - Nx;
				//Buried links
				f_d[ind + 16 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 17 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 8 * ms] + f_d[ind + 9 * ms]
																																												                                                                                                                         + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                                                    + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                                               + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (x[0] == 1 && y[1] == -1) { //NORTH WEST
				Nz = -0.25 * (f_d[ind + 5 * ms] - f_d[ind + 6 * ms]);

				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms];
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms] + Nz;
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms] - Nz;
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms] + Nz;
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms] - Nz;
				//Buried links
				f_d[ind + 7 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 10 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}
			if (x[0] == -1 && y[1] == -1) { //NORTH EAST
				Nz = -0.25 * (f_d[ind + 5 * ms] - f_d[ind + 6 * ms]);

				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms];
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms] + Nz;
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms] - Nz;
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms] + Nz;
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms] - Nz;
				//Buried links
				f_d[ind + 8 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                     + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                           + f_d[ind + 18 * ms]);
				f_d[ind + 9 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                     + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                           + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                     + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                           + f_d[ind + 18 * ms]);

			}
			if (x[0] == 1 && y[1] == 1) { //SOUTH WEST
				Nz = -0.25 * (f_d[ind + 5 * ms] - f_d[ind + 6 * ms]);

				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms];
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms] + Nz;
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms] - Nz;
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms] + Nz;
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms] - Nz;
				//Buried links
				f_d[ind + 8 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                     + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                           + f_d[ind + 18 * ms]);
				f_d[ind + 9 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                     + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                           + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                     + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                                + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                           + f_d[ind + 18 * ms]);

			}
			if (x[0] == -1 && y[1] == 1) { //SOUTH EAST
				Nz = -0.25 * (f_d[ind + 5 * ms] - f_d[ind + 6 * ms]);

				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms];
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms] + Nz;
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms] - Nz;
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms] + Nz;
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms] - Nz;
				//Buried links
				f_d[ind + 7 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				f_d[ind + 10 * ms] = (1. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);
				//Resting link
				f_d[ind + 0 * ms] = (12. / 22.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 12 * ms] + f_d[ind + 13 * ms]
																																												                                                                                                                                                    + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                                               + f_d[ind + 16 * ms] + f_d[ind + 17 * ms]
																																												                                                                                                                                                                                                          + f_d[ind + 18 * ms]);

			}

		}
		if (NoOfDirs == 12) { //CORNERS, TODO
			if (x[0] == 1 && y[1] == 1 && z[2] == 1) { //SOUTH BOTTOM WEST C1
				//		        	printf("i am here C1\n");

				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms];
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms];
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms];
				C18 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                     + f_d[ind + 18 * ms]);
				f_d[ind + 8 * ms] = C18;
				f_d[ind + 9 * ms] = C18;
				f_d[ind + 12 * ms] = C18;
				f_d[ind + 13 * ms] = C18;
				f_d[ind + 16 * ms] = C18;
				f_d[ind + 17 * ms] = C18;
				f_d[ind + 0 * ms] = 12. * C18;
			}
			if (x[0] == -1 && y[1] == 1 && z[2] == 1) { //SOUTH BOTTOM EAST C2
				//		        	printf("i am here C2\n");

				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms];
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms];
				f_d[ind + 15 * ms] = f_d[ind + 18 * ms];
				C27 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                         + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                    + f_d[ind + 18 * ms]);
				f_d[ind + 7 * ms] = C27;
				f_d[ind + 10 * ms] = C27;
				f_d[ind + 11 * ms] = C27;
				f_d[ind + 14 * ms] = C27;
				f_d[ind + 16 * ms] = C27;
				f_d[ind + 17 * ms] = C27;
				f_d[ind + 0 * ms] = 12. * C27;

			}
			if (x[0] == 1 && y[1] == -1 && z[2] == 1) { //NORTH BOTTOM WEST C3
				//		        	printf("i am here C3\n");

				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms];
				f_d[ind + 11 * ms] = f_d[ind + 14 * ms];
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms];
				C36 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                    + f_d[ind + 17 * ms]);
				f_d[ind + 7 * ms] = C36;
				f_d[ind + 10 * ms] = C36;
				f_d[ind + 12 * ms] = C36;
				f_d[ind + 13 * ms] = C36;
				f_d[ind + 15 * ms] = C36;
				f_d[ind + 18 * ms] = C36;
				f_d[ind + 0 * ms] = 12. * C36;
			}
			if (x[0] == -1 && y[1] == -1 && z[2] == 1) { //SOUTH BOTTOM EAST C4
				//		        	printf("i am here C4\n");

				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 5 * ms] = f_d[ind + 6 * ms];
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms];
				f_d[ind + 12 * ms] = f_d[ind + 13 * ms];
				f_d[ind + 16 * ms] = f_d[ind + 17 * ms];
				C45 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                          + f_d[ind + 13 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                     + f_d[ind + 17 * ms]);
				f_d[ind + 8 * ms] = C45;
				f_d[ind + 9 * ms] = C45;
				f_d[ind + 11 * ms] = C45;
				f_d[ind + 14 * ms] = C45;
				f_d[ind + 15 * ms] = C45;
				f_d[ind + 18 * ms] = C45;
				f_d[ind + 0 * ms] = 12. * C45;
			}
			if (x[0] == 1 && y[1] == 1 && z[2] == -1) { //SOUTH TOP WEST C5
				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 7 * ms] = f_d[ind + 10 * ms];
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms];
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms];
				C45 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                          + f_d[ind + 13 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                     + f_d[ind + 17 * ms]);
				f_d[ind + 8 * ms] = C45;
				f_d[ind + 9 * ms] = C45;
				f_d[ind + 11 * ms] = C45;
				f_d[ind + 14 * ms] = C45;
				f_d[ind + 15 * ms] = C45;
				f_d[ind + 18 * ms] = C45;
				f_d[ind + 0 * ms] = 12. * C45;
			}
			if (x[0] == -1 && y[1] == 1 && z[2] == -1) { //SOUTH TOP EAST C6
				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 3 * ms] = f_d[ind + 4 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 8 * ms] = f_d[ind + 9 * ms];
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms];
				f_d[ind + 17 * ms] = f_d[ind + 16 * ms];
				C36 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                         + f_d[ind + 14 * ms] + f_d[ind + 16 * ms]
																																												                                                                                                                                                    + f_d[ind + 17 * ms]);
				f_d[ind + 7 * ms] = C36;
				f_d[ind + 10 * ms] = C36;
				f_d[ind + 12 * ms] = C36;
				f_d[ind + 13 * ms] = C36;
				f_d[ind + 15 * ms] = C36;
				f_d[ind + 18 * ms] = C36;
				f_d[ind + 0 * ms] = 12. * C36;

			}
			if (x[0] == 1 && y[1] == -1 && z[2] == -1) { //NORTH TOP WEST C7
				f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 9 * ms] = f_d[ind + 8 * ms];
				f_d[ind + 13 * ms] = f_d[ind + 12 * ms];
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms];
				C27 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 8 * ms]
																																												                                                                                               + f_d[ind + 9 * ms] + f_d[ind + 12 * ms]
																																												                                                                                                                         + f_d[ind + 13 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                    + f_d[ind + 18 * ms]);
				f_d[ind + 7 * ms] = C27;
				f_d[ind + 10 * ms] = C27;
				f_d[ind + 11 * ms] = C27;
				f_d[ind + 14 * ms] = C27;
				f_d[ind + 16 * ms] = C27;
				f_d[ind + 17 * ms] = C27;
				f_d[ind + 0 * ms] = 12. * C27;

			}
			if (x[0] == -1 && y[1] == -1 && z[2] == -1) { //NORTH TOP EAST C8
				f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
				f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
				f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
				f_d[ind + 10 * ms] = f_d[ind + 7 * ms];
				f_d[ind + 14 * ms] = f_d[ind + 11 * ms];
				f_d[ind + 18 * ms] = f_d[ind + 15 * ms];
				C18 = (1. / 18.)
																																												* (f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
																																												                                           + f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
																																												                                                                     + f_d[ind + 6 * ms] + f_d[ind + 7 * ms]
																																												                                                                                               + f_d[ind + 10 * ms] + f_d[ind + 11 * ms]
																																												                                                                                                                          + f_d[ind + 14 * ms] + f_d[ind + 15 * ms]
																																												                                                                                                                                                     + f_d[ind + 18 * ms]);
				f_d[ind + 8 * ms] = C18;
				f_d[ind + 9 * ms] = C18;
				f_d[ind + 12 * ms] = C18;
				f_d[ind + 13 * ms] = C18;
				f_d[ind + 16 * ms] = C18;
				f_d[ind + 17 * ms] = C18;
				f_d[ind + 0 * ms] = 12. * C18;

			}

		}

	}
	//printf("bci:%d NoOfDirs:%d\n",bci, NoOfDirs);
}

__global__ void gpuBcOutlet2D(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE *f_d,
		FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d, int size) {
	int bci = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d * length_d;

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC_FLUID) {
			switch (outletProfile_d) {
			case OUTLET: //Zou-He
				if ((bcMask_d[bci] & BC_OUTL_E) == BC_OUTL_E) {
					///@todo code: simplify code even though it will change accuracy
					f_d[ind + 3 * ms] =
							f_d[ind + 1 * ms]
							    - 2
							    * ((f_d[ind] + f_d[ind + 2 * ms]
							                       + f_d[ind + 4 * ms]
							                             + 2
							                             * (f_d[ind + 1 * ms]
							                                    + f_d[ind
							                                          + 5
							                                          * ms]
							                                          + f_d[ind
							                                                + 8
							                                                * ms]))
							    		/ (1 - u0_d[ind]))
							    		* u0_d[ind] / 3;

					f_d[ind + 7 * ms] = f_d[ind + 5 * ms]
					                        - ((f_d[ind] + f_d[ind + 2 * ms] + f_d[ind + 4 * ms]
					                                                               + 2
					                                                               * (f_d[ind + 1 * ms]
					                                                                      + f_d[ind + 5 * ms]
					                                                                            + f_d[ind + 8 * ms]))
					                        		/ (1 - u0_d[ind])) * u0_d[ind] / 6;

					f_d[ind + 6 * ms] = f_d[ind + 8 * ms]
					                        - ((f_d[ind] + f_d[ind + 2 * ms] + f_d[ind + 4 * ms]
					                                                               + 2
					                                                               * (f_d[ind + 1 * ms]
					                                                                      + f_d[ind + 5 * ms]
					                                                                            + f_d[ind + 8 * ms]))
					                        		/ (1 - u0_d[ind])) * u0_d[ind] / 6;
				}
				if (bcMask_d[bci] & BC_OUTL_N) {
					///@todo code: fill north-side outlet
				}
				if (bcMask_d[bci] & BC_OUTL_W) {
					///@todo code: fill west-side outlet
				}
				if (bcMask_d[bci] & BC_OUTL_S) {
					///@todo code: fill south-side outlet
				}
				break;
			case OUTLET_SECOND: //open boundary
				if ((bcMask_d[bci] & BC_OUTL_E) == BC_OUTL_E) {
					f_d[ind + ms] = 2 * f_d[ind + ms - 1] - f_d[ind + ms - 2];
					f_d[ind + 5 * ms] = 2 * f_d[ind + 5 * ms - 1]
					                            - f_d[ind + 5 * ms - 2];
					f_d[ind + 8 * ms] = 2 * f_d[ind + 8 * ms - 1]
					                            - f_d[ind + 8 * ms - 2];
				}
				if (bcMask_d[bci] & BC_OUTL_N) {
					///@todo code: fill north-side outlet
				}
				if (bcMask_d[bci] & BC_OUTL_W) {
					///@todo code: fill west-side outlet
				}
				if (bcMask_d[bci] & BC_OUTL_S) {
					///@todo code: fill south-side outlet
				}
				break;
			case OUTLET_FIRST: //first order
				if ((bcMask_d[bci] & BC_OUTL_E) == BC_OUTL_E) {
					f_d[ind + ms] = f_d[ind - 1 + ms];
					f_d[ind + 5 * ms] = f_d[ind - 1 + 5 * ms];
					f_d[ind + 8 * ms] = f_d[ind - 1 + 8 * ms];
				}
				if (bcMask_d[bci] & BC_OUTL_N) {
					///@todo code: fill north-side outlet
				}
				if (bcMask_d[bci] & BC_OUTL_W) {
					///@todo code: fill west-side outlet
				}
				if (bcMask_d[bci] & BC_OUTL_S) {
					///@todo code: fill south-side outlet
				}
				break;
			}
		}
	}
}

__global__ void gpuBcOutlet3D(int *bcIdx_d, unsigned long long *bcMask_d,
		FLOAT_TYPE *f_d, FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *w_d, FLOAT_TYPE *rho_d,
		int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int bci = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																									+ threadIdx.x;
	int ms = depth_d * length_d;
	int n = length_d;
	int l = length_d * depth_d;

	FLOAT_TYPE r, u, v, w;

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC3D_FLUID) {
			switch (outletProfile_d) {
			case OUTLET: //Zou-He
				if ((bcMask_d[bci] & BC3D_OUTL_1) == BC3D_OUTL_1) {
					printf(
							"Zou-He for Outlet not implemented, go for open boundary 2nd order\n");

					///@todo code: simplify code even though it will change accuracy
					//            f_d[ind+3*ms] = f_d[ind+1*ms]
					//                          - 2*( ( f_d[ind]
					//                                + f_d[ind+2*ms]
					//                                + f_d[ind+4*ms]
					//                                + 2*( f_d[ind+1*ms]
					//                                    + f_d[ind+5*ms]
					//                                    + f_d[ind+8*ms])
					//                                )/(1-u0_d[ind])
					//                              )*u0_d[ind]/3;
					//
					//            f_d[ind+7*ms] = f_d[ind+5*ms]
					//                          - ( ( f_d[ind]
					//                              + f_d[ind+2*ms]
					//                              + f_d[ind+4*ms]
					//                              + 2*( f_d[ind+1*ms]
					//                                  + f_d[ind+5*ms]
					//                                  + f_d[ind+8*ms])
					//                              )/(1-u0_d[ind])
					//                            )*u0_d[ind]/6;
					//
					//            f_d[ind+6*ms] = f_d[ind+8*ms]
					//                          - ( ( f_d[ind]
					//                              + f_d[ind+2*ms]
					//                              + f_d[ind+4*ms]
					//                              + 2*( f_d[ind+1*ms]
					//                                  + f_d[ind+5*ms]
					//                                  + f_d[ind+8*ms])
					//                              )/(1-u0_d[ind])
					//                            )*u0_d[ind]/6;
				}
				if ((bcMask_d[bci] & BC3D_OUTL_3) == BC3D_OUTL_3) {
					///@todo code: fill north-side outlet
				}
				if ((bcMask_d[bci] & BC3D_OUTL_2) == BC3D_OUTL_2) {
					///@todo code: fill west-side outlet
				}
				if ((bcMask_d[bci] & BC3D_OUTL_4) == BC3D_OUTL_4) {
					///@todo code: fill south-side outlet
				}
				if ((bcMask_d[bci] & BC3D_OUTL_5) == BC3D_OUTL_5) {
					///@todo code:
				}
				if ((bcMask_d[bci] & BC3D_OUTL_6) == BC3D_OUTL_6) {
					///@todo code:
				}
				break;
			case OUTLET_SECOND: //open boundary
				if ((bcMask_d[bci] & BC3D_OUTL_1) == BC3D_OUTL_1) {
					//				printf("i am here");
					f_d[ind + 2 * ms] = 2 * f_d[ind + 2 * ms - 1]
					                            - f_d[ind + 2 * ms - 2];
					f_d[ind + 8 * ms] = 2 * f_d[ind + 8 * ms - 1]
					                            - f_d[ind + 8 * ms - 2];
					f_d[ind + 10 * ms] = 2 * f_d[ind + 10 * ms - 1]
					                             - f_d[ind + 10 * ms - 2];
					f_d[ind + 12 * ms] = 2 * f_d[ind + 12 * ms - 1]
					                             - f_d[ind + 12 * ms - 2];
					f_d[ind + 14 * ms] = 2 * f_d[ind + 14 * ms - 1]
					                             - f_d[ind + 14 * ms - 2];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_3) == BC3D_OUTL_3) {
					f_d[ind + 4 * ms] = 2 * f_d[ind + 4 * ms + n]
					                            - f_d[ind + 4 * ms + 2 * n];
					f_d[ind + 9 * ms] = 2 * f_d[ind + 9 * ms + n]
					                            - f_d[ind + 9 * ms + 2 * n];
					f_d[ind + 10 * ms] = 2 * f_d[ind + 10 * ms + n]
					                             - f_d[ind + 10 * ms + 2 * n];
					f_d[ind + 16 * ms] = 2 * f_d[ind + 16 * ms + n]
					                             - f_d[ind + 16 * ms + 2 * n];
					f_d[ind + 18 * ms] = 2 * f_d[ind + 18 * ms + n]
					                             - f_d[ind + 18 * ms + 2 * n];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_2) == BC3D_OUTL_2) {
					f_d[ind + 1 * ms] = 2 * f_d[ind + 1 * ms + 1]
					                            - f_d[ind + 1 * ms + 2];
					f_d[ind + 7 * ms] = 2 * f_d[ind + 7 * ms + 1]
					                            - f_d[ind + 7 * ms + 2];
					f_d[ind + 9 * ms] = 2 * f_d[ind + 9 * ms + 1]
					                            - f_d[ind + 9 * ms + 2];
					f_d[ind + 11 * ms] = 2 * f_d[ind + 11 * ms + 1]
					                             - f_d[ind + 11 * ms + 2];
					f_d[ind + 13 * ms] = 2 * f_d[ind + 13 * ms + 1]
					                             - f_d[ind + 13 * ms + 2];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_4) == BC3D_OUTL_4) {
					f_d[ind + 3 * ms] = 2 * f_d[ind + 3 * ms + n]
					                            - f_d[ind + 3 * ms + 2 * n];
					f_d[ind + 7 * ms] = 2 * f_d[ind + 7 * ms + n]
					                            - f_d[ind + 7 * ms + 2 * n];
					f_d[ind + 8 * ms] = 2 * f_d[ind + 8 * ms + n]
					                            - f_d[ind + 8 * ms + 2 * n];
					f_d[ind + 15 * ms] = 2 * f_d[ind + 15 * ms + n]
					                             - f_d[ind + 15 * ms + 2 * n];
					f_d[ind + 17 * ms] = 2 * f_d[ind + 17 * ms + n]
					                             - f_d[ind + 17 * ms + 2 * n];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_5) == BC3D_OUTL_5) {
					f_d[ind + 6 * ms] = 2 * f_d[ind + 6 * ms - l]
					                            - f_d[ind + 6 * ms - 2 * l];
					f_d[ind + 13 * ms] = 2 * f_d[ind + 13 * ms - l]
					                             - f_d[ind + 13 * ms - 2 * l];
					f_d[ind + 14 * ms] = 2 * f_d[ind + 14 * ms - l]
					                             - f_d[ind + 14 * ms - 2 * l];
					f_d[ind + 17 * ms] = 2 * f_d[ind + 17 * ms - l]
					                             - f_d[ind + 17 * ms - 2 * l];
					f_d[ind + 18 * ms] = 2 * f_d[ind + 18 * ms - l]
					                             - f_d[ind + 18 * ms - 2 * l];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_6) == BC3D_OUTL_6) {
					f_d[ind + 5 * ms] = 2 * f_d[ind + 5 * ms + l]
					                            - f_d[ind + 5 * ms + 2 * l];
					f_d[ind + 11 * ms] = 2 * f_d[ind + 11 * ms + l]
					                             - f_d[ind + 11 * ms + 2 * l];
					f_d[ind + 12 * ms] = 2 * f_d[ind + 12 * ms + l]
					                             - f_d[ind + 12 * ms + 2 * l];
					f_d[ind + 15 * ms] = 2 * f_d[ind + 15 * ms + l]
					                             - f_d[ind + 15 * ms + 2 * l];
					f_d[ind + 16 * ms] = 2 * f_d[ind + 16 * ms + l]
					                             - f_d[ind + 16 * ms + 2 * l];
				}
				break;
			case OUTLET_FIRST: //first order
				if ((bcMask_d[bci] & BC3D_OUTL_1) == BC3D_OUTL_1) {
					f_d[ind + 2 * ms] = f_d[ind + 2 * ms - 1];
					f_d[ind + 8 * ms] = f_d[ind + 8 * ms - 1];
					f_d[ind + 10 * ms] = f_d[ind + 10 * ms - 1];
					f_d[ind + 12 * ms] = f_d[ind + 12 * ms - 1];
					f_d[ind + 14 * ms] = f_d[ind + 14 * ms - 1];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_3) == BC3D_OUTL_3) {
					f_d[ind + 4 * ms] = f_d[ind + 4 * ms - n];
					f_d[ind + 9 * ms] = f_d[ind + 9 * ms - n];
					f_d[ind + 10 * ms] = f_d[ind + 10 * ms - n];
					f_d[ind + 16 * ms] = f_d[ind + 16 * ms - n];
					f_d[ind + 18 * ms] = f_d[ind + 18 * ms - n];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_2) == BC3D_OUTL_2) {
					f_d[ind + 1 * ms] = f_d[ind + 1 * ms + 1];
					f_d[ind + 7 * ms] = f_d[ind + 7 * ms + 1];
					f_d[ind + 9 * ms] = f_d[ind + 9 * ms + 1];
					f_d[ind + 11 * ms] = f_d[ind + 11 * ms + 1];
					f_d[ind + 13 * ms] = f_d[ind + 13 * ms + 1];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_4) == BC3D_OUTL_4) {
					f_d[ind + 3 * ms] = f_d[ind + 3 * ms + n];
					f_d[ind + 7 * ms] = f_d[ind + 7 * ms + n];
					f_d[ind + 8 * ms] = f_d[ind + 8 * ms + n];
					f_d[ind + 15 * ms] = f_d[ind + 15 * ms + n];
					f_d[ind + 17 * ms] = f_d[ind + 17 * ms + n];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_5) == BC3D_OUTL_5) {
					f_d[ind + 6 * ms] = f_d[ind + 6 * ms - l];
					f_d[ind + 13 * ms] = f_d[ind + 13 * ms - l];
					f_d[ind + 14 * ms] = f_d[ind + 14 * ms - l];
					f_d[ind + 17 * ms] = f_d[ind + 17 * ms - l];
					f_d[ind + 18 * ms] = f_d[ind + 18 * ms - l];
				}
				if ((bcMask_d[bci] & BC3D_OUTL_6) == BC3D_OUTL_6) {
					f_d[ind + 5 * ms] = f_d[ind + 5 * ms + l];
					f_d[ind + 11 * ms] = f_d[ind + 11 * ms + l];
					f_d[ind + 12 * ms] = f_d[ind + 12 * ms + l];
					f_d[ind + 15 * ms] = f_d[ind + 15 * ms + l];
					f_d[ind + 16 * ms] = f_d[ind + 16 * ms + l];
				}
				break;
			case OUTLET_HEmodel: //
				u = u_d[ind];
				v = v_d[ind];
				w = w_d[ind];
				r = rho_d[ind];
				//					printf("%f %f %f %f\n", u, v,w,r);
				if ((bcMask_d[bci] & BC3D_OUTL_1) == BC3D_OUTL_1) {
					f_d[ind + 2 * ms] = feqc3D(u, cx3D_d[ 2 ], v, cy3D_d[ 2 ], w, cz3D_d[ 2 ], r, w3D_d[ 2 ]);
					f_d[ind + 8 * ms] = feqc3D(u, cx3D_d[ 8 ], v, cy3D_d[ 8 ], w, cz3D_d[ 8 ], r, w3D_d[ 8 ]);
					f_d[ind + 10 * ms] = feqc3D(u, cx3D_d[10 ], v, cy3D_d[10 ], w, cz3D_d[10 ], r, w3D_d[10 ]);
					f_d[ind + 12 * ms] = feqc3D(u, cx3D_d[12 ], v, cy3D_d[12 ], w, cz3D_d[12 ], r, w3D_d[12 ]);
					f_d[ind + 14 * ms] = feqc3D(u, cx3D_d[14 ], v, cy3D_d[14 ], w, cz3D_d[14 ], r, w3D_d[14 ]);
				}
				if ((bcMask_d[bci] & BC3D_OUTL_3) == BC3D_OUTL_3) {
					//TODO
				}
				if ((bcMask_d[bci] & BC3D_OUTL_2) == BC3D_OUTL_2) {
					//TODO
				}
				if ((bcMask_d[bci] & BC3D_OUTL_4) == BC3D_OUTL_4) {
					//TODO
				}
				if ((bcMask_d[bci] & BC3D_OUTL_5) == BC3D_OUTL_5) {
					//TODO
				}
				if ((bcMask_d[bci] & BC3D_OUTL_6) == BC3D_OUTL_6) {
					//TODO
				}
				break;
			}
		}
	}
}

__global__ void gpuBcPeriodic2D(int *bcIdx_d, int *bcMask_d,
		FLOAT_TYPE* r_f_d,FLOAT_TYPE* b_f_d, int size, int *orientation_d, int test_case, FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d, FLOAT_TYPE *rho_d,
		FLOAT_TYPE *u_d, FLOAT_TYPE *v_d) {
	int ind = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = depth_d * length_d;
	int offsetY = length_d * (depth_d - 1);
	int offsetX = length_d - 1;
	if (ind < ms) {

		int ori = orientation_d[ind];
		switch(ori){
		case 1: //NORTH
			if(test_case == 2){
				FLOAT_TYPE u_temp = u_d[ind], v_temp = v_d[ind];
				FLOAT_TYPE r_temp = (1. / (1. + v_temp)) * (r_f_d[ind] + r_f_d[ind + 1 * ms] + r_f_d[ind + 3 * ms]
				                                                                                     + 2 * (r_f_d[ind + 2 * ms] + r_f_d[ind + 6 * ms] + r_f_d[ind + 5 * ms]) );
				FLOAT_TYPE b_temp = (1. / (1. + v_temp)) * (b_f_d[ind] + b_f_d[ind + 1 * ms] + b_f_d[ind + 3 * ms]
				                                                                                     + 2 * (b_f_d[ind + 2 * ms] + b_f_d[ind + 6 * ms] + b_f_d[ind + 5 * ms]) );

				r_f_d[ind + 4 * ms] = r_f_d[ind + 2 * ms] - (2. / 3.) * r_temp*v_temp;
				r_f_d[ind + 7 * ms] = r_f_d[ind + 5 * ms] + 0.5 * (r_f_d[ind + 1 * ms] - r_f_d[ind + 3 * ms]) - (1. / 6.) * (r_temp*v_temp) - 0.5 * (r_temp*u_temp);
				r_f_d[ind + 8 * ms] = r_f_d[ind + 6 * ms] + 0.5 * (r_f_d[ind + 3 * ms] - r_f_d[ind + 1 * ms]) - (1. / 6.) * (r_temp*v_temp) + 0.5 * (r_temp*u_temp);

				b_f_d[ind + 4 * ms] = b_f_d[ind + 2 * ms] - (2. / 3.) * b_temp*v_temp;
				b_f_d[ind + 7 * ms] = b_f_d[ind + 5 * ms] + 0.5 * (b_f_d[ind + 1 * ms] - b_f_d[ind + 3 * ms]) - (1. / 6.) * (b_temp*v_temp) - 0.5 * (b_temp*u_temp);
				b_f_d[ind + 8 * ms] = b_f_d[ind + 6 * ms] + 0.5 * (b_f_d[ind + 3 * ms] - b_f_d[ind + 1 * ms]) - (1. / 6.) * (b_temp*v_temp) + 0.5 * (b_temp*u_temp);

				r_rho_d[ind] = r_temp;
				b_rho_d[ind] = b_temp;
				rho_d[ind] = r_temp + b_temp;
			}
			else if(test_case == 6){
				r_f_d[ind + 4 * ms] = r_f_d[ind + 2 * ms];
				r_f_d[ind + 7 * ms] = r_f_d[ind + 5 * ms];
				r_f_d[ind + 8 * ms] = r_f_d[ind + 6 * ms];

				b_f_d[ind + 4 * ms] = b_f_d[ind + 2 * ms];
				b_f_d[ind + 7 * ms] = b_f_d[ind + 5 * ms];
				b_f_d[ind + 8 * ms] = b_f_d[ind + 6 * ms];
			}
			else {
				r_f_d[ind + 4 * ms] = r_f_d[ind + 4 * ms - offsetY];
				r_f_d[ind + 7 * ms] = r_f_d[ind + 7 * ms - offsetY];
				r_f_d[ind + 8 * ms] = r_f_d[ind + 8 * ms - offsetY];

				b_f_d[ind + 4 * ms] = b_f_d[ind + 4 * ms - offsetY];
				b_f_d[ind + 7 * ms] = b_f_d[ind + 7 * ms - offsetY];
				b_f_d[ind + 8 * ms] = b_f_d[ind + 8 * ms - offsetY];
			}
			break;
		case 2: //SOUTH
			if(test_case == 2){
				FLOAT_TYPE u_temp = u_d[ind], v_temp = v_d[ind];

				FLOAT_TYPE r_temp = (1. / (1. - v_temp)) *
						(r_f_d[ind] + r_f_d[ind + 1 * ms] + r_f_d[ind + 3 * ms] + 2 * (r_f_d[ind + 4 * ms] + r_f_d[ind + 7 * ms] + r_f_d[ind + 8 * ms]));
				FLOAT_TYPE b_temp = (1. / (1. - v_temp)) *
						(b_f_d[ind] + b_f_d[ind + 1 * ms] + b_f_d[ind + 3 * ms] + 2 * (b_f_d[ind + 4 * ms] + b_f_d[ind + 7 * ms] + b_f_d[ind + 8 * ms]));

				r_f_d[ind + 2 * ms] = r_f_d[ind + 4 * ms] + (2. / 3.) * r_temp*v_temp;
				r_f_d[ind + 5 * ms] = r_f_d[ind + 7 * ms] - 0.5 * (r_f_d[ind + 1 * ms] - r_f_d[ind + 3 * ms]) + (1. / 6.) * (r_temp*v_temp) + 0.5 * (r_temp*u_temp);
				r_f_d[ind + 6 * ms] = r_f_d[ind + 8 * ms] + 0.5 * (r_f_d[ind + 1 * ms] - r_f_d[ind + 3 * ms]) + (1. / 6.) * (r_temp*v_temp) - 0.5 * (r_temp*u_temp);

				b_f_d[ind + 2 * ms] = b_f_d[ind + 4 * ms] + (2. / 3.) * b_temp*v_temp;
				b_f_d[ind + 5 * ms] = b_f_d[ind + 7 * ms] - 0.5 * (b_f_d[ind + 1 * ms] - b_f_d[ind + 3 * ms]) + (1. / 6.) * (b_temp*v_temp) + 0.5 * (b_temp*u_temp);
				b_f_d[ind + 6 * ms] = b_f_d[ind + 8 * ms] + 0.5 * (b_f_d[ind + 1 * ms] - b_f_d[ind + 3 * ms]) + (1. / 6.) * (b_temp*v_temp) - 0.5 * (b_temp*u_temp);

				r_rho_d[ind] = r_temp;
				b_rho_d[ind] = b_temp;
				rho_d[ind] = r_temp + b_temp;
			}
			else if(test_case == 6){
				r_f_d[ind + 2 * ms] = r_f_d[ind + 4 * ms];
				r_f_d[ind + 5 * ms] = r_f_d[ind + 7 * ms];
				r_f_d[ind + 6 * ms] = r_f_d[ind + 8 * ms];

				b_f_d[ind + 2 * ms] = b_f_d[ind + 4 * ms];
				b_f_d[ind + 5 * ms] = b_f_d[ind + 7 * ms];
				b_f_d[ind + 6 * ms] = b_f_d[ind + 8 * ms];
			}
			else{
				r_f_d[ind + 2 * ms] = r_f_d[ind + 2 * ms + offsetY];
				r_f_d[ind + 5 * ms] = r_f_d[ind + 5 * ms + offsetY];
				r_f_d[ind + 6 * ms] = r_f_d[ind + 6 * ms + offsetY];

				b_f_d[ind + 2 * ms] = b_f_d[ind + 2 * ms + offsetY];
				b_f_d[ind + 5 * ms] = b_f_d[ind + 5 * ms + offsetY];
				b_f_d[ind + 6 * ms] = b_f_d[ind + 6 * ms + offsetY];
			}
			break;
		case 3: //EAST
			r_f_d[ind + 3 * ms] = r_f_d[ind + 3 * ms - offsetX];
			r_f_d[ind + 7 * ms] = r_f_d[ind + 7 * ms - offsetX];
			r_f_d[ind + 6 * ms] = r_f_d[ind + 6 * ms - offsetX];

			b_f_d[ind + 3 * ms] = b_f_d[ind + 3 * ms - offsetX];
			b_f_d[ind + 7 * ms] = b_f_d[ind + 7 * ms - offsetX];
			b_f_d[ind + 6 * ms] = b_f_d[ind + 6 * ms - offsetX];
			break;
		case 4: //WEST
			r_f_d[ind + 1 * ms] = r_f_d[ind + 1 * ms + offsetX];
			r_f_d[ind + 5 * ms] = r_f_d[ind + 5 * ms + offsetX];
			r_f_d[ind + 8 * ms] = r_f_d[ind + 8 * ms + offsetX];

			b_f_d[ind + 1 * ms] = b_f_d[ind + 1 * ms + offsetX];
			b_f_d[ind + 5 * ms] = b_f_d[ind + 5 * ms + offsetX];
			b_f_d[ind + 8 * ms] = b_f_d[ind + 8 * ms + offsetX];
			break;
		case 5: //NE
			r_f_d[ind + 3 * ms] = r_f_d[ind + 3 * ms - offsetX];
			r_f_d[ind + 4 * ms] = r_f_d[ind + 4 * ms - offsetX];
			r_f_d[ind + 7 * ms] = r_f_d[ind + 7 * ms - offsetX];

			b_f_d[ind + 3 * ms] = b_f_d[ind + 3 * ms - offsetX];
			b_f_d[ind + 4 * ms] = b_f_d[ind + 4 * ms - offsetX];
			b_f_d[ind + 7 * ms] = b_f_d[ind + 7 * ms - offsetX];
			break;
		case 6: //NW
			r_f_d[ind + 1 * ms] = r_f_d[ind + 1 * ms + offsetX];
			r_f_d[ind + 4 * ms] = r_f_d[ind + 4 * ms + offsetX];
			r_f_d[ind + 8 * ms] = r_f_d[ind + 8 * ms + offsetX];

			b_f_d[ind + 1 * ms] = b_f_d[ind + 1 * ms + offsetX];
			b_f_d[ind + 4 * ms] = b_f_d[ind + 4 * ms + offsetX];
			b_f_d[ind + 8 * ms] = b_f_d[ind + 8 * ms + offsetX];
			break;
		case 7: //SE
			r_f_d[ind + 2 * ms] = r_f_d[ind + 2 * ms - offsetX];
			r_f_d[ind + 3 * ms] = r_f_d[ind + 3 * ms - offsetX];
			r_f_d[ind + 6 * ms] = r_f_d[ind + 6 * ms - offsetX];

			b_f_d[ind + 2 * ms] = b_f_d[ind + 2 * ms - offsetX];
			b_f_d[ind + 3 * ms] = b_f_d[ind + 3 * ms - offsetX];
			b_f_d[ind + 6 * ms] = b_f_d[ind + 6 * ms - offsetX];
			break;
		case 8: //SW
			r_f_d[ind + 2 * ms] = r_f_d[ind + 2 * ms + offsetX];
			r_f_d[ind + 1 * ms] = r_f_d[ind + 1 * ms + offsetX];
			r_f_d[ind + 5 * ms] = r_f_d[ind + 5 * ms + offsetX];

			b_f_d[ind + 2 * ms] = b_f_d[ind + 2 * ms + offsetX];
			b_f_d[ind + 1 * ms] = b_f_d[ind + 1 * ms + offsetX];
			b_f_d[ind + 5 * ms] = b_f_d[ind + 5 * ms + offsetX];
			break;
		}
	}
}

__global__ void gpuBcPeriodic3D(int *bcIdx_d, unsigned long long *bcMask_d,
		FLOAT_TYPE* f_d, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int bci = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																									+ threadIdx.x;
	int ms = depth_d * length_d * height_d;
	int offsetX = length_d - 1; //to get from west to east pairs of periodic nodes
	int offsetY = length_d * (depth_d - 1); //to get from north to south pairs of periodic nodes
	int offsetZ = length_d * depth_d * (height_d - 1); //to get from north to south pairs of periodic nodes

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC3D_FLUID) {
			if (bcMask_d[bci] & BC3D_CYCL_B) {
				/*WEST*/if ((bcMask_d[bci] & BC3D_CYCL_2)) {
					f_d[ind + 1 * ms] = f_d[ind + 1 * ms + offsetX];
					f_d[ind + 7 * ms] = f_d[ind + 7 * ms + offsetX];
					f_d[ind + 9 * ms] = f_d[ind + 9 * ms + offsetX];
					f_d[ind + 11 * ms] = f_d[ind + 11 * ms + offsetX];
					f_d[ind + 13 * ms] = f_d[ind + 13 * ms + offsetX];
				}
				/*EAST*/if ((bcMask_d[bci] & BC3D_CYCL_1)) {
					f_d[ind + 2 * ms] = f_d[ind + 2 * ms - offsetX];
					f_d[ind + 8 * ms] = f_d[ind + 8 * ms - offsetX];
					f_d[ind + 10 * ms] = f_d[ind + 10 * ms - offsetX];
					f_d[ind + 12 * ms] = f_d[ind + 12 * ms - offsetX];
					f_d[ind + 14 * ms] = f_d[ind + 14 * ms - offsetX];
				}

				/*SOUTH*/if ((bcMask_d[bci] & BC3D_CYCL_4)) {
					f_d[ind + 3 * ms] = f_d[ind + 3 * ms + offsetY];
					f_d[ind + 7 * ms] = f_d[ind + 7 * ms + offsetY];
					f_d[ind + 8 * ms] = f_d[ind + 8 * ms + offsetY];
					f_d[ind + 15 * ms] = f_d[ind + 15 * ms + offsetY];
					f_d[ind + 17 * ms] = f_d[ind + 17 * ms + offsetY];
				}
				/*NORTH*/if ((bcMask_d[bci] & BC3D_CYCL_3)) {
					f_d[ind + 4 * ms] = f_d[ind + 4 * ms - offsetY];
					f_d[ind + 9 * ms] = f_d[ind + 9 * ms - offsetY];
					f_d[ind + 10 * ms] = f_d[ind + 10 * ms - offsetY];
					f_d[ind + 16 * ms] = f_d[ind + 16 * ms - offsetY];
					f_d[ind + 18 * ms] = f_d[ind + 18 * ms - offsetY];
				}
				/*BOTTOM*/if ((bcMask_d[bci] & BC3D_CYCL_6)) {
					f_d[ind + 5 * ms] = f_d[ind + 5 * ms + offsetZ];
					f_d[ind + 11 * ms] = f_d[ind + 11 * ms + offsetZ];
					f_d[ind + 12 * ms] = f_d[ind + 12 * ms + offsetZ];
					f_d[ind + 15 * ms] = f_d[ind + 15 * ms + offsetZ];
					f_d[ind + 16 * ms] = f_d[ind + 16 * ms + offsetZ];
				}
				/*TOP*/if ((bcMask_d[bci] & BC3D_CYCL_5)) {
					f_d[ind + 6 * ms] = f_d[ind + 6 * ms - offsetZ];
					f_d[ind + 13 * ms] = f_d[ind + 13 * ms - offsetZ];
					f_d[ind + 14 * ms] = f_d[ind + 14 * ms - offsetZ];
					f_d[ind + 17 * ms] = f_d[ind + 17 * ms - offsetZ];
					f_d[ind + 18 * ms] = f_d[ind + 18 * ms - offsetZ];
				}
				if (bcMask_d[bci] & BC3D_CORNER) {
					printf("CORNER\n");
				}
			}
		}
	}
}
__global__ void gpuBcSymm3D(int *bcIdx_d, unsigned long long *bcMask_d,
		FLOAT_TYPE* f_d, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int bci = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																									+ threadIdx.x;
	int ms = depth_d * length_d * height_d;

	if (bci < size) {
		int ind = bcIdx_d[bci];
		if (bcMask_d[bci] & BC3D_FLUID) {
			if (bcMask_d[bci] & BC3D_SYMP_B) {
				/*WEST*/if ((bcMask_d[bci] & BC3D_SYMP_2)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 2))
								== ( BC3D_SYMP_2))) {
					//	        		printf("i am here Per W\n");
					f_d[ind + 1 * ms] = f_d[ind + 2 * ms];
					f_d[ind + 7 * ms] = f_d[ind + 8 * ms];
					f_d[ind + 9 * ms] = f_d[ind + 10 * ms];
					f_d[ind + 11 * ms] = f_d[ind + 12 * ms];
					f_d[ind + 13 * ms] = f_d[ind + 14 * ms];
				}
				/*EAST*/if ((bcMask_d[bci] & BC3D_SYMP_1)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 1))
								== ( BC3D_SYMP_1))) {
					//		printf("i am here Per E\n");
					f_d[ind + 2 * ms] = f_d[ind + 1 * ms];
					f_d[ind + 8 * ms] = f_d[ind + 7 * ms];
					f_d[ind + 10 * ms] = f_d[ind + 9 * ms] ;
					f_d[ind + 12 * ms] = f_d[ind + 11 * ms];
					f_d[ind + 14 * ms] = f_d[ind + 13 * ms];
				}

				/*SOUTH*/if ((bcMask_d[bci] & BC3D_SYMP_4)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 4))
								== ( BC3D_SYMP_4))) {
					//    		printf("i am here Symm S\n");
					f_d[ind + 3 * ms] = f_d[ind + 4*ms];
					f_d[ind + 7 * ms] = f_d[ind + 9*ms];
					f_d[ind + 8 * ms] = f_d[ind + 10*ms];
					f_d[ind + 15 * ms] = f_d[ind + 16*ms];
					f_d[ind + 17 * ms] = f_d[ind + 18 * ms];
				}
				/*NORTH*/if ((bcMask_d[bci] & BC3D_SYMP_3)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 3))
								== ( BC3D_SYMP_3))) {
					//    		printf("i am here Symm N\n");
					f_d[ind + 4 * ms] = f_d[ind + 3 * ms];
					f_d[ind + 9 * ms] = f_d[ind + 7 * ms];
					f_d[ind + 10 * ms] = f_d[ind + 8* ms];
					f_d[ind + 16 * ms] = f_d[ind + 15 * ms];
					f_d[ind + 18 * ms] = f_d[ind + 17 * ms];
				}
				/*BOTTOM*/if ((bcMask_d[bci] & BC3D_SYMP_6)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 6))
								== ( BC3D_SYMP_6))) {
					//   		printf("i am here Per B\n");
					f_d[ind + 5 * ms] = f_d[ind + 6*ms];
					f_d[ind + 11 * ms] = f_d[ind + 13*ms];
					f_d[ind + 12 * ms] = f_d[ind + 14*ms];
					f_d[ind + 15 * ms] = f_d[ind + 17*ms];
					f_d[ind + 16 * ms] = f_d[ind + 18*ms];
				}
				/*TOP*/if ((bcMask_d[bci] & BC3D_SYMP_5)
						&& ((bcMask_d[bci]
						              & BC3D_MASK((unsigned long long)BC3D_ALL, 5))
								== ( BC3D_SYMP_5))) {
					//        		printf("i am here Per T\n");//TODO
					f_d[ind + 6 * ms] = f_d[ind + 5 * ms];
					f_d[ind + 13 * ms] = f_d[ind + 11 * ms];
					f_d[ind + 14 * ms] = f_d[ind + 12 * ms];
					f_d[ind + 17 * ms] = f_d[ind + 15 * ms];
					f_d[ind + 18 * ms] = f_d[ind + 16 * ms];
				}
			}
		}
	}
}

