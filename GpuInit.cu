#include <stdio.h>
#include <complex.h>
#include <cuComplex.h>
#include "GpuFunctions.h"
#include "CellFunctions.h"
#include "ArrayUtils.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include <cmath>
#include "cuda.h"
#include "math.h"

__constant__ InletProfile inletProfile_d;
__constant__ BoundaryType boundaryType_d;
__constant__ OutletProfile outletProfile_d;
__constant__ int dlBoundaryId_d;
__constant__ int cx2D_d[9];
__constant__ int cy2D_d[9];
__constant__ int length_d;
__constant__ int depth_d;
__constant__ int height_d;
__constant__ int c2D_d[9];
__constant__ int opp2D_d[9];
__constant__ FLOAT_TYPE delta_d;
__constant__ FLOAT_TYPE w2D_d[9];
__constant__ FLOAT_TYPE omega_d;
__constant__ FLOAT_TYPE omegaA_d;
__constant__ FLOAT_TYPE rhoIn_d;
__constant__ FLOAT_TYPE uIn_d;
__constant__ FLOAT_TYPE vIn_d;
__constant__ FLOAT_TYPE minInletCoordY_d;
__constant__ FLOAT_TYPE maxInletCoordY_d;
__constant__ FLOAT_TYPE velMomMap2D_d[81];
__constant__ FLOAT_TYPE momCollMtx2D_d[81];
//#### 3D d3q19 ####//
__constant__ int cx3D_d[19];
__constant__ int cy3D_d[19];
__constant__ int cz3D_d[19];
__constant__ int c3D_d[19];
__constant__ int opp3D_d[19];
__constant__ FLOAT_TYPE w3D_d[19];
__constant__ FLOAT_TYPE minInletCoordZ_d;
__constant__ FLOAT_TYPE maxInletCoordZ_d;
__constant__ FLOAT_TYPE velMomMap3D_d[361];
__constant__ FLOAT_TYPE momCollMtx3D_d[361];
__constant__ FLOAT_TYPE wIn_d;
__constant__ FLOAT_TYPE g_d;

//COLOR GRADIENT //
__constant__ FLOAT_TYPE beta_d;
__constant__ FLOAT_TYPE r_alpha_d;
__constant__ FLOAT_TYPE b_alpha_d;
__constant__ FLOAT_TYPE bubble_radius_d;
__constant__ FLOAT_TYPE A_d;
__constant__ FLOAT_TYPE r_density_d;
__constant__ FLOAT_TYPE b_density_d;
__constant__ bool external_force_d;

//COLOR GRADIENT 2D//
__constant__ FLOAT_TYPE control_param_d;
__constant__ FLOAT_TYPE phi_d[9];
__constant__ FLOAT_TYPE teta_d[9];
__constant__ FLOAT_TYPE chi_d[9];
__constant__ FLOAT_TYPE psi_d[9];
__constant__ FLOAT_TYPE w_pert_d[9];
__constant__ FLOAT_TYPE g_limit_d;
__constant__ FLOAT_TYPE c_norms_d[9];
__constant__ FLOAT_TYPE cg_w_d[9];
__constant__ FLOAT_TYPE hocg_w_d[25];
__constant__ int hocg_cx_d[25];
__constant__ int hocg_cy_d[25];

//COLOR GRADIENT 3D//
__constant__ FLOAT_TYPE r_viscosity_d;
__constant__ FLOAT_TYPE b_viscosity_d;
__constant__ FLOAT_TYPE c_norms3D_d[19];
__constant__ FLOAT_TYPE w_pert3D_d[19];
__constant__ FLOAT_TYPE phi3D_d[19];
__constant__ FLOAT_TYPE teta3D_d[19];
__constant__ FLOAT_TYPE chi3D_d[19];
__constant__ FLOAT_TYPE psi3D_d[19];
__constant__ FLOAT_TYPE cg_w3D_d[19];
__constant__ FLOAT_TYPE hocg_w3D_d[105];
__constant__ int hocg_cx3D_d[105];
__constant__ int hocg_cy3D_d[105];
__constant__ int hocg_cz3D_d[105];
__constant__ int hoc3D_d[105];

__host__ void initConstants2D(Arguments *args,
		FLOAT_TYPE maxInletCoordY, FLOAT_TYPE minInletCoordY,
		FLOAT_TYPE delta, int m, int n) {
	//CONSTANT LATTICE QUANTITIES d2q9
	int s = m * n;
	FLOAT_TYPE w2D[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36.,
			1. / 36., 1. / 36., 1. / 36. };
	int opp2D[9] = { 0, 3 * s, 4 * s, 1 * s, 2 * s, 7 * s, 8 * s, 5 * s, 6 * s };
	int cx2D[9] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
	int cy2D[9] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
	int c2D[9] = { 0, -1, -1 * n, 1, n, -1 * n - 1, -1 * n + 1, n + 1, n - 1 };

	// Calculate collision freq
	FLOAT_TYPE omega = 1.0 / (3. * args->viscosity + 0.5);
	FLOAT_TYPE omegaA = 8 * (2 - omega) / (8 - omega);

	cudaMemcpyToSymbol(cx2D_d, cx2D, 9 * sizeof(int));
	cudaMemcpyToSymbol(cy2D_d, cy2D, 9 * sizeof(int));
	cudaMemcpyToSymbol(c2D_d, c2D, 9 * sizeof(int));
	cudaMemcpyToSymbol(opp2D_d, opp2D, 9 * sizeof(int));

	cudaMemcpyToSymbol(outletProfile_d, &args->outletProfile,
			sizeof(OutletProfile));
	cudaMemcpyToSymbol(boundaryType_d, &args->boundaryType,
			sizeof(BoundaryType));
	cudaMemcpyToSymbol(dlBoundaryId_d, &args->boundaryId, sizeof(int));

	cudaMemcpyToSymbol(depth_d, &m, sizeof(int));
	cudaMemcpyToSymbol(length_d, &n, sizeof(int));
	cudaMemcpyToSymbol(w2D_d, w2D, 9 * sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(omega_d, &omega, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(omegaA_d, &omegaA, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(delta_d, &delta, sizeof(FLOAT_TYPE));

	cudaMemcpyToSymbol(inletProfile_d, &args->inletProfile,
			sizeof(InletProfile));
	cudaMemcpyToSymbol(rhoIn_d, &args->rho, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(uIn_d, &args->u, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(vIn_d, &args->v, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(minInletCoordY_d, &minInletCoordY, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(maxInletCoordY_d, &maxInletCoordY, sizeof(FLOAT_TYPE));

	// Initialize variables for MRT Collision model, if used
	if (args->collisionModel == MRT) {
		FLOAT_TYPE *velMomMap2D = createHostArrayFlt(81);
		FLOAT_TYPE *momCollMtx2D = createHostArrayFlt(81);
		MRTInitializer2D(velMomMap2D, momCollMtx2D, omega);

		cudaMemcpyToSymbol(velMomMap2D_d, velMomMap2D, 81 * sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(momCollMtx2D_d, momCollMtx2D,
				81 * sizeof(FLOAT_TYPE));
	}

	cudaMemcpyToSymbol(g_d, &args->g, sizeof(FLOAT_TYPE));

	if (args->multiPhase){

		cudaMemcpyToSymbol(control_param_d, &args->control_param, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(beta_d, &args->beta, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(A_d, &args->A, sizeof(FLOAT_TYPE));

		FLOAT_TYPE aux1 = 1.0 / 5.0;
		FLOAT_TYPE aux2 = 1.0 /20.0;
		FLOAT_TYPE phi[9] = {0.0, aux1, aux1, aux1, aux1, aux2, aux2,aux2, aux2};
		cudaMemcpyToSymbol(phi_d, phi, 9 * sizeof(FLOAT_TYPE));
		aux1 = -1.0 / 5.0;
		aux2 = -1.0 / 20.0;
		FLOAT_TYPE teta[9] = {1.0, aux1, aux1, aux1, aux1, aux2, aux2,aux2, aux2};
		cudaMemcpyToSymbol(teta_d, teta, 9 * sizeof(FLOAT_TYPE));
		aux1 = -1.0 / 6.0;
		aux2 = 1.0 / 12.0;
		FLOAT_TYPE chi[9] = {-8.0 / 3.0, aux1, aux1, aux1, aux1, aux2, aux2,aux2, aux2};
//		FLOAT_TYPE chi[9] = {0.0};
		cudaMemcpyToSymbol(chi_d, chi, 9 * sizeof(FLOAT_TYPE));
		aux1 = 1.0 / 2.0;
		aux2 = 1.0 / 8.0;
		FLOAT_TYPE psi[9] = {0.0, aux1, aux1, aux1, aux1, aux2, aux2,aux2, aux2};
//		FLOAT_TYPE psi[9] = {0.0};
		cudaMemcpyToSymbol(psi_d, psi, 9 * sizeof(FLOAT_TYPE));

		FLOAT_TYPE w_pert[9] = {-4.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
				5.0 / 108.0, 5.0 / 108.0, 5.0 / 108.0, 5.0 / 108.0};
		cudaMemcpyToSymbol(w_pert_d, w_pert, 9 * sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(g_limit_d, &args->g_limit, sizeof(FLOAT_TYPE));

		FLOAT_TYPE c_norms[9];
		for(int i = 0; i < 9;i++){
			c_norms[i] = sqrt(cx2D[i] * cx2D[i] + cy2D[i] * cy2D[i]);
		}
		cudaMemcpyToSymbol(c_norms_d, c_norms, 9 * sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(r_density_d, &args->r_density, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_density_d, &args->b_density, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(r_alpha_d, &args->r_alpha, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_alpha_d, &args->b_alpha, sizeof(FLOAT_TYPE));
		args->bubble_radius /= n;
		cudaMemcpyToSymbol(bubble_radius_d, &args->bubble_radius, sizeof(FLOAT_TYPE));
//		FLOAT_TYPE st_predicted = (2.0/9.0)*(1.0+1.0/args->gamma)/(0.5*(r_omega+b_omega))*0.5*args->r_density*(args->r_A+args->b_A);
//		cudaMemcpyToSymbol(st_predicted_d, &st_predicted, sizeof(FLOAT_TYPE));

		FLOAT_TYPE cg_w[9] = {0., 4. / 12., 4. / 12., 4. / 12., 4. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12.};
		cudaMemcpyToSymbol(cg_w_d, cg_w, 9 * sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(r_viscosity_d, &args->r_viscosity, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_viscosity_d, &args->b_viscosity, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(external_force_d, &args->external_force, sizeof(bool));

		FLOAT_TYPE hocg_w[25] = {0., 960. / 5040., 960. / 5040., 960. / 5040., 960. / 5040., 448. / 5040., 448. / 5040.,
				448. / 5040., 448. / 5040., 84. / 5040., 32. / 5040., 1. / 5040., 32. / 5040., 84. / 5040., 32. / 5040.,
				1. / 5040., 32. / 5040., 84. / 5040., 32. / 5040., 1. / 5040., 32. / 5040., 84. / 5040., 32. / 5040., 1. / 5040., 32. / 5040.};
		cudaMemcpyToSymbol(hocg_w_d, hocg_w, 25 * sizeof(FLOAT_TYPE));
		int hocg_cx[25] = {0,1,0,-1,0,1,-1,-1,1,0,1,2,2,2,2,2,1,0,-1,-2,-2,-2,-2,-2,-1};
		cudaMemcpyToSymbol(hocg_cx_d, hocg_cx, 25 * sizeof(int));
		int hocg_cy[25] = {0,0,1,0,-1,1,1,-1,-1,2,2,2,1,0,-1,-2,-2,-2,-2,-2,-1,0,1,2,2};
		cudaMemcpyToSymbol(hocg_cy_d, hocg_cy, 25 * sizeof(int));

		cudaMemcpyToSymbol(r_viscosity_d, &args->r_viscosity, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_viscosity_d, &args->b_viscosity, sizeof(FLOAT_TYPE));
	}
}

__host__ void initHOColorGradient(int *color_gradient_directions, int n, int m){
	int jn = m-1;
	int js = 0;
	int ie = n-1;
	int iw = 0;
	for(int j = 0; j < m;j++){
		for(int i = 0; i < n; i++){

			if(j == m-1  && i != 0 && i != n-1) //NORTH
				color_gradient_directions[j * n + i] = 1;
			else if(j == 0 && i != 0 && i != n-1 ) //SOUTH
				color_gradient_directions[j * n + i] = 2;
			else if(i == n - 1  && j != 0 && j != m-1 ) //EAST
				color_gradient_directions[j * n + i] = 3;
			else if(i == 0 && j != 0 && j != m-1 ) //WEST
				color_gradient_directions[j * n + i] = 4;
			else if(j == 1 || j == m-2  || i == 1 || i == n-2) //Neighbour boundary
				color_gradient_directions[j * n + i] = 11;
		}
	}

	//CORNERS
	//Probably should handle corners
	color_gradient_directions[jn * n + ie] = 5; //NE
	color_gradient_directions[jn * n + iw] = 6; //NW
	color_gradient_directions[js * n + ie] = 7; //SE
	color_gradient_directions[js * n + iw] = 8; //SW
}

__host__ void initColorGradient(int *color_gradient_directions, int n, int m){
	//NORTH is 1, SOUTH 2, EAST 3, WEST 4, 0 is okay
	int jn = m-1;
	int js = 0;
	int ie = n-1;
	int iw = 0;
	for(int j = 0; j < m;j++){
		for(int i = 0; i < n; i++){

			if(j == m-1  && i != 0 && i != n-1) //NORTH
				color_gradient_directions[j * n + i] = 1;
			else if(j == 0 && i != 0 && i != n-1 ) //SOUTH
				color_gradient_directions[j * n + i] = 2;
			else if(i == n - 1  && j != 0 && j != m-1 ) //EAST
				color_gradient_directions[j * n + i] = 3;
			else if(i == 0 && j != 0 && j != m-1 ) //WEST
				color_gradient_directions[j * n + i] = 4;
		}
	}

	//CORNERS
	color_gradient_directions[jn * n + ie] = 5; //NE
	color_gradient_directions[jn * n + iw] = 6; //NW
	color_gradient_directions[js * n + ie] = 7; //SE
	color_gradient_directions[js * n + iw] = 8; //SW
}

__host__ void initColorGradient3D(int *color_gradient_directions, int n, int m, int h){
	//NORTH is 1, SOUTH 2, EAST 3, WEST 4, FRONT 5, BACK 6, 0 is okay
	int ms = n * m;

	for(int k = 0; k < h; k++){
		for(int j = 0; j < m;j++){
			for(int i = 0; i < n; i++){
				if(j != m-1 && j != 0 && i != 0 && i != n-1 && k != 0 && k != h-1)
					color_gradient_directions[k * ms + j * n + i] = 0;
				else if(j == m-1 && i != 0 && i != n-1 && k != 0 && k != h-1) //NORTH
					color_gradient_directions[k * ms + j * n + i] = 1;
				else if(j == 0 && i != 0 && i != n-1 && k != 0 && k != h-1) //SOUTH
					color_gradient_directions[k * ms + j * n + i] = 2;
				else if(i == n - 1 && j != 0 && j != m-1 && k != 0 && k != h-1) //EAST
					color_gradient_directions[k * ms + j * n + i] = 3;
				else if(i == 0 && j != 0 && j != m-1 && k != 0 && k != h-1 ) //WEST
					color_gradient_directions[k * ms + j * n + i] = 4;
				else if(k == h-1 && j != 0 && j != m-1 && i != 0 && i != n-1) //FRONT
					color_gradient_directions[k * ms + j * n + i] = 5;
				else if(k == 0 && j != 0 && j != m-1 && i != 0 && i != n-1)  //BACK
					color_gradient_directions[k * ms + j * n + i] = 6;
				else { //edges and corners
					color_gradient_directions[k * ms + j * n + i] = -1; // corners
				}
			}
		}
	}
}

__host__ void initHOColorGradient3D(int *color_gradient_directions, int n, int m, int h){

	//NORTH is 1, SOUTH 2, EAST 3, WEST 4, FRONT 5, BACK 6, 0 is okay
	int ms = n * m;

	for(int k = 0; k < h; k++){
		for(int j = 0; j < m;j++){
			for(int i = 0; i < n; i++){
				if(j == 1 || j == m-2  || i == 1 || i == n-2|| k == 1 || k == h-2) //Neighbour boundary
					color_gradient_directions[k * ms + j * n + i] = 11;
				else if(j != m-1 && j != 0 && i != 0 && i != n-1 && k != 0 && k != h-1)
					color_gradient_directions[k * ms + j * n + i] = 0;
				else if(j == m-1 && i != 0 && i != n-1 && k != 0 && k != h-1) //NORTH
					color_gradient_directions[k * ms + j * n + i] = 1;
				else if(j == 0 && i != 0 && i != n-1 && k != 0 && k != h-1) //SOUTH
					color_gradient_directions[k * ms + j * n + i] = 2;
				else if(i == n - 1 && j != 0 && j != m-1 && k != 0 && k != h-1) //EAST
					color_gradient_directions[k * ms + j * n + i] = 3;
				else if(i == 0 && j != 0 && j != m-1 && k != 0 && k != h-1 ) //WEST
					color_gradient_directions[k * ms + j * n + i] = 4;
				else if(k == h-1 && j != 0 && j != m-1 && i != 0 && i != n-1) //FRONT
					color_gradient_directions[k * ms + j * n + i] = 5;
				else if(k == 0 && j != 0 && j != m-1 && i != 0 && i != n-1)  //BACK
					color_gradient_directions[k * ms + j * n + i] = 6;
				else { //edges and corners
					color_gradient_directions[k * ms + j * n + i] = -1; // corners
				}
			}
		}
	}
}

__global__ void initCGBubble(FLOAT_TYPE *x_d, FLOAT_TYPE *y_d, FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *r_f_d,
		FLOAT_TYPE *b_f_d, FLOAT_TYPE *f_d, int test_case){
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int ms = length_d * depth_d;

	if(index < ms){
		FLOAT_TYPE aux1, aux2;
		int index_x, index_y;
		switch (test_case) {
		case 1: //steady bubble
			if( sqrt( (x_d[index]-0.5) * (x_d[index]-0.5) + (y_d[index]-0.5)*(y_d[index]-0.5)) <= bubble_radius_d){
				aux1 = (1 - r_alpha_d) / 5.0;
				aux2 = (1 - r_alpha_d) / 20.0;
				r_rho_d[index] = r_density_d;
				r_f_d[index + 0 * ms] = r_density_d * r_alpha_d;
				r_f_d[index + 1 * ms] = r_density_d * aux1;
				r_f_d[index + 2 * ms] = r_density_d * aux1;
				r_f_d[index + 3 * ms] = r_density_d * aux1;
				r_f_d[index + 4 * ms] = r_density_d * aux1;
				r_f_d[index + 5 * ms] = r_density_d * aux2;
				r_f_d[index + 6 * ms] = r_density_d * aux2;
				r_f_d[index + 7 * ms] = r_density_d * aux2;
				r_f_d[index + 8 * ms] = r_density_d * aux2;
			}
			else {
				aux1 = (1 - b_alpha_d) / 5.0;
				aux2 = (1 - b_alpha_d) / 20.0;
				b_rho_d[index] = b_density_d;
				b_f_d[index + 0 * ms] = b_density_d * b_alpha_d;
				b_f_d[index + 1 * ms] = b_density_d * aux1;
				b_f_d[index + 2 * ms] = b_density_d * aux1;
				b_f_d[index + 3 * ms] = b_density_d * aux1;
				b_f_d[index + 4 * ms] = b_density_d * aux1;
				b_f_d[index + 5 * ms] = b_density_d * aux2;
				b_f_d[index + 6 * ms] = b_density_d * aux2;
				b_f_d[index + 7 * ms] = b_density_d * aux2;
				b_f_d[index + 8 * ms] = b_density_d * aux2;
			}
			break;
		case 2: // coulette
			if(y_d[index] > 0.5){
				aux1 = (1 - r_alpha_d) / 5.0;
				aux2 = (1 - r_alpha_d) / 20.0;
				r_rho_d[index] = r_density_d;
				r_f_d[index + 0 * ms] = r_density_d * r_alpha_d;
				r_f_d[index + 1 * ms] = r_density_d * aux1;
				r_f_d[index + 2 * ms] = r_density_d * aux1;
				r_f_d[index + 3 * ms] = r_density_d * aux1;
				r_f_d[index + 4 * ms] = r_density_d * aux1;
				r_f_d[index + 5 * ms] = r_density_d * aux2;
				r_f_d[index + 6 * ms] = r_density_d * aux2;
				r_f_d[index + 7 * ms] = r_density_d * aux2;
				r_f_d[index + 8 * ms] = r_density_d * aux2;
			}
			else {
				aux1 = (1 - b_alpha_d) / 5.0;
				aux2 = (1 - b_alpha_d) / 20.0;
				b_rho_d[index] = b_density_d;
				b_f_d[index + 0 * ms] = b_density_d * b_alpha_d;
				b_f_d[index + 1 * ms] = b_density_d * aux1;
				b_f_d[index + 2 * ms] = b_density_d * aux1;
				b_f_d[index + 3 * ms] = b_density_d * aux1;
				b_f_d[index + 4 * ms] = b_density_d * aux1;
				b_f_d[index + 5 * ms] = b_density_d * aux2;
				b_f_d[index + 6 * ms] = b_density_d * aux2;
				b_f_d[index + 7 * ms] = b_density_d * aux2;
				b_f_d[index + 8 * ms] = b_density_d * aux2;
			}
			break;
		case 3: //square
			index_x = index % length_d;
			index_y = (index - index_x) / length_d;

			if( index_x < 0.75 * length_d && index_x > 0.25 * length_d && index_y < 0.75 * depth_d && index_y > 0.25 * depth_d){
				aux1 = (1 - r_alpha_d) / 5.0;
				aux2 = (1 - r_alpha_d) / 20.0;
				r_rho_d[index] = r_density_d;
				r_f_d[index + 0 * ms] = r_density_d * r_alpha_d;
				r_f_d[index + 1 * ms] = r_density_d * aux1;
				r_f_d[index + 2 * ms] = r_density_d * aux1;
				r_f_d[index + 3 * ms] = r_density_d * aux1;
				r_f_d[index + 4 * ms] = r_density_d * aux1;
				r_f_d[index + 5 * ms] = r_density_d * aux2;
				r_f_d[index + 6 * ms] = r_density_d * aux2;
				r_f_d[index + 7 * ms] = r_density_d * aux2;
				r_f_d[index + 8 * ms] = r_density_d * aux2;
			}
			else {
				aux1 = (1 - b_alpha_d) / 5.0;
				aux2 = (1 - b_alpha_d) / 20.0;
				b_rho_d[index] = b_density_d;
				b_f_d[index + 0 * ms] = b_density_d * b_alpha_d;
				b_f_d[index + 1 * ms] = b_density_d * aux1;
				b_f_d[index + 2 * ms] = b_density_d * aux1;
				b_f_d[index + 3 * ms] = b_density_d * aux1;
				b_f_d[index + 4 * ms] = b_density_d * aux1;
				b_f_d[index + 5 * ms] = b_density_d * aux2;
				b_f_d[index + 6 * ms] = b_density_d * aux2;
				b_f_d[index + 7 * ms] = b_density_d * aux2;
				b_f_d[index + 8 * ms] = b_density_d * aux2;
			}
			break;
		case 4: //coalescence
			if( sqrt( (x_d[index]-0.5) * (x_d[index]-0.5) + (y_d[index]-0.5 + bubble_radius_d)*(y_d[index]-0.5 + bubble_radius_d)) <= bubble_radius_d ||
					sqrt( (x_d[index]-0.5) * (x_d[index]-0.5) + (y_d[index]-0.5 - bubble_radius_d)*(y_d[index]-0.5 - bubble_radius_d)) <= bubble_radius_d	){
				aux1 = (1 - r_alpha_d) / 5.0;
				aux2 = (1 - r_alpha_d) / 20.0;
				r_rho_d[index] = r_density_d;
				r_f_d[index + 0 * ms] = r_density_d * r_alpha_d;
				r_f_d[index + 1 * ms] = r_density_d * aux1;
				r_f_d[index + 2 * ms] = r_density_d * aux1;
				r_f_d[index + 3 * ms] = r_density_d * aux1;
				r_f_d[index + 4 * ms] = r_density_d * aux1;
				r_f_d[index + 5 * ms] = r_density_d * aux2;
				r_f_d[index + 6 * ms] = r_density_d * aux2;
				r_f_d[index + 7 * ms] = r_density_d * aux2;
				r_f_d[index + 8 * ms] = r_density_d * aux2;
			}
			else {
				aux1 = (1 - b_alpha_d) / 5.0;
				aux2 = (1 - b_alpha_d) / 20.0;
				b_rho_d[index] = b_density_d;
				b_f_d[index + 0 * ms] = b_density_d * b_alpha_d;
				b_f_d[index + 1 * ms] = b_density_d * aux1;
				b_f_d[index + 2 * ms] = b_density_d * aux1;
				b_f_d[index + 3 * ms] = b_density_d * aux1;
				b_f_d[index + 4 * ms] = b_density_d * aux1;
				b_f_d[index + 5 * ms] = b_density_d * aux2;
				b_f_d[index + 6 * ms] = b_density_d * aux2;
				b_f_d[index + 7 * ms] = b_density_d * aux2;
				b_f_d[index + 8 * ms] = b_density_d * aux2;
			}
			break;
		case 5: //oscilating
			index_x = index % length_d;
			index_y = (index - index_x) / length_d;
			if( ( ((index_x - 0.5 * length_d) / (0.125 * length_d)) * ((index_x - 0.5 * length_d) / (0.125 * length_d)) +
					((index_y - 0.5 * depth_d) / (0.1875 * depth_d)) * ((index_y - 0.5 * depth_d) / (0.1875 * depth_d)) <= 1 )){
				aux1 = (1 - r_alpha_d) / 5.0;
				aux2 = (1 - r_alpha_d) / 20.0;
				r_rho_d[index] = r_density_d;
				r_f_d[index + 0 * ms] = r_density_d * r_alpha_d;
				r_f_d[index + 1 * ms] = r_density_d * aux1;
				r_f_d[index + 2 * ms] = r_density_d * aux1;
				r_f_d[index + 3 * ms] = r_density_d * aux1;
				r_f_d[index + 4 * ms] = r_density_d * aux1;
				r_f_d[index + 5 * ms] = r_density_d * aux2;
				r_f_d[index + 6 * ms] = r_density_d * aux2;
				r_f_d[index + 7 * ms] = r_density_d * aux2;
				r_f_d[index + 8 * ms] = r_density_d * aux2;
			}
			else {
				aux1 = (1 - b_alpha_d) / 5.0;
				aux2 = (1 - b_alpha_d) / 20.0;
				b_rho_d[index] = b_density_d;
				b_f_d[index + 0 * ms] = b_density_d * b_alpha_d;
				b_f_d[index + 1 * ms] = b_density_d * aux1;
				b_f_d[index + 2 * ms] = b_density_d * aux1;
				b_f_d[index + 3 * ms] = b_density_d * aux1;
				b_f_d[index + 4 * ms] = b_density_d * aux1;
				b_f_d[index + 5 * ms] = b_density_d * aux2;
				b_f_d[index + 6 * ms] = b_density_d * aux2;
				b_f_d[index + 7 * ms] = b_density_d * aux2;
				b_f_d[index + 8 * ms] = b_density_d * aux2;
			}
			break;
		case 6:
			if( y_d[index] > (2.0 + 0.1 * cos( 2*M_PI*x_d[index]))){
				aux1 = (1 - r_alpha_d) / 5.0;
				aux2 = (1 - r_alpha_d) / 20.0;
				r_rho_d[index] = r_density_d;
				r_f_d[index + 0 * ms] = r_density_d * r_alpha_d;
				r_f_d[index + 1 * ms] = r_density_d * aux1;
				r_f_d[index + 2 * ms] = r_density_d * aux1;
				r_f_d[index + 3 * ms] = r_density_d * aux1;
				r_f_d[index + 4 * ms] = r_density_d * aux1;
				r_f_d[index + 5 * ms] = r_density_d * aux2;
				r_f_d[index + 6 * ms] = r_density_d * aux2;
				r_f_d[index + 7 * ms] = r_density_d * aux2;
				r_f_d[index + 8 * ms] = r_density_d * aux2;
			}
			else {
				aux1 = (1 - b_alpha_d) / 5.0;
				aux2 = (1 - b_alpha_d) / 20.0;
				b_rho_d[index] = b_density_d;
				b_f_d[index + 0 * ms] = b_density_d * b_alpha_d;
				b_f_d[index + 1 * ms] = b_density_d * aux1;
				b_f_d[index + 2 * ms] = b_density_d * aux1;
				b_f_d[index + 3 * ms] = b_density_d * aux1;
				b_f_d[index + 4 * ms] = b_density_d * aux1;
				b_f_d[index + 5 * ms] = b_density_d * aux2;
				b_f_d[index + 6 * ms] = b_density_d * aux2;
				b_f_d[index + 7 * ms] = b_density_d * aux2;
				b_f_d[index + 8 * ms] = b_density_d * aux2;
			}
			break;
		default:
			break;
		}
		// initialise density
		rho_d[index] = r_rho_d[index] + b_rho_d[index];
		f_d[index] = r_f_d[index] + b_f_d[index];
		f_d[index + 1 * ms] = r_f_d[index + 1 * ms] + b_f_d[index + 1 * ms];
		f_d[index + 2 * ms] = r_f_d[index + 2 * ms] + b_f_d[index + 2 * ms];
		f_d[index + 3 * ms] = r_f_d[index + 3 * ms] + b_f_d[index + 3 * ms];
		f_d[index + 4 * ms] = r_f_d[index + 4 * ms] + b_f_d[index + 4 * ms];
		f_d[index + 5 * ms] = r_f_d[index + 5 * ms] + b_f_d[index + 5 * ms];
		f_d[index + 6 * ms] = r_f_d[index + 6 * ms] + b_f_d[index + 6 * ms];
		f_d[index + 7 * ms] = r_f_d[index + 7 * ms] + b_f_d[index + 7 * ms];
		f_d[index + 8 * ms] = r_f_d[index + 8 * ms] + b_f_d[index + 8 * ms];
	}
}

__global__ void initCGBubble3D(FLOAT_TYPE *x_d, FLOAT_TYPE *y_d, FLOAT_TYPE *z_d, FLOAT_TYPE *r_rho_d, FLOAT_TYPE *b_rho_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *r_f_d,
		FLOAT_TYPE *b_f_d, FLOAT_TYPE *f_d, int test_case){
	int index =  (blockIdx.x + blockIdx.y * gridDim.x) * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x) + threadIdx.x;
	int ms = depth_d*length_d*height_d;
	if(index < ms){
		int index_x, index_y, index_z, temp;
		switch (test_case) {
		case 1: //steady bubble
			if( sqrt( (x_d[index]-0.5) * (x_d[index]-0.5) + (y_d[index]-0.5)*(y_d[index]-0.5) + (z_d[index]-0.5)*(z_d[index]-0.5)) <= bubble_radius_d){

				r_rho_d[index] = r_density_d;
				// initialise density
				rho_d[index] = r_density_d;
				//r_f_d[index] = 0; First r_f_d is 0
				f_d[index + 1 * ms] = r_f_d[index + 1 * ms] = r_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = r_f_d[index + 2 * ms] = r_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = r_f_d[index + 3 * ms] = r_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = r_f_d[index + 4 * ms] = r_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = r_f_d[index + 5 * ms] = r_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = r_f_d[index + 6 * ms] = r_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = r_f_d[index + 7 * ms] = r_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = r_f_d[index + 8 * ms] = r_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = r_f_d[index + 9 * ms] = r_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = r_f_d[index + 10 * ms] = r_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = r_f_d[index + 11 * ms] = r_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = r_f_d[index + 12 * ms] = r_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = r_f_d[index + 13 * ms] = r_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = r_f_d[index + 14 * ms] = r_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = r_f_d[index + 15 * ms] = r_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = r_f_d[index + 16 * ms] = r_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = r_f_d[index + 17 * ms] = r_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = r_f_d[index + 18 * ms] = r_density_d * phi3D_d[18];
			}
			else {
				b_rho_d[index] = b_density_d;
				// initialise density
				rho_d[index] = b_density_d;
				f_d[index + 1 * ms] = b_f_d[index + 1 * ms] = b_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = b_f_d[index + 2 * ms] = b_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = b_f_d[index + 3 * ms] = b_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = b_f_d[index + 4 * ms] = b_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = b_f_d[index + 5 * ms] = b_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = b_f_d[index + 6 * ms] = b_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = b_f_d[index + 7 * ms] = b_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = b_f_d[index + 8 * ms] = b_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = b_f_d[index + 9 * ms] = b_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = b_f_d[index + 10 * ms] = b_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = b_f_d[index + 11 * ms] = b_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = b_f_d[index + 12 * ms] = b_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = b_f_d[index + 13 * ms] = b_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = b_f_d[index + 14 * ms] = b_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = b_f_d[index + 15 * ms] = b_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = b_f_d[index + 16 * ms] = b_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = b_f_d[index + 17 * ms] = b_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = b_f_d[index + 18 * ms] = b_density_d * phi3D_d[18];
			}
			break;
		case 2: //SQUARE
			index_x = index % length_d;
			temp = (index - index_x) / length_d;
			index_y = temp % depth_d;
			index_z = (temp - index_y) / depth_d;

			if( index_x < 0.75 * length_d && index_x > 0.25 * length_d && index_y < 0.75 * depth_d && index_y > 0.25 * depth_d &&
					index_z < 0.75 * height_d && index_z > 0.25 * height_d){
				r_rho_d[index] = r_density_d;
				// initialise density
				rho_d[index] = r_density_d;
				//r_f_d[index] = 0; First r_f_d is 0
				f_d[index + 1 * ms] = r_f_d[index + 1 * ms] = r_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = r_f_d[index + 2 * ms] = r_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = r_f_d[index + 3 * ms] = r_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = r_f_d[index + 4 * ms] = r_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = r_f_d[index + 5 * ms] = r_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = r_f_d[index + 6 * ms] = r_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = r_f_d[index + 7 * ms] = r_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = r_f_d[index + 8 * ms] = r_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = r_f_d[index + 9 * ms] = r_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = r_f_d[index + 10 * ms] = r_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = r_f_d[index + 11 * ms] = r_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = r_f_d[index + 12 * ms] = r_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = r_f_d[index + 13 * ms] = r_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = r_f_d[index + 14 * ms] = r_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = r_f_d[index + 15 * ms] = r_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = r_f_d[index + 16 * ms] = r_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = r_f_d[index + 17 * ms] = r_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = r_f_d[index + 18 * ms] = r_density_d * phi3D_d[18];
			}
			else {
				b_rho_d[index] = b_density_d;
				// initialise density
				rho_d[index] = b_density_d;
				f_d[index + 1 * ms] = b_f_d[index + 1 * ms] = b_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = b_f_d[index + 2 * ms] = b_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = b_f_d[index + 3 * ms] = b_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = b_f_d[index + 4 * ms] = b_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = b_f_d[index + 5 * ms] = b_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = b_f_d[index + 6 * ms] = b_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = b_f_d[index + 7 * ms] = b_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = b_f_d[index + 8 * ms] = b_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = b_f_d[index + 9 * ms] = b_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = b_f_d[index + 10 * ms] = b_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = b_f_d[index + 11 * ms] = b_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = b_f_d[index + 12 * ms] = b_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = b_f_d[index + 13 * ms] = b_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = b_f_d[index + 14 * ms] = b_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = b_f_d[index + 15 * ms] = b_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = b_f_d[index + 16 * ms] = b_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = b_f_d[index + 17 * ms] = b_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = b_f_d[index + 18 * ms] = b_density_d * phi3D_d[18];
			}
			break;
		case 3: //coalescence
			if(sqrt( (x_d[index]-0.5) * (x_d[index]-0.5) + (y_d[index]-0.5 + bubble_radius_d)*(y_d[index]-0.5 + bubble_radius_d) + (z_d[index]-0.5)*(z_d[index]-0.5)) <= bubble_radius_d ||
					sqrt( (x_d[index]-0.5) * (x_d[index]-0.5) + (y_d[index]-0.5 - bubble_radius_d)*(y_d[index]-0.5 - bubble_radius_d) + (z_d[index]-0.5)*(z_d[index]-0.5)) <= bubble_radius_d	){
				r_rho_d[index] = r_density_d;
				// initialise density
				rho_d[index] = r_density_d;
				//r_f_d[index] = 0; First r_f_d is 0
				f_d[index + 1 * ms] = r_f_d[index + 1 * ms] = r_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = r_f_d[index + 2 * ms] = r_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = r_f_d[index + 3 * ms] = r_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = r_f_d[index + 4 * ms] = r_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = r_f_d[index + 5 * ms] = r_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = r_f_d[index + 6 * ms] = r_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = r_f_d[index + 7 * ms] = r_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = r_f_d[index + 8 * ms] = r_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = r_f_d[index + 9 * ms] = r_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = r_f_d[index + 10 * ms] = r_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = r_f_d[index + 11 * ms] = r_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = r_f_d[index + 12 * ms] = r_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = r_f_d[index + 13 * ms] = r_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = r_f_d[index + 14 * ms] = r_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = r_f_d[index + 15 * ms] = r_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = r_f_d[index + 16 * ms] = r_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = r_f_d[index + 17 * ms] = r_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = r_f_d[index + 18 * ms] = r_density_d * phi3D_d[18];
			}
			else {
				b_rho_d[index] = b_density_d;
				// initialise density
				rho_d[index] = b_density_d;
				f_d[index + 1 * ms] = b_f_d[index + 1 * ms] = b_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = b_f_d[index + 2 * ms] = b_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = b_f_d[index + 3 * ms] = b_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = b_f_d[index + 4 * ms] = b_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = b_f_d[index + 5 * ms] = b_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = b_f_d[index + 6 * ms] = b_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = b_f_d[index + 7 * ms] = b_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = b_f_d[index + 8 * ms] = b_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = b_f_d[index + 9 * ms] = b_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = b_f_d[index + 10 * ms] = b_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = b_f_d[index + 11 * ms] = b_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = b_f_d[index + 12 * ms] = b_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = b_f_d[index + 13 * ms] = b_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = b_f_d[index + 14 * ms] = b_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = b_f_d[index + 15 * ms] = b_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = b_f_d[index + 16 * ms] = b_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = b_f_d[index + 17 * ms] = b_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = b_f_d[index + 18 * ms] = b_density_d * phi3D_d[18];
			}
			break;
		case 4: //Couette
			if(y_d[index] > 0.5){
				r_rho_d[index] = r_density_d;
				// initialise density
				rho_d[index] = r_density_d;
				//r_f_d[index] = 0; First r_f_d is 0
				f_d[index + 1 * ms] = r_f_d[index + 1 * ms] = r_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = r_f_d[index + 2 * ms] = r_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = r_f_d[index + 3 * ms] = r_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = r_f_d[index + 4 * ms] = r_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = r_f_d[index + 5 * ms] = r_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = r_f_d[index + 6 * ms] = r_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = r_f_d[index + 7 * ms] = r_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = r_f_d[index + 8 * ms] = r_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = r_f_d[index + 9 * ms] = r_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = r_f_d[index + 10 * ms] = r_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = r_f_d[index + 11 * ms] = r_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = r_f_d[index + 12 * ms] = r_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = r_f_d[index + 13 * ms] = r_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = r_f_d[index + 14 * ms] = r_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = r_f_d[index + 15 * ms] = r_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = r_f_d[index + 16 * ms] = r_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = r_f_d[index + 17 * ms] = r_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = r_f_d[index + 18 * ms] = r_density_d * phi3D_d[18];
			}
			else {
				b_rho_d[index] = b_density_d;
				// initialise density
				rho_d[index] = b_density_d;
				f_d[index + 1 * ms] = b_f_d[index + 1 * ms] = b_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = b_f_d[index + 2 * ms] = b_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = b_f_d[index + 3 * ms] = b_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = b_f_d[index + 4 * ms] = b_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = b_f_d[index + 5 * ms] = b_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = b_f_d[index + 6 * ms] = b_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = b_f_d[index + 7 * ms] = b_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = b_f_d[index + 8 * ms] = b_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = b_f_d[index + 9 * ms] = b_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = b_f_d[index + 10 * ms] = b_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = b_f_d[index + 11 * ms] = b_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = b_f_d[index + 12 * ms] = b_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = b_f_d[index + 13 * ms] = b_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = b_f_d[index + 14 * ms] = b_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = b_f_d[index + 15 * ms] = b_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = b_f_d[index + 16 * ms] = b_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = b_f_d[index + 17 * ms] = b_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = b_f_d[index + 18 * ms] = b_density_d * phi3D_d[18];
			}
			break;
		case 5: //RT

			if( y_d[index] > (1 + 0.1 * cos( 2*M_PI*x_d[index]))){
				r_rho_d[index] = r_density_d;
				// initialise density
				rho_d[index] = r_density_d;
				//r_f_d[index] = 0; First r_f_d is 0
				f_d[index + 1 * ms] = r_f_d[index + 1 * ms] = r_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = r_f_d[index + 2 * ms] = r_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = r_f_d[index + 3 * ms] = r_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = r_f_d[index + 4 * ms] = r_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = r_f_d[index + 5 * ms] = r_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = r_f_d[index + 6 * ms] = r_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = r_f_d[index + 7 * ms] = r_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = r_f_d[index + 8 * ms] = r_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = r_f_d[index + 9 * ms] = r_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = r_f_d[index + 10 * ms] = r_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = r_f_d[index + 11 * ms] = r_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = r_f_d[index + 12 * ms] = r_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = r_f_d[index + 13 * ms] = r_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = r_f_d[index + 14 * ms] = r_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = r_f_d[index + 15 * ms] = r_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = r_f_d[index + 16 * ms] = r_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = r_f_d[index + 17 * ms] = r_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = r_f_d[index + 18 * ms] = r_density_d * phi3D_d[18];
			}
			else {
				b_rho_d[index] = b_density_d;
				// initialise density
				rho_d[index] = b_density_d;
				f_d[index + 1 * ms] = b_f_d[index + 1 * ms] = b_density_d * phi3D_d[1];
				f_d[index + 2 * ms] = b_f_d[index + 2 * ms] = b_density_d * phi3D_d[2];
				f_d[index + 3 * ms] = b_f_d[index + 3 * ms] = b_density_d * phi3D_d[3];
				f_d[index + 4 * ms] = b_f_d[index + 4 * ms] = b_density_d * phi3D_d[4];
				f_d[index + 5 * ms] = b_f_d[index + 5 * ms] = b_density_d * phi3D_d[5];
				f_d[index + 6 * ms] = b_f_d[index + 6 * ms] = b_density_d * phi3D_d[6];
				f_d[index + 7 * ms] = b_f_d[index + 7 * ms] = b_density_d * phi3D_d[7];
				f_d[index + 8 * ms] = b_f_d[index + 8 * ms] = b_density_d * phi3D_d[8];
				f_d[index + 9 * ms] = b_f_d[index + 9 * ms] = b_density_d * phi3D_d[9];
				f_d[index + 10 * ms] = b_f_d[index + 10 * ms] = b_density_d * phi3D_d[10];
				f_d[index + 11 * ms] = b_f_d[index + 11 * ms] = b_density_d * phi3D_d[11];
				f_d[index + 12 * ms] = b_f_d[index + 12 * ms] = b_density_d * phi3D_d[12];
				f_d[index + 13 * ms] = b_f_d[index + 13 * ms] = b_density_d * phi3D_d[13];
				f_d[index + 14 * ms] = b_f_d[index + 14 * ms] = b_density_d * phi3D_d[14];
				f_d[index + 15 * ms] = b_f_d[index + 15 * ms] = b_density_d * phi3D_d[15];
				f_d[index + 16 * ms] = b_f_d[index + 16 * ms] = b_density_d * phi3D_d[16];
				f_d[index + 17 * ms] = b_f_d[index + 17 * ms] = b_density_d * phi3D_d[17];
				f_d[index + 18 * ms] = b_f_d[index + 18 * ms] = b_density_d * phi3D_d[18];
			}
			break;
		default:
			break;
		}
	}
}

__host__ void initConstants3D(Arguments *args,
		FLOAT_TYPE maxInletCoordY, FLOAT_TYPE minInletCoordY,
		FLOAT_TYPE maxInletCoordZ, FLOAT_TYPE minInletCoordZ,
		FLOAT_TYPE delta, int m, int n, int h) {
	//CONSTANT LATTICE QUANTITIES d3q19
	//     D3Q19 LATTICE CONFIGURATION

	// 	  Z= 1 LAYER				// 	  Z=0 LAYER					// 	  Z= -1 LAYER
	//                15       		//        8       3       7		//                17
	//                |     		//          \     |     /		//                |
	//                |   			//            \   |   /		    //                |
	//                | 			//              \ | /			//                |
	//       12 - - - 5 - - - 11	//        2 - - - 0 - - - 1		//        14 - - -6 - - - 13
	//  |             | 			//              / | \			//                |
	//  |             |  			//            /   |   \			//                |
	//  |----> x      |     		//          /     |     \		//                |
	//                16   			//        10      4      9		//                18

	int s = m * n * h;
	FLOAT_TYPE w[19] =
	{ 1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1.
			/ 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1.
			/ 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1.
			/ 36. };

	//  			      0   1   2   3   4   5   6    7    8   9  10  11   12   13   14   15    16   17   18
	int opp3D[19] = { 0 * s, 2 * s, 1 * s, 4 * s, 3 * s, 6 * s, 5 * s, 10 * s, 9
			* s, 8 * s, 7 * s, 14 * s,
			//    					12   13   14   15    16   17   18
			13 * s, 12 * s, 11 * s, 18 * s, 17 * s, 16 * s, 15 * s };

	//					   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
	int cx3D[19] = { 0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0,  0 };
	int cy3D[19] = { 0, 0,  0, 1, -1, 0,  0, 1,  1, -1, -1, 0,  0,  0,  0, 1, -1,  1, -1 };
	int cz3D[19] = { 0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1,  1, -1, -1, 1,  1, -1, -1 };
	//   				   0   1   2   3   4   5   6     7      8     9    10    11     12     13
	int c3D[19] = { 0, -1, 1, -1 * n, n, -m * n, +m * n, -1 * n - 1, -1 * n + 1,
			n - 1, n + 1, -m * n - 1, -m * n + 1, +m * n - 1,
			//   	    			14      15       16    17     18
			m * n + 1, -m * n - n, -m * n + n, m * n - n, m * n + n };

	// Calculate collision freq
	FLOAT_TYPE omega = 1.0 / (3. * args->viscosity + 0.5);
	FLOAT_TYPE omegaA = 8 * (2 - omega) / (8 - omega);

	cudaMemcpyToSymbol(cx3D_d, cx3D, 19 * sizeof(int));
	cudaMemcpyToSymbol(cy3D_d, cy3D, 19 * sizeof(int));
	cudaMemcpyToSymbol(cz3D_d, cz3D, 19 * sizeof(int));
	cudaMemcpyToSymbol(c3D_d, c3D, 19 * sizeof(int));
	cudaMemcpyToSymbol(opp3D_d, opp3D, 19 * sizeof(int));
	cudaMemcpyToSymbol(w3D_d, w, 19 * sizeof(FLOAT_TYPE));

	cudaMemcpyToSymbol(g_d, &args->g, sizeof(FLOAT_TYPE));

	cudaMemcpyToSymbol(outletProfile_d, &args->outletProfile,
			sizeof(OutletProfile));
	cudaMemcpyToSymbol(boundaryType_d, &args->boundaryType,
			sizeof(BoundaryType));
	cudaMemcpyToSymbol(dlBoundaryId_d, &args->boundaryId, sizeof(int));

	cudaMemcpyToSymbol(depth_d, &m, sizeof(int));
	cudaMemcpyToSymbol(length_d, &n, sizeof(int));
	cudaMemcpyToSymbol(height_d, &h, sizeof(int));
	cudaMemcpyToSymbol(omega_d, &omega, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(omegaA_d, &omegaA, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(delta_d, &delta, sizeof(FLOAT_TYPE));

	cudaMemcpyToSymbol(inletProfile_d, &args->inletProfile,
			sizeof(InletProfile));
	cudaMemcpyToSymbol(rhoIn_d, &args->rho, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(uIn_d, &args->u, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(vIn_d, &args->v, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(wIn_d, &args->w, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(minInletCoordY_d, &minInletCoordY, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(maxInletCoordZ_d, &maxInletCoordZ, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(minInletCoordY_d, &minInletCoordY, sizeof(FLOAT_TYPE));
	cudaMemcpyToSymbol(maxInletCoordZ_d, &maxInletCoordZ, sizeof(FLOAT_TYPE));

	// Initialize variables for MRT Collision model, if used
	if (args->collisionModel == MRT) {
		FLOAT_TYPE *velMomMap3D = createHostArrayFlt(361);
		FLOAT_TYPE *momCollMtx3D = createHostArrayFlt(361);

		MRTInitializer3D(velMomMap3D, momCollMtx3D, omega);

		cudaMemcpyToSymbol(velMomMap3D_d, velMomMap3D,
				361 * sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(momCollMtx3D_d, momCollMtx3D,
				361 * sizeof(FLOAT_TYPE));
	}

	if(args->multiPhase){

		cudaMemcpyToSymbol(control_param_d, &args->control_param, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(beta_d, &args->beta, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(A_d, &args->A, sizeof(FLOAT_TYPE));

		FLOAT_TYPE w_pert3D[19] = {-2.0 / 9.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
				1.0 / 54.0, 1.0 / 54.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0,
				1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0, 1.0 / 27.0};
		cudaMemcpyToSymbol(w_pert3D_d, w_pert3D, 19 * sizeof(FLOAT_TYPE));

		FLOAT_TYPE phi3D[19] = {0., 1. / 12., 1. / 12., 1. / 12., 1. / 12.,
				1. / 12., 1. / 12., 1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24.,
				1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24.};
		cudaMemcpyToSymbol(phi3D_d, phi3D, 19 * sizeof(FLOAT_TYPE));

		FLOAT_TYPE teta3D[19] = {1., -1. / 12., -1. / 12., -1. / 12., -1. / 12.,
				-1. / 12., -1. / 12., -1. / 24., -1. / 24., -1. / 24., -1. / 24., -1. / 24.,
				-1. / 24., -1. / 24., -1. / 24., -1. / 24., -1. / 24., -1. / 24., -1. / 24.};
		cudaMemcpyToSymbol(teta3D_d, teta3D, 19 * sizeof(FLOAT_TYPE));

		FLOAT_TYPE chi3D[19] = {-5.0 / 2.0, -1. / 6., -1. / 6., -1. / 6., -1. / 6.,
				-1. / 6., -1. / 6., 1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24.,
				1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24., 1. / 24.};

//		FLOAT_TYPE chi3D[19] = {0.0};
		cudaMemcpyToSymbol(chi3D_d, chi3D, 19 * sizeof(FLOAT_TYPE));

		FLOAT_TYPE psi3D[19] = {0., 1. / 4., 1. / 4., 1. / 4., 1. / 4.,
				1. / 4., 1. / 4., 1. / 8., 1. / 8., 1. / 8., 1. / 8., 1. / 8.,
				1. / 8., 1. / 8., 1. / 8., 1. / 8., 1. / 8., 1. / 8., 1. / 8.};
//		FLOAT_TYPE psi3D[19] = {0.0};
		cudaMemcpyToSymbol(psi3D_d, psi3D, 19 * sizeof(FLOAT_TYPE));

		FLOAT_TYPE cg_w3D[19] = {0., 1. / 6., 1. / 6., 1. / 6., 1. / 6.,
				1. / 6., 1. / 6., 1. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12.,
				1. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12., 1. / 12.};
		cudaMemcpyToSymbol(cg_w3D_d, cg_w3D, 19 * sizeof(FLOAT_TYPE));

		cudaMemcpyToSymbol(g_limit_d, &args->g_limit, sizeof(FLOAT_TYPE));

		FLOAT_TYPE c_norms3D[19];
		for(int i = 0; i < 19;i++){
			c_norms3D[i] = sqrt(cx3D[i] * cx3D[i] + cy3D[i] * cy3D[i] + cz3D[i] * cz3D[i]);
		}
		cudaMemcpyToSymbol(c_norms3D_d, c_norms3D, 19 * sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(r_density_d, &args->r_density, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_density_d, &args->b_density, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(r_alpha_d, &args->r_alpha, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_alpha_d, &args->b_alpha, sizeof(FLOAT_TYPE));
		args->bubble_radius /= n;
		cudaMemcpyToSymbol(bubble_radius_d, &args->bubble_radius, sizeof(FLOAT_TYPE));

		cudaMemcpyToSymbol(r_viscosity_d, &args->r_viscosity, sizeof(FLOAT_TYPE));
		cudaMemcpyToSymbol(b_viscosity_d, &args->b_viscosity, sizeof(FLOAT_TYPE));

		FLOAT_TYPE hocg_w[105];
		hocg_w[0] = 0.0;
		for(int i = 1; i < 7;i++)
			hocg_w[i] = 4.0 / 45.0;
		for(int i = 7; i < 19;i++)
			hocg_w[i] = 1.0 / 21.0;
		for(int i = 19; i < 27;i++)
			hocg_w[i] = 2.0 / 105.0;
		for(int i = 27; i < 33;i++)
			hocg_w[i] = 5.0 / 504.0;
		for(int i = 33; i < 57;i++)
			hocg_w[i] = 1.0 / 315.0;
		for(int i = 57; i < 81;i++)
			hocg_w[i] = 1.0 / 630.0;
		for(int i = 81; i < 105;i++)
			hocg_w[i] = 1.0 / 5040.0;


		cudaMemcpyToSymbol(hocg_w3D_d, hocg_w, 105 * sizeof(FLOAT_TYPE));
		int hocg_cx[105] = {0,1,0,0,-1,0,0,1,1,0,-1,1,-1,1,0,0,-1,-1,0,1,-1,1,1,-1,-1,1,
				-1,2,0,0,-2,0,0,2,1,2,1,0,0,-2,2,-1,1,-2,-1,-2,2,-2,-1,1,-1,0,0,0,0,0,0,
				2,1,1,-2,2,2,-2,-2,2,-2,-1,1,1,-1,-1,1,-1,-1,1,1,-1,-1,1,-1,2,2,1,-2,2,2,
				-2,-2,2,-2,-2,2,2,-2,-2,2,-2,-1,1,1,-1,-1,1,-1};
		cudaMemcpyToSymbol(hocg_cx3D_d, hocg_cx, 105 * sizeof(int));
		int hocg_cy[105] = {0,0,1,0,0,-1,0,1,0,1,1,-1,0,0,-1,1,-1,0,-1,1,1,-1,1,-1,1,-1,
				-1,0,2,0,0,-2,0,1,2,0,0,2,1,1,-1,2,-2,-1,-2,0,0,0,0,0,0,-2,2,-2,-1,1,-1,
				1,2,1,1,-1,1,-1,1,-1,-1,2,-2,2,-2,2,-2,-2,1,-1,1,-1,1,-1,-1,2,1,2,2,-2,2
				,-2,2,-2,-2,1,-1,1,-1,1,-1,-1,2,-2,2,-2,2,-2,-2};
		cudaMemcpyToSymbol(hocg_cy3D_d, hocg_cy, 105 * sizeof(int));
		int hocg_cz[105] = {0,0,0,1,0,0,-1,0,1,1,0,0,1,-1,1,-1,0,-1,-1,1,1,1,-1,1,-1,-1,
				-1,0,0,2,0,0,-2,0,0,1,2,1,2,0,0,0,0,0,0,1,-1,-1,2,-2,-2,1,-1,-1,2,-2,-2,
				1,1,2,1,1,-1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,2,2,-2,2,-2,-2,-2,1,2,2,1,1,-1,
				1,-1,-1,-1,2,2,-2,2,-2,-2,-2,2,2,-2,2,-2,-2,-2};
		cudaMemcpyToSymbol(hocg_cz3D_d, hocg_cz, 105 * sizeof(int));

		int hoc3D[105];
		int ms = n * m;
		for(int i = 0; i < 105; i++){
			hoc3D[i] = hocg_cz[i] * ms + hocg_cy[i] * n + hocg_cx[i];
		}
		cudaMemcpyToSymbol(hoc3D_d, hoc3D, 105 * sizeof(int));
	}
}

__global__ void gpuInitInletProfile2D(FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d,
		FLOAT_TYPE *y_d, int size) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	FLOAT_TYPE inletLenghth2 = 0.;

	if (idx < size) {
		inletLenghth2 = (maxInletCoordY_d - minInletCoordY_d)
																																																																																												* (maxInletCoordY_d - minInletCoordY_d);
		u0_d[idx] = 4 * 1.5 * uIn_d * (y_d[idx] - minInletCoordY_d)
																																																																																												* (maxInletCoordY_d - y_d[idx]) / inletLenghth2;
		v0_d[idx] = vIn_d;
	}
}

__global__ void gpuInitInletProfile3D(FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d,
		FLOAT_TYPE *w0_d, FLOAT_TYPE *y_d, FLOAT_TYPE *z_d, int size) {
	int blockId = blockIdx.x + blockIdx.y * gridDim.x;
	int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x)
																																																																																											+ threadIdx.x;
	FLOAT_TYPE rad = 0.;
	FLOAT_TYPE Tta = 0.;
	FLOAT_TYPE eta = 0.;
	FLOAT_TYPE U = 0.;	//Non dymensional U
	FLOAT_TYPE a1 = 0.296;
	FLOAT_TYPE c[5] = { 0.0, -0.046, 0.0, 0.0, 0.0 };  //Epsilon 1==square duct

	if (idx < size) {
		rad = sqrtf((powf(y_d[idx], 2)) + (powf(z_d[idx], 2))); //transform to polar coord
		Tta = atanf(z_d[idx] / y_d[idx]);					//radius and angle

		//        Computation of the velocity profile
		eta = rad / maxInletCoordZ_d;
		//		printf("checkComment: z, %f \n", z_d[idx]);				//checkComment
		printf("checkComment: y, %f \n", y_d[idx]);				//checkComment
		printf("checkComment: eta, %f \n", eta);				//checkComment

		//        for ( i=0; i<6; i++){  UNFOLDED ->
		U = 1 - (0.25 * eta * eta) / a1
				+ (c[0] / a1) * powf(eta, (2. * 1.)) * cos(2. * 1. * Tta)
				+ (c[1] / a1) * powf(eta, (2. * 2.)) * cos(2. * 2. * Tta)
				+ (c[2] / a1) * powf(eta, (2. * 3.)) * cos(2. * 3. * Tta)
				+ (c[3] / a1) * powf(eta, (2. * 4.)) * cos(2. * 4. * Tta)
				+ (c[4] / a1) * powf(eta, (2. * 5.)) * cos(2. * 5. * Tta);
		//        }
		//        printf("checkComment: eta, %f \n",maxInletCoordZ_d);//checkComment
		//        %We adjust the profile so not backflow is found
		//        U(U<0)=0;
		U = (U < 0.0) ? 0.0 : U;
		// TODO: steps 1 and 2
		//        1) We adjust the profile so no-slip velocity is imposed
		//        U(:,1)=0; U(:,end)=0; U(1,:)=0; U(end,:)=0;
		//        2) We adjust the profile to the U_B
		//        Umax = uIn_d*(2*maxInletCoordY_d)*(2*maxInletCoordZ_d)/trapz(y,trapz(z,U,2)); trapz is trapezoidal rule in matlab, check http://csweb.cs.wfu.edu/bigiron/LittleFE-CUDA-TrapezoidalRule/build/html/cudaAlg.html
		//        u0_d[idx] = Umax*U;
		//        Rough estimate in square ducts Umax = 2*uIn_d
		u0_d[idx] = 2 * uIn_d * U;
		printf("checkComment: U, %f \n", U);					//checkComment

	}
}

__host__ int initBoundaryConditions2D(int *bcNodeIdX, int *bcNodeIdY,
		FLOAT_TYPE *q, int *bcBoundId, int *fluid,
		FLOAT_TYPE *bcX, FLOAT_TYPE *bcY, FLOAT_TYPE *nodeX, FLOAT_TYPE *nodeY,
		int *latticeId, int *stream, int *bcType, int *bcMask, int *bcIdx,
		int *mask, FLOAT_TYPE delta, int m, int n, int size) {
	int bci; //bound_dim
	int ms = m * n;
	FLOAT_TYPE dx, dy;
	int opp[9] = { 0, 3 * ms, 4 * ms, 1 * ms, 2 * ms, 7 * ms, 8 * ms, 5 * ms, 6
			* ms };
	FLOAT_TYPE qLat[9] = { 0., 1., 1., 1., 1., sqrt(2), sqrt(2), sqrt(2), sqrt(
			2) };

	for (bci = 0; bci < size; ++bci) {
		int ind = bcNodeIdX[bci] + bcNodeIdY[bci] * n; //node_dim
		int dir = latticeId[bci];
		bcIdx[ind] = ind;
		bcMask[ind] |= BC_MASK(bcType[bci], dir); //was BC_ID
		bcMask[ind] |= BC_F(fluid[ind]);

		if (!(bcMask[ind] & BC_B(BC_ALL))) // if position 16 &17 have not at least one 1
		{							// from a first read, they are set here to
			bcMask[ind] |= BC_B(bcType[bci]);  // the bcType
		} else if ((bcMask[ind] & BC_B(BC_ALL)) != BC_B(bcType[bci])) {
			bcMask[ind] &= ~(BC_B(BC_ALL)); //delete previous
			bcMask[ind] |= BC_WALL_B;
			bcMask[ind] |= BC_CORNER;
			// printf("corner %d\n", ind);
		}

		bcMask[ind] |= BOUND_ID(bcBoundId[bci]); //BoundId must be between 0..255
		dx = bcX[bci] - nodeX[ind];
		dy = bcY[bci] - nodeY[ind];
		q[ind * 8 + dir - 1] = sqrt(dx * dx + dy * dy) / (delta * qLat[dir]); //q = |AC|/|AB| A,B nodos C boundary
		//q_d[ind+(dir-1)*ms] = sqrt( dx*dx + dy*dy ) / (*delta_d * qLat_d[dir]);
		stream[ind + opp[dir] - ms] = 0; //no stream on opposite of the (wall) BC
		//printf("ind:%d(%d,%d) str:%d\n", ind, bcNodeIdX[bci], bcNodeIdY[bci], ind + opp[dir] - ms);

		mask[ind] = 1;
	}
	for (bci = 0; bci < size; ++bci) {
		int ind = bcNodeIdX[bci] + bcNodeIdY[bci] * n; //node_dim
		//corners
		if (bcMask[ind] & BC_CORNER) {
			if (bcMask[ind] & BC_E(BC_ALL) && bcMask[ind] & BC_N(BC_ALL)) {
				bcMask[ind] &= ~(BC_NE(BC_ALL));
				bcMask[ind] |= BC_WALL_NE;
			}
			if (bcMask[ind] & BC_W(BC_ALL) && bcMask[ind] & BC_N(BC_ALL)) {
				bcMask[ind] &= ~(BC_NW(BC_ALL));
				bcMask[ind] |= BC_WALL_NW;
			}
			if (bcMask[ind] & BC_W(BC_ALL) && bcMask[ind] & BC_S(BC_ALL)) {
				bcMask[ind] &= ~(BC_SW(BC_ALL));
				bcMask[ind] |= BC_WALL_SW;
			}
			if (bcMask[ind] & BC_E(BC_ALL) && bcMask[ind] & BC_S(BC_ALL)) {
				bcMask[ind] &= ~(BC_SE(BC_ALL));
				bcMask[ind] |= BC_WALL_SE;
			}
		}
	}
	int i, num = 0;
	for (i = 0; i < m * n; ++i) {
		num += mask[i];
	}
	return num;
}

__host__ void collapseBc2D(int *bcIdx, int *bcIdxCollapsed_d, int *bcMask,
		int *bcMaskCollapsed_d,
		FLOAT_TYPE *q, FLOAT_TYPE *qCollapsed_d, int *mask, int m, int n,
		int size) {
	int *bcIdxCollapsed = (int*) malloc(size * sizeof(int));
	int *bcMaskCollapsed = (int*) malloc(size * sizeof(int));
	FLOAT_TYPE *QCollapsed = (FLOAT_TYPE*) malloc(
			size * 8 * sizeof(FLOAT_TYPE));

	int flyId = 0;
	int i, j;
	for (i = 0; i < n * m; ++i) {
		if (mask[i]) {
			bcIdxCollapsed[flyId] = bcIdx[i];
			bcMaskCollapsed[flyId] = bcMask[i];
			for (j = 0; j < 8; ++j) {
				QCollapsed[flyId * 8 + j] = q[i * 8 + j];
			}
			flyId++;
		}
	}

	cudaMemcpy(bcIdxCollapsed_d, bcIdxCollapsed, size * sizeof(int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(bcMaskCollapsed_d, bcMaskCollapsed, size * sizeof(int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(qCollapsed_d, QCollapsed, 8 * size * sizeof(FLOAT_TYPE),
			cudaMemcpyHostToDevice);
	free(bcIdxCollapsed);
	free(bcMaskCollapsed);
	free(QCollapsed);
}

__host__ int initBoundaryConditions3D(int *bcNodeIdX, int *bcNodeIdY,
		int *bcNodeIdZ, FLOAT_TYPE *q, int *bcBoundId, int *nodeType,
		FLOAT_TYPE *bcX, FLOAT_TYPE *bcY, FLOAT_TYPE *bcZ,
		FLOAT_TYPE *nodeX, FLOAT_TYPE *nodeY, FLOAT_TYPE *nodeZ, int *latticeId,
		bool *stream, int *bcType, unsigned long long *bcMask, int *bcIdx,
		int *mask, FLOAT_TYPE delta, int m, int n, int h, int numConns,
		int CurvedBCs) {
	int bci; //bound_dim
	int ms = m * n * h;
	FLOAT_TYPE dx, dy, dz;
	int opp[19] = { 0, 2 * ms, 1 * ms, 4 * ms, 3 * ms, 6 * ms, 5 * ms, 10 * ms,
			9 * ms, 8 * ms, 7 * ms, 14 * ms, 13 * ms, 12 * ms, 11 * ms, 18 * ms,
			17 * ms, 16 * ms, 15 * ms };
	FLOAT_TYPE qLat[19] = { 0., 1., 1., 1., 1., 1., 1., sqrt(2.), sqrt(2.),
			sqrt(2.), sqrt(2.), sqrt(2.), sqrt(2.), sqrt(2.), sqrt(2.), sqrt(
					2.), sqrt(2.), sqrt(2.), sqrt(2.) };

	for (bci = 0; bci < numConns; ++bci) {
		int ind = bcNodeIdX[bci] + bcNodeIdY[bci] * n + bcNodeIdZ[bci] * n * m; //node_dim
		//         printf("node_dim %d %d %d %d %d \n", bcNodeIdX[bci],bcNodeIdY[bci],bcNodeIdZ[bci], latticeId[bci], nodeType[ind]);
		int dir = latticeId[bci];
		bcIdx[ind] = ind;
		bcMask[ind] |= BC3D_MASK((unsigned long long )bcType[bci], dir);
		//        printf("nodetype: %d ind: %d\n",nodeType[ind], ind );
		//        printf("bcMask[ind] |= BC3D_MASK((unsigned long long)bcType[bci], dir); %#016lX\n", (bcMask[ind] |= BC3D_MASK((unsigned long long)bcType[bci], dir)) );
		bcMask[ind] |= (unsigned long long) BC3D_F(
				(unsigned long long )nodeType[ind]);

		if (!(bcMask[ind]
		             & (unsigned long long) BC3D_B((unsigned long long)BC3D_ALL))) {
			bcMask[ind] |= (unsigned long long) BC3D_B(
					(unsigned long long )bcType[bci]);
		} else if ((bcMask[ind] & BC3D_B((unsigned long long)BC3D_ALL))
				!= BC3D_B((unsigned long long )bcType[bci])) {

			bcMask[ind] &= ~(unsigned long long) (BC3D_B(
					(unsigned long long)BC3D_ALL)); //delete previous
			bcMask[ind] |= (unsigned long long) BC3D_WALL_B;
			bcMask[ind] |= (unsigned long long) BC3D_CORNER;
			//             printf("corner %d\n", ind);
		}

		//[ind] |= BOUND3D_ID(bcBoundId[bci]); //BoundId must be between 0..255
		dx = bcX[bci] - nodeX[ind];
		dy = bcY[bci] - nodeY[ind];
		dz = bcZ[bci] - nodeZ[ind];
		if (CurvedBCs == (BoundaryType) CURVED) {

			q[ind * 18 + dir - 1] = sqrt(dx * dx + dy * dy + dz * dz)
																																																																																													/ (delta * qLat[dir]); //q = |AC|/|AB| A,B nodos C boundary
		}
		//        if(ind == 30757)        printf("q[ind*18 + dir-1]: %f\n", q[ind*18 + dir-1]);
		//q_d[ind+(dir-1)*ms] = sqrt( dx*dx + dy*dy ) / (*delta_d * qLat_d[dir]);
		stream[ind + opp[dir] - ms] = 0; //no stream on opposite of the (wall) BC3D
		//printf("ind:%d(%d,%d) str:%d\n", ind, bcNodeIdX[bci], bcNodeIdY[bci], ind + opp[dir] - ms);

		mask[ind] = 1;
	}
	for (bci = 0; bci < numConns; ++bci) {
		int ind = bcNodeIdX[bci] + bcNodeIdY[bci] * n + bcNodeIdZ[bci] * n * m;
		; //node_dim
		if (bcMask[ind] & BC3D_CORNER) {
			for (int dir_i = 0; dir_i < 19; dir_i++) {
				if (bcMask[ind] & BC3D_MASK((unsigned long long)BC3D_ALL, dir_i)) {
					bcMask[ind] &= ~(unsigned long long) BC3D_MASK(
							(unsigned long long)BC3D_ALL, dir_i);
					bcMask[ind] |= BC3D_MASK((unsigned long long)BC3D_WALL,
							dir_i);

				}
			}

		}

	}
	int i, num = 0;
	for (i = 0; i < m * n * h; ++i) {
		num += mask[i];
	}
	return num;

}
__host__ void collapseBc3D(int *bcIdx, int *bcIdxCollapsed_d,
		unsigned long long *bcMask, unsigned long long *bcMaskCollapsed_d,
		FLOAT_TYPE *q, FLOAT_TYPE *qCollapsed_d, int *mask, int m, int n, int h,
		int size, int CurvedBCs) {
	int *bcIdxCollapsed = (int*) malloc(size * sizeof(int));
	unsigned long long *bcMaskCollapsed = (unsigned long long*) malloc(
			size * sizeof(unsigned long long));
	FLOAT_TYPE *QCollapsed;
	if (CurvedBCs == (BoundaryType) CURVED)
		QCollapsed = (FLOAT_TYPE*) malloc(size * 18 * sizeof(FLOAT_TYPE));

	int flyId = 0;
	int i, j;
	for (i = 0; i < n * m * h; ++i) {
		if (mask[i]) {
			bcIdxCollapsed[flyId] = bcIdx[i];
			bcMaskCollapsed[flyId] = bcMask[i];
			if (CurvedBCs == (BoundaryType) CURVED) {
				for (j = 0; j < 18; ++j) {
					QCollapsed[flyId * 18 + j] = q[i * 18 + j];
				}
			}
			flyId++;
		}
	}

	cudaMemcpy(bcIdxCollapsed_d, bcIdxCollapsed, size * sizeof(int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(bcMaskCollapsed_d, bcMaskCollapsed,
			size * sizeof(unsigned long long), cudaMemcpyHostToDevice);
	free(bcIdxCollapsed);
	free(bcMaskCollapsed);
	if (CurvedBCs == (BoundaryType) CURVED) {
		cudaMemcpy(qCollapsed_d, QCollapsed, 18 * size * sizeof(FLOAT_TYPE),
				cudaMemcpyHostToDevice);
		free(QCollapsed);
	}

}
