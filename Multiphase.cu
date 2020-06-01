#include "Multiphase.h"
#include <math.h>
#include <stdio.h>
#include "ArrayUtils.h"

void mp2DColl(int n, int m, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE r_omega, FLOAT_TYPE b_omega, FLOAT_TYPE control_param, FLOAT_TYPE del,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE r_A,  FLOAT_TYPE b_A, FLOAT_TYPE *r_fPert, FLOAT_TYPE *b_fPert,
		FLOAT_TYPE *weight, int *cx, int *cy, FLOAT_TYPE r_nu, FLOAT_TYPE b_nu){

	FLOAT_TYPE cu1, cu2, r_CollPert, b_CollPert;
	FLOAT_TYPE cosin;
	FLOAT_TYPE color_gradient_norm;
	FLOAT_TYPE k_r, k_b, k_k;
	FLOAT_TYPE norm_c;
	FLOAT_TYPE prod_c_g;
	FLOAT_TYPE r_pert, b_pert;
	FLOAT_TYPE r_feq, b_feq;
	FLOAT_TYPE fn05, aux1, mean_nu, omega_eff;
	FLOAT_TYPE cg_w[9] = {0.0, 4.0/12.0, 4.0/12.0, 4.0/12.0, 4.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0};
	int index, index9, temp_index;
	for (int j=0; j < m; j++){
		for (int i=0;i < n; i++){
			index = j*n + i;


			color_gradient[index * 2] = 0;
			color_gradient[index * 2 + 1] = 0;
			for (int k=0; k<9; k++){

				// calculate color gradient - 4th order

				if (i!=0 && j!=0 && i!=(n-1) && j!=(m-1)){ // Interior points - In the boundary it is calculated by "mirroring" the density
					temp_index = (j + cy[k]) * n + i + cx[k];
					color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k] * cg_w[k];
					color_gradient[index * 2 + 1] += (r_rho[temp_index] - b_rho[temp_index]) * cy[k] * cg_w[k];
				}
				else if (j==(m-1) && i!=0 && i!=(n-1)) {// north boundary
					temp_index = (j - abs(cy[k])) * n + i + cx[k];
					color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k] * cg_w[k];
					color_gradient[index * 2 + 1] = 0;
				}
				else if (j==0 && i!=0 && i!=(n-1)){  // south boundary
					temp_index = (j + abs(cy[k])) * n + i + cx[k];
					color_gradient[index * 2] += (r_rho[temp_index] - b_rho[temp_index]) * cx[k] * cg_w[k];
					color_gradient[index * 2 + 1] = 0;
				}
				else if (i==(n-1) && j!=0 && j!=(m-1)){  // east boundary
					temp_index = (j + cy[k]) * n + i - abs(cx[k]);
					color_gradient[index * 2] = 0;
					color_gradient[index * 2 + 1] +=  (r_rho[temp_index] - b_rho[temp_index]) * cy[k] * cg_w[k];
				}
				else if (i==0 && j!=0 && j!=(m-1)){ //  west boundary
					temp_index = (j + cy[k]) * n + i + abs(cx[k]);
					color_gradient[index * 2] = 0;
					color_gradient[index * 2 + 1] +=  (r_rho[temp_index] - b_rho[temp_index]) * cy[k] * cg_w[k];
				}
			}

			aux1 = r_rho[index]/(rho[index]*r_nu) + b_rho[index]/(rho[index]*b_nu);
			mean_nu = 1.0/aux1;
			omega_eff = 1.0/(3.0*mean_nu+0.5);

			// relaxation parameter to choose a proper omega at the interface
			//			if (r_omega != b_omega){
			//				chi=(r_rho[index] - b_rho[index])/rho[index];
			//				if(chi >= -control_param && chi <= control_param){
			//					if (chi > del)
			//						r_omega_temp=r_omega;
			//					else if (chi <= del && chi > 0)
			//						r_omega_temp=a1 + a2 * chi + a3 * chi * chi;
			//					else if (chi <= 0 && chi >= -del)
			//						r_omega_temp=a1 + a4 * chi + a5 * chi * chi;
			//					else if (chi < -del)
			//						r_omega_temp=b_omega;
			//
			//					b_omega_temp = r_omega_temp;
			//				}
			//				else{
			//					r_omega_temp = r_omega;
			//					b_omega_temp = b_omega;
			//				}
			//			}
			//			else
			//			{
			//				r_omega_temp=r_omega;
			//				b_omega_temp=r_omega_temp;
			//			}
			//			printf("romega "FLOAT_FORMAT" bomega"FLOAT_FORMAT" \n", r_omega_temp, b_omega_temp);
			cu1 = u[index]*u[index] + v[index]*v[index];

			// invariable quantities
			color_gradient_norm = sqrt(pow(color_gradient[index * 2],2) + pow(color_gradient[index * 2 + 1],2));
			k_r=r_rho[index]/rho[index];
			k_b=b_rho[index]/rho[index];
			k_k= beta * r_rho[index] * b_rho[index]/(pow(rho[index],2));
			for (int k=0;k<9;k++){
				if (color_gradient_norm > g_limit){
					prod_c_g=cx[k]*color_gradient[index * 2]+cy[k]*color_gradient[index * 2 + 1];
					if (k!=0){
						norm_c= sqrt(pow(cx[k],2)+pow(cy[k],2));
						cosin= prod_c_g / (color_gradient_norm*norm_c);
					}
					else
						cosin=0.0;
					// calculate perturbation terms

					r_pert=0.5*r_A*color_gradient_norm*(weight[k]*pow((prod_c_g/color_gradient_norm),2)-w_pert[k]);
					b_pert=0.5*b_A*color_gradient_norm*(weight[k]*pow((prod_c_g/color_gradient_norm),2)-w_pert[k]);

				}
				else{
					// ther perturbation terms are null
					r_pert=0.0;
					b_pert=0.0;
				}


				cu2 = u[index]*cx[k] + v[index]*cy[k];
				// calculate equilibrium distribution function
				r_feq = r_rho[index] * (r_phi[k] + weight[k] * (3 * cu2 + 4.5 * cu2 * cu2 - 1.5 * cu1));
				b_feq = b_rho[index] * (b_phi[k] + weight[k] * (3 * cu2 + 4.5 * cu2 * cu2 - 1.5 * cu1));

				index9 = i + j * n + k * m * n;
				// calculate updated distribution function
				r_CollPert = omega_eff*r_feq + (1-omega_eff)*r_f[index9]+r_pert;
				b_CollPert = omega_eff*b_feq + (1-omega_eff)*b_f[index9]+b_pert;

				fn05 = r_CollPert + b_CollPert;
				//				// perform recolor step
				r_fPert[index9]=k_r*fn05+k_k*cosin*(r_rho[index]*r_phi[k]+b_rho[index]*b_phi[k]);
				b_fPert[index9]=k_b*fn05-k_k*cosin*(r_rho[index]*r_phi[k]+b_rho[index]*b_phi[k]);
			}
		}
	}
}

void mp3DColl(int n, int m, int h, FLOAT_TYPE *rho, FLOAT_TYPE *u,
		FLOAT_TYPE *v, FLOAT_TYPE *w, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE *w_pert, FLOAT_TYPE *color_gradient,
		FLOAT_TYPE beta, FLOAT_TYPE g_limit,  FLOAT_TYPE A, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl,
		FLOAT_TYPE *weight, int *cx, int *cy, int *cz, FLOAT_TYPE *f, FLOAT_TYPE r_nu, FLOAT_TYPE b_nu, FLOAT_TYPE r_alpha,
		FLOAT_TYPE b_alpha, FLOAT_TYPE *chi, FLOAT_TYPE *phi, FLOAT_TYPE *psi, FLOAT_TYPE *teta, FLOAT_TYPE *cg_w){

	FLOAT_TYPE cu1, cu2, f_CollPert;
	FLOAT_TYPE cosin;
	FLOAT_TYPE color_gradient_norm;
	FLOAT_TYPE k_r, k_b, k_k;
	FLOAT_TYPE norm_c;
	FLOAT_TYPE prod_c_g;
	FLOAT_TYPE pert;
	FLOAT_TYPE f_eq;
	int index, index9, temp_index;
	FLOAT_TYPE grad_rho_x, grad_rho_y, grad_rho_z;
	FLOAT_TYPE G[9];
	FLOAT_TYPE H[9];
	FLOAT_TYPE prod_u_grad_rho, aux1, mean_nu, omega_eff, mean_alpha, TC;

	for(int k = 0; k < h; k++){
		for (int j=0; j < m; j++){
			for (int i=0;i < n; i++){
				index = k*n*m + j*n + i;

				color_gradient[index * 3] = 0;
				color_gradient[index * 3 + 1] = 0;
				color_gradient[index * 3 + 2] = 0;
				grad_rho_x = 0.0;
				grad_rho_y = 0.0;
				grad_rho_z = 0.0;
				for (int dir=0; dir < 19; dir++){

					// calculate color gradient - 4th order

					if (i!=0 && j!=0 && k != 0 && i!=(n-1) && j!=(m-1) && k != (h-1)){ // Interior points - In the boundary it is calculated by "mirroring" the density
						temp_index = (k + cz[dir]) * n * m + (j + cy[dir]) * n + i + cx[dir];
						grad_rho_x += rho[temp_index] * cx[dir] * cg_w[dir];
						grad_rho_y += rho[temp_index] * cy[dir] * cg_w[dir];
						grad_rho_z += rho[temp_index] * cz[dir] * cg_w[dir];
						color_gradient[index * 3] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cx[dir] * cg_w[dir];
						color_gradient[index * 3 + 1] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cy[dir] * cg_w[dir];
						color_gradient[index * 3 + 2] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cz[dir] * cg_w[dir];
					}
					else if (j==(m-1) && i!=0 && i!=(n-1)) {// north boundary
						temp_index = (k + cz[dir]) * n * m + (j - abs(cy[dir])) * n + i + cx[dir];
						grad_rho_x += rho[temp_index] * cx[dir] * cg_w[dir];
						grad_rho_y = 0;
						grad_rho_z += rho[temp_index] * cz[dir] * cg_w[dir];
						color_gradient[index * 3] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cx[dir] * cg_w[dir];
						color_gradient[index * 3 + 1] = 0;
						color_gradient[index * 3 + 2] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cz[dir] * cg_w[dir];
					}
					else if (j==0 && i!=0 && i!=(n-1)){  // south boundary
						temp_index = (k + cz[dir]) * n * m + (j + abs(cy[dir])) * n + i + cx[dir];
						grad_rho_x += rho[temp_index] * cx[dir] * cg_w[dir];
						grad_rho_y = 0;
						grad_rho_z += rho[temp_index] * cz[dir] * cg_w[dir];
						color_gradient[index * 3] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cx[dir] * cg_w[dir];
						color_gradient[index * 3 + 1] = 0;
						color_gradient[index * 3 + 2] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cz[dir] * cg_w[dir];
					}
					else if (i==(n-1) && j!=0 && j!=(m-1)){  // east boundary
						temp_index = (k + cz[dir]) * n * m + (j + cy[dir]) * n + i - abs(cx[dir]);
						grad_rho_x = 0;
						grad_rho_y += rho[temp_index] * cy[dir] * cg_w[dir];
						grad_rho_z += rho[temp_index] * cz[dir] * cg_w[dir];
						color_gradient[index * 3] = 0;
						color_gradient[index * 3 + 1] +=  ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cy[dir] * cg_w[dir];
						color_gradient[index * 3 + 2] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cz[dir] * cg_w[dir];
					}
					else if (i==0 && j!=0 && j!=(m-1)){ //  west boundary
						temp_index = (k + cz[dir]) * n * m + (j + cy[dir]) * n + i + abs(cx[dir]);
						grad_rho_x = 0;
						grad_rho_y += rho[temp_index] * cy[dir] * cg_w[dir];
						grad_rho_z += rho[temp_index] * cz[dir] * cg_w[dir];
						color_gradient[index * 3] = 0;
						color_gradient[index * 3 + 1] +=  ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cy[dir] * cg_w[dir];
						color_gradient[index * 3 + 2] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cz[dir] * cg_w[dir];
					}
					else if (k==(h-1) && j!=0 && j!=(m-1)){ //  front boundary
						temp_index = (k - abs(cz[dir])) * n * m + (j + cy[dir]) * n + i + abs(cx[dir]);
						grad_rho_x += rho[temp_index] * cx[dir] * cg_w[dir];
						grad_rho_y += rho[temp_index] * cy[dir] * cg_w[dir];
						grad_rho_z = 0;
						color_gradient[index * 3] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cx[dir] * cg_w[dir];
						color_gradient[index * 3 + 1] +=  ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cy[dir] * cg_w[dir];
						color_gradient[index * 3 + 2] = 0;
					}
					else if (k==0 && j!=0 && j!=(m-1)){ //  back boundary
						temp_index = (k + abs(cz[dir])) * n * m + (j + cy[dir]) * n + i + abs(cx[dir]);
						grad_rho_x += rho[temp_index] * cx[dir] * cg_w[dir];
						grad_rho_y += rho[temp_index] * cy[dir] * cg_w[dir];
						grad_rho_z = 0;
						color_gradient[index * 3] += ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cx[dir] * cg_w[dir];
						color_gradient[index * 3 + 1] +=  ((r_rho[temp_index] - b_rho[temp_index]) / rho[temp_index]) * cy[dir] * cg_w[dir];
						color_gradient[index * 3 + 2] = 0;
					}
				}



				G[0] = 2.0 * u[index] * grad_rho_x;
				G[1] = u[index]*grad_rho_y + v[index]*grad_rho_x;
				G[2] = u[index]*grad_rho_z + w[index]*grad_rho_x;
				G[3] = G[1];
				G[4] = 2.0*v[index]*grad_rho_y;
				G[5] = v[index]*grad_rho_z + w[index]*grad_rho_y;
				G[6] = G[2];
				G[7] = G[5];
				G[8] = 2.0*w[index]*grad_rho_z;

				prod_u_grad_rho = u[index]*grad_rho_x + v[index]*grad_rho_y + w[index]*grad_rho_z;

				cu1 = u[index]*u[index] + v[index]*v[index] + w[index] * w[index];

				// invariable quantities
				color_gradient_norm = sqrt(pow(color_gradient[index * 3],2) + pow(color_gradient[index * 3 + 1],2) + pow(color_gradient[index * 3 + 2],2));
				k_r=r_rho[index]/rho[index];
				k_b=b_rho[index]/rho[index];
				k_k= beta * r_rho[index] * b_rho[index]/rho[index];

				aux1 = r_rho[index]/(rho[index]*r_nu) + b_rho[index]/(rho[index]*b_nu);
				mean_nu = 1.0/aux1;

				omega_eff = 1.0/(3.0*mean_nu+0.5);

				mean_alpha = r_alpha*r_rho[index]/rho[index] + b_alpha*b_rho[index]/rho[index];

				for (int dir=0;dir<19;dir++){

					if (color_gradient_norm > g_limit){
						prod_c_g=cx[dir]*color_gradient[index * 3]+cy[dir]*color_gradient[index * 3 + 1] + cz[dir] * color_gradient[index * 3 + 2];
						if (dir!=0){
							norm_c= sqrt(pow(cx[dir],2)+pow(cy[dir],2) + pow(cz[dir],2));
							cosin= prod_c_g / (color_gradient_norm*norm_c);
						}
						else
							cosin=0.0;
						// calculate perturbation terms

						pert=0.5*A*color_gradient_norm*(weight[dir]*pow((prod_c_g/color_gradient_norm),2)-w_pert[dir]);

					}
					else{
						// ther perturbation terms are null
						pert = 0.0;
					}
					// Auxiliar tensor: diadic product of the speed velcity:
					//[cx,cy,cx]*[cx cy cz]
					H[0] = cx[dir]*cx[dir];
					H[1] = cx[dir]*cy[dir];
					H[2] = cx[dir]*cz[dir];
					H[3] = H[1];
					H[4] = cy[dir]*cy[dir];
					H[5] = cy[dir]*cz[dir];
					H[6] = H[2];
					H[7] = H[5];
					H[8] = cz[dir]*cz[dir];

					//Tensor contraction
					TC = 0;
					for(int l = 0; l < 9; l++){
						TC += G[l] * H[l];
					}



					cu2 = u[index]*cx[dir] + v[index]*cy[dir] + w[index]*cz[dir];
					// calculate equilibrium distribution function
					f_eq = mean_alpha*(chi[dir]*prod_u_grad_rho + psi[dir]*TC) + rho[index] * ( phi[dir] + teta[dir]*mean_alpha + weight[dir] * (3*cu2+4.5*cu2*cu2-1.5*cu1));

					index9 = index + dir * n * m * h;
					// calculate updated distribution function
					f_CollPert = omega_eff*f_eq + (1-omega_eff)*f[index9] + pert;

					r_fColl[index9] = k_r*f_CollPert + k_k*cosin*(phi[dir]+teta[dir]*mean_alpha);
					b_fColl[index9] = k_b*f_CollPert + k_k*cosin*(phi[dir]+teta[dir]*mean_alpha);
				}
			}
		}
	}
}

void createBubble(FLOAT_TYPE *x, FLOAT_TYPE *y,int n, int m, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *r_phi, FLOAT_TYPE *b_phi, FLOAT_TYPE *rho, int test_case) {
	int i, j, k;
	int index, index2;
	switch (test_case) {
	case 1:
		for (j=0; j < m; j++){
			for(i = 0; i < n; i++){
				index = j * n + i;

				if( sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5),2)) <= radius ){
					r_rho[index] = r_density;
					for (k=0; k < 9; k++){
						// initialise distribution function with small, non-zero values
						index2 = i + j * n + k * m * n;
						r_f[index2] = r_rho[index] * r_phi[k];
					}
				}
				else {
					b_rho[index]=b_density;
					for (k=0; k < 9; k++){
						// initialise distribution function with small, non-zero values
						index2 = i + j * n + k * m * n;
						b_f[index2]   = b_rho[index]*b_phi[k];
					}
				}
				// initialise density
				rho[index] = r_rho[index]+b_rho[index];
			}
		}
		break;
	case 2:
		for (j=0; j < m; j++){
			for(i = 0; i < n; i++){
				index = j * n + i;
				if( y[index] > 0.5 ){
					//r_rho[index] = r_density;
					for (k=0; k < 9; k++){
						r_rho[index] = r_density;
						// initialise distribution function with small, non-zero values
						index2 = i + j * n + k * m * n;
						r_f[index2] = r_density * r_phi[k];
					}
				}
				else {
					//b_rho[index]=b_density;
					for (k=0; k < 9; k++){
						b_rho[index] = b_density;
						// initialise distribution function with small, non-zero values
						index2 = i + j * n + k * m * n;
						b_f[index2]   = b_density*b_phi[k];
					}
				}
				// initialise density
				rho[index] = r_rho[index]+b_rho[index];
			}
		}
		break;
	case 3:
		for (j=0; j < m; j++){
			for(i = 0; i < n; i++){
				index = j * n + i;

				if( sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5 + radius),2)) <= radius || sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5 - radius),2)) <= radius){
					r_rho[index] = r_density;
					for (k=0; k < 9; k++){
						// initialise distribution function with small, non-zero values
						index2 = i + j * n + k * m * n;
						r_f[index2] = r_rho[index] * r_phi[k];
					}
				}
				else {
					b_rho[index]=b_density;
					for (k=0; k < 9; k++){
						// initialise distribution function with small, non-zero values
						index2 = i + j * n + k * m * n;
						b_f[index2]   = b_rho[index]*b_phi[k];
					}
				}
				// initialise density
				rho[index] = r_rho[index]+b_rho[index];
			}
		}
		break;
	default:
		break;
	}

}

void createBubble3D(FLOAT_TYPE *x, FLOAT_TYPE *y, FLOAT_TYPE *z, int n, int m, int h, FLOAT_TYPE radius, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho,
		FLOAT_TYPE r_density, FLOAT_TYPE b_density, FLOAT_TYPE *phi, FLOAT_TYPE *rho, FLOAT_TYPE *f) {
	int i, j, k, dir;
	int index, index2;
	for(k = 0; k < h; k++){
		for (j=0; j < m; j++){
			for(i = 0; i < n; i++){
				index = k * m * n + j * n + i;

				if( sqrt( pow((x[index]-0.5), 2) + pow((y[index]-0.5),2) + pow((z[index]-0.5),2)) <= radius ){
					r_rho[index] = r_density;
					for (dir=0; dir < 19; dir++){
						// initialise distribution function with small, non-zero values
						index2 = index + dir * h*m*n;
						r_f[index2] = r_rho[index] * phi[dir];
					}
				}
				else {
					b_rho[index]=b_density;
					for (dir=0; dir < 19; dir++){
						// initialise distribution function with small, non-zero values
						index2 = index + dir * h*m*n;
						b_f[index2]   = b_rho[index]*phi[dir];
					}
				}
				for(dir = 0; dir < 19; dir++){
					index2 = index + dir * h*m*n;
					f[index2] = b_f[index2] + r_f[index2];
				}
				// initialise density
				rho[index] = r_rho[index]+b_rho[index];
			}
		}
	}
}

void updateMacroMP(int n, int m, FLOAT_TYPE *u, FLOAT_TYPE *v,FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *rho, FLOAT_TYPE control_param,
		FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE st_predicted, int test_case){

	int index_aux1=0;
	int index_aux2=0;
	FLOAT_TYPE p_in=0.0;
	FLOAT_TYPE p_out=0.0;
	FLOAT_TYPE u_cum, v_cum;
	FLOAT_TYPE r_sum, b_sum;
	FLOAT_TYPE chi;
	int index, index9;
	int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	FLOAT_TYPE st_laplace;
	int start_index, end_index;
	switch(test_case){
	case 1:
		start_index = 0;
		end_index = m;
		break;
	case 2:
		start_index = 1;
		end_index = m-1;
		break;
	default:
		start_index = 0;
		end_index = m;
		break;
	}
	for (int j=start_index; j < end_index;j++){
		for (int i=0; i < n; i++){


			// densities
			index = j * n + i;
			r_sum = 0.0;
			b_sum = 0.0;
			for(int k = 0; k < 9; k++){
				index9 = i + j * n + k * m * n;
				r_sum += r_f[index9];
				b_sum += b_f[index9];
			}
			r_rho[index] = r_sum;
			b_rho[index]= b_sum;
			rho[index] = r_rho[index]+b_rho[index];

			// p_in and p_out for the surface tension
			chi=(r_rho[index]-b_rho[index])/rho[index];
			if (chi >= control_param){
				index_aux1++;
				p_in += r_rho[index];
			}
			else if (chi <= -control_param){
				index_aux2++;
				p_out+=b_rho[index];
			}

			// auxiliar variables
			u_cum=0.0;
			v_cum=0.0;
			// velocities
			for (int k=0; k < 9; k++){
				index9 = i + j * n + k * m * n;
				u_cum += (r_f[index9]+b_f[index9])*cx[k];
				v_cum += (r_f[index9]+b_f[index9])*cy[k];
			}
			u[index]   = u_cum/rho[index];
			v[index]  = v_cum/rho[index];

		}
	}

	// Calculate surface tension
	p_in=(3.0/5.0)*(1.0-r_alpha)*p_in/index_aux1;      // pressure average inside the bubble
	p_out=(3.0/5.0)*(1.0-b_alpha)*p_out/index_aux2;   // pressure average outside the bubble
	st_laplace=bubble_radius*(p_in-p_out);

	st_error[iteration]=abs(st_predicted-st_laplace)/(st_predicted)*100.0;
}

void updateMacroMP3D(int n, int m, int h, FLOAT_TYPE *u, FLOAT_TYPE *v, FLOAT_TYPE *w,FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *rho, FLOAT_TYPE control_param,
		FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE st_predicted, int *cx, int *cy, int *cz, FLOAT_TYPE *f){

	int index_aux1=0;
	int index_aux2=0;
	FLOAT_TYPE p_in=0.0;
	FLOAT_TYPE p_out=0.0;
	FLOAT_TYPE u_cum, v_cum, w_cum;
	FLOAT_TYPE r_sum, b_sum;
	FLOAT_TYPE chi;
	int index, index9;
	FLOAT_TYPE st_laplace;
	// Density and Velocity
	for(int k = 1; k < h - 1; k++){
		for (int j=1; j < m - 1;j++){
			for (int i=1; i < n - 1; i++){
				// auxiliar variables
				u_cum=0.0;
				v_cum=0.0;
				w_cum=0.0;
				// densities
				index = k * m * n + j * n + i;
				r_sum = 0.0;
				b_sum = 0.0;
				for(int dir = 0; dir < 19; dir++){
					r_sum += r_f[index + dir * m * n * h];
					b_sum += b_f[index + dir * m * n * h];
					f[index + dir * m * n * h] = r_f[index + dir * m * n * h] + b_f[index + dir * m * n * h];
				}
				r_rho[index] = r_sum;
				b_rho[index]= b_sum;
				rho[index] = r_rho[index]+b_rho[index];

				// p_in and p_out for the surface tension
				chi=(r_rho[index]-b_rho[index])/rho[index];
				if (chi >= control_param){
					index_aux1++;
					p_in += r_rho[index];
				}
				else if (chi <= -control_param){
					index_aux2++;
					p_out+=b_rho[index];
				}

				// velocities
				for (int dir=0; dir < 19; dir++){
					index9 = index + dir * m * n * h;
					u_cum += (r_f[index9]+b_f[index9])*cx[dir];
					v_cum += (r_f[index9]+b_f[index9])*cy[dir];
					w_cum += (r_f[index9]+b_f[index9])*cz[dir];
				}
				u[index]   = u_cum/rho[index];
				v[index]  = v_cum/rho[index];
				w[index]  = w_cum/rho[index];

			}
		}
	}

	// Calculate surface tension
	p_in=(3.0/5.0)*(1.0-r_alpha)*p_in/index_aux1;      // pressure average inside the bubble
	p_out=(3.0/5.0)*(1.0-b_alpha)*p_out/index_aux2;   // pressure average outside the bubble
	st_laplace=bubble_radius*(p_in-p_out);

	st_error[iteration]=abs(st_predicted-st_laplace)/(st_predicted)*100.0;
}

void peridicBoundaries(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *u, FLOAT_TYPE *v, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE *rho, int test_case){

	int index_end, index_start;
	int jn = m-1;
	int js = 0;
	int ie = n-1;
	int iw = 0;
	FLOAT_TYPE r_temp, b_temp, u_temp, v_temp;


	if(test_case == 2){
		for (int i=1; i < n-1; i++){
			// north boundary
			index_end = jn * n + i;
			index_start = js * n + i;

			u_temp = u[index_end];
			v_temp = v[index_end];
			r_temp = (1. / (1. + v_temp)) * (r_f[index_end] + r_f[index_end + 1 * m * n] + r_f[index_end + 3 * m * n]
			                                                                                   + 2 * (r_f[index_end + 2 * m * n] + r_f[index_end + 6 * m * n] + r_f[index_end + 5 * m * n]) );
			b_temp = (1. / (1. + v_temp)) * (b_f[index_end] + b_f[index_end + 1 * m * n] + b_f[index_end + 3 * m * n]
			                                                                                   + 2 * (b_f[index_end + 2 * m * n] + b_f[index_end + 6 * m * n] + b_f[index_end + 5 * m * n]) );

			r_f[index_end + 4 * m * n] = r_f[index_end + 2 * m * n] - (2. / 3.) * r_temp*v_temp;
			r_f[index_end + 7 * m * n] = r_f[index_end + 5 * m * n] + 0.5 * (r_f[index_end + 1 * m * n] - r_f[index_end + 3 * m * n]) - (1. / 6.) * (r_temp*v_temp) - 0.5 * (r_temp*u_temp);
			r_f[index_end + 8 * m * n] = r_f[index_end + 6 * m * n] + 0.5 * (r_f[index_end + 3 * m * n] - r_f[index_end + 1 * m * n]) - (1. / 6.) * (r_temp*v_temp) + 0.5 * (r_temp*u_temp);

			b_f[index_end + 4 * m * n] = b_f[index_end + 2 * m * n] - (2. / 3.) * b_temp*v_temp;
			b_f[index_end + 7 * m * n] = b_f[index_end + 5 * m * n] + 0.5 * (b_f[index_end + 1 * m * n] - b_f[index_end + 3 * m * n]) - (1. / 6.) * (b_temp*v_temp) - 0.5 * (b_temp*u_temp);
			b_f[index_end + 8 * m * n] = b_f[index_end + 6 * m * n] + 0.5 * (b_f[index_end + 3 * m * n] - b_f[index_end + 1 * m * n]) - (1. / 6.) * (b_temp*v_temp) + 0.5 * (b_temp*u_temp);

			r_rho[index_end] = r_temp;
			b_rho[index_end] = b_temp;
			rho[index_end] = r_temp + b_temp;

			//south boundary
			u_temp = u[index_start];
			v_temp = v[index_start];

			r_temp = (1. / (1. - v_temp)) * (r_f[index_start] + r_f[index_start + 1 * m * n] + r_f[index_start + 3 * m * n]
			                                                                                       + 2 * (r_f[index_start + 4 * m * n] + r_f[index_start + 7 * m * n] + r_f[index_start + 8 * m * n]));
			b_temp = (1. / (1. - v_temp)) * (b_f[index_start] + b_f[index_start + 1 * m * n] + b_f[index_start + 3 * m * n]
			                                                                                       + 2 * (b_f[index_start + 4 * m * n] + b_f[index_start + 7 * m * n] + b_f[index_start + 8 * m * n]));

			r_f[index_start + 2 * m * n] = r_f[index_start + 4 * m * n] + (2. / 3.) * r_temp*v_temp;
			r_f[index_start + 5 * m * n] = r_f[index_start + 7 * m * n] - 0.5 * (r_f[index_start + 1 * m * n] - r_f[index_start + 3 * m * n]) + (1. / 6.) * (r_temp*v_temp) + 0.5 * (r_temp*u_temp);
			r_f[index_start + 6 * m * n] = r_f[index_start + 8 * m * n] + 0.5 * (r_f[index_start + 1 * m * n] - r_f[index_start + 3 * m * n]) + (1. / 6.) * (r_temp*v_temp) - 0.5 * (r_temp*u_temp);

			b_f[index_start + 2 * m * n] = b_f[index_start + 4 * m * n] + (2. / 3.) * b_temp*v_temp;
			b_f[index_start + 5 * m * n] = b_f[index_start + 7 * m * n] - 0.5 * (b_f[index_start + 1 * m * n] - b_f[index_start + 3 * m * n]) + (1. / 6.) * (b_temp*v_temp) + 0.5 * (b_temp*u_temp);
			b_f[index_start + 6 * m * n] = b_f[index_start + 8 * m * n] + 0.5 * (b_f[index_start + 1 * m * n] - b_f[index_start + 3 * m * n]) + (1. / 6.) * (b_temp*v_temp) - 0.5 * (b_temp*u_temp);

			r_rho[index_start] = r_temp;
			b_rho[index_start] = b_temp;
			rho[index_start] = r_temp + b_temp;
		}
	}
	else{
		for (int i=1; i < n-1; i++){
			// north boundary
			index_end = jn * n + i;
			index_start = js * n + i;

			r_f[index_end + 4 * m * n] = r_f[index_start + 4 * m * n];
			r_f[index_end + 7 * m * n] = r_f[index_start + 7 * m * n];
			r_f[index_end + 8 * m * n] = r_f[index_start + 8 * m * n];

			b_f[index_end + 4 * m * n] = b_f[index_start + 4 * m * n];
			b_f[index_end + 7 * m * n] = b_f[index_start + 7 * m * n];
			b_f[index_end + 8 * m * n] = b_f[index_start + 8 * m * n];

			//south boundary
			r_f[index_start + 2 * m * n] = r_f[index_end + 2 * m * n];
			r_f[index_start + 5 * m * n] = r_f[index_end + 5 * m * n];
			r_f[index_start + 6 * m * n] = r_f[index_end + 6 * m * n];

			b_f[index_start + 2 * m * n] = b_f[index_end + 2 * m * n];
			b_f[index_start + 5 * m * n] = b_f[index_end + 5 * m * n];
			b_f[index_start + 6 * m * n] = b_f[index_end + 6 * m * n];
		}
	}

	for (int j=1; j < m-1; j++){
		// east boundary
		index_end = j*n + ie;
		index_start = j*n + iw;

		r_f[index_end + 3 * m * n] = r_f[index_start + 3 * m * n];
		r_f[index_end + 7 * m * n] = r_f[index_start + 7 * m * n];
		r_f[index_end + 6 * m * n] = r_f[index_start + 6 * m * n];

		b_f[index_end + 3 * m * n] = b_f[index_start + 3 * m * n];
		b_f[index_end + 7 * m * n] = b_f[index_start + 7 * m * n];
		b_f[index_end + 6 * m * n] = b_f[index_start + 6 * m * n];

		// west boundary
		r_f[index_start + 1 * m * n] = r_f[index_end + 1 * m * n];
		r_f[index_start + 5 * m * n] = r_f[index_end + 5 * m * n];
		r_f[index_start + 8 * m * n] = r_f[index_end + 8 * m * n];

		b_f[index_start + 1 * m * n] = b_f[index_end + 1 * m * n];
		b_f[index_start + 5 * m * n] = b_f[index_end + 5 * m * n];
		b_f[index_start + 8 * m * n] = b_f[index_end + 8 * m * n];

	}

	// north-east corner
	r_f[(jn*n+ie) + 3 * m * n] = r_f[(jn*n+iw) + 3 * m * n];
	r_f[(jn*n+ie) + 4 * m * n] = r_f[(jn*n+iw) + 4 * m * n];
	r_f[(jn*n+ie) + 7 * m * n] = r_f[(jn*n+iw) + 7 * m * n];

	b_f[(jn*n+ie) + 3 * m * n] = b_f[(jn*n+iw) + 3 * m * n];
	b_f[(jn*n+ie) + 4 * m * n] = b_f[(jn*n+iw) + 4 * m * n];
	b_f[(jn*n+ie) + 7 * m * n] = b_f[(jn*n+iw) + 7 * m * n];

	// north-west corner
	r_f[(jn*n+iw) + 1 * m * n] = r_f[(jn*n+ie) + 1 * m * n];
	r_f[(jn*n+iw) + 4 * m * n] = r_f[(jn*n+ie) + 4 * m * n];
	r_f[(jn*n+iw) + 8 * m * n] = r_f[(jn*n+ie) + 8 * m * n];

	b_f[(jn*n+iw) + 1 * m * n] = b_f[(jn*n+ie) + 1 * m * n];
	b_f[(jn*n+iw) + 4 * m * n] = b_f[(jn*n+ie) + 4 * m * n];
	b_f[(jn*n+iw) + 8 * m * n] = b_f[(jn*n+ie) + 8 * m * n];

	// south-east corner
	r_f[(js*n+ie) + 2 * m * n] = r_f[(js*n+iw) + 2 * m * n];
	r_f[(js*n+ie) + 3 * m * n] = r_f[(js*n+iw) + 3 * m * n];
	r_f[(js*n+ie) + 6 * m * n] = r_f[(js*n+iw) + 6 * m * n];

	b_f[(js*n+ie) + 2 * m * n] = b_f[(js*n+iw) + 2 * m * n];
	b_f[(js*n+ie) + 3 * m * n] = b_f[(js*n+iw) + 3 * m * n];
	b_f[(js*n+ie) + 6 * m * n] = b_f[(js*n+iw) + 6 * m * n];

	// south-west corner
	r_f[(js*n+iw) + 2 * m * n] = r_f[(js*n+ie) + 2 * m * n];
	r_f[(js*n+iw) + 1 * m * n] = r_f[(js*n+ie) + 1 * m * n];
	r_f[(js*n+iw) + 5 * m * n] = r_f[(js*n+ie) + 5 * m * n];

	b_f[(js*n+iw) + 2 * m * n] = b_f[(js*n+ie) + 2 * m * n];
	b_f[(js*n+iw) + 1 * m * n] = b_f[(js*n+ie) + 1 * m * n];
	b_f[(js*n+iw) + 5 * m * n] = b_f[(js*n+ie) + 5 * m * n];
}

void peridicBoundaries3D(int n, int m, int h, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho){

	int index_end, index_start;
	int jn = m-1;
	int js = 0;
	int ie = n-1;
	int iw = 0;
	int kf = h-1;
	int kb = 0;
	int ms = m * n * h;
	for(int k = 0; k < h - 1; k++){
		for (int i=0; i < n; i++){
			// north boundary
			index_end = k * n * m + jn * n + i;
			index_start = k * n * m + js * n + i;

			r_f[index_end + 4 * ms] = r_f[index_start + 4 * ms];
			r_f[index_end + 9 * ms] = r_f[index_start + 9 * ms];
			r_f[index_end + 10 * ms] = r_f[index_start + 10 * ms];
			r_f[index_end + 16 * ms] = r_f[index_start + 16 * ms];
			r_f[index_end + 18 * ms] = r_f[index_start + 18 * ms];

			b_f[index_end + 4 * ms] = b_f[index_start + 4 * ms];
			b_f[index_end + 9 * ms] = b_f[index_start + 9 * ms];
			b_f[index_end + 10 * ms] = b_f[index_start + 10 * ms];
			b_f[index_end + 16 * ms] = b_f[index_start + 16 * ms];
			r_f[index_end + 18 * ms] = r_f[index_start + 18 * ms];


			//south boundary
			r_f[index_start + 3 * ms] = r_f[index_end + 3 * ms];
			r_f[index_start + 7 * ms] = r_f[index_end + 7 * ms];
			r_f[index_start + 8 * ms] = r_f[index_end + 8 * ms];
			r_f[index_start + 15 * ms] = r_f[index_end + 15 * ms];
			r_f[index_start + 17 * ms] = r_f[index_end + 17 * ms];

			b_f[index_start + 3 * ms] = b_f[index_end + 3 * ms];
			b_f[index_start + 7 * ms] = b_f[index_end + 7 * ms];
			b_f[index_start + 8 * ms] = b_f[index_end + 8 * ms];
			b_f[index_start + 15 * ms] = b_f[index_end + 15 * ms];
			b_f[index_start + 17 * ms] = b_f[index_end + 17 * ms];
		}

		for (int j=1; j < m-1; j++){
			// east boundary
			index_end = k * m * n + j*n + ie;
			index_start = k * m * n + j*n + iw;

			r_f[index_end + 2 * ms] = r_f[index_start + 2 * ms];
			r_f[index_end + 8 * ms] = r_f[index_start + 8 * ms];
			r_f[index_end + 10 * ms] = r_f[index_start + 10 * ms];
			r_f[index_end + 12 * ms] = r_f[index_start + 12 * ms];
			r_f[index_end + 14 * ms] = r_f[index_start + 14 * ms];

			b_f[index_end + 2 * ms] = b_f[index_start + 2 * ms];
			b_f[index_end + 8 * ms] = b_f[index_start + 8 * ms];
			b_f[index_end + 10 * ms] = b_f[index_start + 10 * ms];
			b_f[index_end + 12 * ms] = b_f[index_start + 12 * ms];
			b_f[index_end + 14 * ms] = b_f[index_start + 14 * ms];

			// west boundary
			r_f[index_start + 1 * ms] = r_f[index_end + 1 * ms];
			r_f[index_start + 7 * ms] = r_f[index_end + 7 * ms];
			r_f[index_start + 9 * ms] = r_f[index_end + 9 * ms];
			r_f[index_start + 11 * ms] = r_f[index_end + 11 * ms];
			r_f[index_start + 13 * ms] = r_f[index_end + 13 * ms];

			b_f[index_start + 1 * ms] = b_f[index_end + 1 * ms];
			b_f[index_start + 7 * ms] = b_f[index_end + 7 * ms];
			b_f[index_start + 9 * ms] = b_f[index_end + 9 * ms];
			r_f[index_start + 11 * ms] = r_f[index_end + 11 * ms];
			r_f[index_start + 13 * ms] = r_f[index_end + 13 * ms];

		}
	}

	for(int j = 0; j < m; j++){
		for(int i = 0; i < n; i++){
			index_end = kf * m * n + j*n + i;
			index_start = kb * m * n + j*n + i;

			//Front boundary
			r_f[index_end + 6 * ms] = r_f[index_start + 6 * ms];
			r_f[index_end + 13 * ms] = r_f[index_start + 13 * ms];
			r_f[index_end + 14 * ms] = r_f[index_start + 14 * ms];
			r_f[index_end + 17 * ms] = r_f[index_start + 17 * ms];
			r_f[index_end + 18 * ms] = r_f[index_start + 18 * ms];

			b_f[index_end + 6 * ms] = b_f[index_start + 6 * ms];
			b_f[index_end + 13 * ms] = b_f[index_start + 13 * ms];
			b_f[index_end + 14 * ms] = b_f[index_start + 14 * ms];
			b_f[index_end + 17 * ms] = b_f[index_start + 17 * ms];
			b_f[index_end + 18 * ms] = b_f[index_start + 18 * ms];

			// back boundary
			r_f[index_start + 5 * ms] = r_f[index_end + 5 * ms];
			r_f[index_start + 11 * ms] = r_f[index_end + 11 * ms];
			r_f[index_start + 12 * ms] = r_f[index_end + 12 * ms];
			r_f[index_start + 15 * ms] = r_f[index_end + 15 * ms];
			r_f[index_start + 16 * ms] = r_f[index_end + 16 * ms];

			b_f[index_start + 5 * ms] = b_f[index_end + 5 * ms];
			b_f[index_start + 11 * ms] = b_f[index_end + 11 * ms];
			b_f[index_start + 12 * ms] = b_f[index_end + 12 * ms];
			r_f[index_start + 15 * ms] = r_f[index_end + 15 * ms];
			r_f[index_start + 16 * ms] = r_f[index_end + 16 * ms];
		}
	}
}

void streamMP(int n, int m, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl){
	// stream on interior first
	int index,i,j;
	for (j=1;j < m-1;j++){
		for (i=1; i < n-1; i++){
			index = j*n+i;
			r_f[index] = r_fColl[index];
			r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
			r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
			r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
			r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
			r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];
			r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];
			r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];
			r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

			b_f[index] = b_fColl[index];
			b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
			b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
			b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
			b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
			b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];
			b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];
			b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];
			b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];
		}
	}
	for (i=1; i < n-1; i++){
		//north boundary
		j = m-1;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
		r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
		r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
		r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];
		r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
		b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
		b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
		b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];
		b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];

		//South boundary
		j = 0;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
		r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
		r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
		r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];
		r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
		b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
		b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
		b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];
		b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];
	}

	for (j=1;j < m-1;j++){
		//east
		i = n-1;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
		r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
		r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
		r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];
		r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
		b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
		b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
		b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];
		b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];

		//west
		i = 0;
		index = j*n+i;

		r_f[index] = r_fColl[index];
		r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
		r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
		r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
		r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];
		r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];

		b_f[index] = b_fColl[index];
		b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
		b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
		b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
		b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];
		b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];
	}

	// north-east corner
	i=n-1; j=m-1;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
	r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
	r_f[index + 5 * m * n] = r_fColl[((j-1) * n + i - 1) + 5 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
	b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
	b_f[index + 5 * m * n] = b_fColl[((j-1) * n + i - 1) + 5 * m * n];

	//north-west corner
	i=0; j=m-1;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 2 * m * n] = r_fColl[((j-1) * n + i) + 2 * m * n];
	r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
	r_f[index + 6 * m * n] = r_fColl[((j-1) * n + i + 1) + 6 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 2 * m * n] = b_fColl[((j-1) * n + i) + 2 * m * n];
	b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
	b_f[index + 6 * m * n] = b_fColl[((j-1) * n + i + 1) + 6 * m * n];

	// south-east corner
	i=n-1; j=0;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 1 * m * n] = r_fColl[(j*n+i-1) + 1 * m * n];
	r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
	r_f[index + 8 * m * n] = r_fColl[((j+1) * n + i - 1) + 8 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 1 * m * n] = b_fColl[(j*n+i-1) + 1 * m * n];
	b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
	b_f[index + 8 * m * n] = b_fColl[((j+1) * n + i - 1) + 8 * m * n];

	// south-west corner
	i=0; j=0;
	index = j*n+i;

	r_f[index] = r_fColl[index];
	r_f[index + 3 * m * n] = r_fColl[(j*n + i + 1) + 3 * m * n];
	r_f[index + 4 * m * n] = r_fColl[((j+1) * n + i) + 4 * m * n];
	r_f[index + 7 * m * n] = r_fColl[((j+1) * n + i + 1) + 7 * m * n];

	b_f[index] = b_fColl[index];
	b_f[index + 3 * m * n] = b_fColl[(j*n + i + 1) + 3 * m * n];
	b_f[index + 4 * m * n] = b_fColl[((j+1) * n + i) + 4 * m * n];
	b_f[index + 7 * m * n] = b_fColl[((j+1) * n + i + 1) + 7 * m * n];

}

void streamMP3D(int n, int m, int h, FLOAT_TYPE *r_f, FLOAT_TYPE *b_f, FLOAT_TYPE *r_fColl, FLOAT_TYPE *b_fColl, bool *stream){
	// stream on interior first
	int index,i,j,k;
	int ms = m*n*h;
	int c3D[19] = { 0, -1, 1, -1 * n, n, -m * n, +m * n, -1 * n - 1, -1 * n + 1,
			n - 1, n + 1, -m * n - 1, -m * n + 1, +m * n - 1,
			m * n + 1, -m * n - n, -m * n + n, m * n - n, m * n + n };

	for(k = 1; k < h - 1; k++){
		for (j=1;j < m-1;j++){
			for (i=1; i < n-1; i++){
				index = k * m * n + j*n+i;
				r_f[index] = r_fColl[index];
				b_f[index] = b_fColl[index];
				for(int dir = 1; dir < 19; dir++){
					r_f[index + dir * ms] = (stream[index+	(dir-1) * ms]	==	1)	?	r_fColl[index + dir * ms + c3D[dir]]:	r_f[index + dir * ms];
					b_f[index + dir * ms] = (stream[index+	(dir-1) * ms]	==	1)	?	b_fColl[index + dir * ms + c3D[dir]]:	b_f[index + dir * ms];
				}
			}
		}
	}
}

void resetArrays(FLOAT_TYPE *color_gradient, int n, int m){
	for(int i = 0; i < m * n *2; i++){
		color_gradient[i] = 0.0;
	}
}

FLOAT_TYPE* convertArray(int n, int m, FLOAT_TYPE *arr){
	FLOAT_TYPE *result = createHostArrayFlt(n*m, ARRAY_NONE);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			result[j*n+i] = arr[i*m+j];
		}
	}

	return result;
}

void updateSurfaceTension(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, FLOAT_TYPE control_param,
		FLOAT_TYPE st_predicted, FLOAT_TYPE *st_error, int iteration, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha, FLOAT_TYPE bubble_radius, int n, int m){

	FLOAT_TYPE chi, p_in = 0.0, p_out = 0.0, st_laplace, r;
	int index, index_aux1 = 0, index_aux2 = 0;
	for(int j = 1; j < m-1; j++){
		for(int i = 1; i < n-1; i++){
			index = j * n + i;
			r = r_rho[index] + b_rho[index];
			// p_in and p_out for the surface tension
			chi=(r_rho[index]-b_rho[index])/r;
			if (chi >= control_param){
				index_aux1++;
				p_in += r_rho[index];
			}
			else if (chi <= -control_param){
				index_aux2++;
				p_out+=b_rho[index];
			}
		}
	}

	p_in=(3.0/5.0)*(1.0-r_alpha)*p_in/index_aux1;      // pressure average inside the bubble
	p_out=(3.0/5.0)*(1.0-b_alpha)*p_out/index_aux2;   // pressure average outside the bubble
	st_laplace=bubble_radius*(p_in-p_out);

	st_error[iteration]=abs(st_predicted-st_laplace)/(st_predicted)*100.0;
}

FLOAT_TYPE calculateSurfaceTension(FLOAT_TYPE p_in_mean, FLOAT_TYPE p_out_mean, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha,
		FLOAT_TYPE bubble_radius, FLOAT_TYPE st_predicted){

	FLOAT_TYPE st_laplace;
	p_in_mean=(3.0/5.0)*(1.0-r_alpha)*p_in_mean;      // pressure average inside the bubble
	p_out_mean=(3.0/5.0)*(1.0-b_alpha)*p_out_mean;   // pressure average outside the bubble
	st_laplace=bubble_radius * (p_in_mean-p_out_mean);

	return abs(st_predicted-st_laplace)/(st_predicted)*100.0;
}

FLOAT_TYPE calculateSurfaceTension3D(FLOAT_TYPE p_in_mean, FLOAT_TYPE p_out_mean, FLOAT_TYPE r_alpha, FLOAT_TYPE b_alpha,
		FLOAT_TYPE bubble_radius, FLOAT_TYPE st_predicted){

	FLOAT_TYPE st_laplace;
	p_in_mean=(1.0/2.0)*(1.0-r_alpha)*p_in_mean;      // pressure average inside the bubble
	p_out_mean=(1.0/2.0)*(1.0-b_alpha)*p_out_mean;   // pressure average outside the bubble
	st_laplace=bubble_radius * (p_in_mean-p_out_mean);

	return abs(st_predicted-st_laplace)/(st_predicted)*100.0;
}

void validateCoalescenceCase(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, int n, int m, FLOAT_TYPE radius, int h){
	int j,k;
	if(m % 2 == 0)
		j = m/2;
	else
		j = (m+1) / 2;
	if(h % 2 == 0)
		k = h / 2;
	else
		k = (h+1) / 2;

	FLOAT_TYPE rho;
	FLOAT_TYPE aux = 0.0;
	int ms = m * n;
	for(int i = 0; i < n; i++){
		rho = r_rho[k * ms + j * n + i] + b_rho[k * ms + j * n + i];
		if((r_rho[k * ms + j * n + i] - b_rho[k * ms + j * n + i]) / rho > 0.0){
			aux += 1 * (r_rho[k * ms + j * n + i] - b_rho[k * ms + j * n + i]) / rho;
		}
	}

	FLOAT_TYPE predicted;
	if(h > 0)
		predicted = radius * pow(2.0, (1.0 / 3.0));
	else
		predicted = radius * sqrt(2.0);

	aux /= n * 2.0;
	printf("Predicted radius: "FLOAT_FORMAT" Final radius: "FLOAT_FORMAT"  (physical units)\n", predicted, aux);
	printf("Radius error % = "FLOAT_FORMAT" \n", (abs(aux - predicted) / predicted) * 100.0);
}

void analyticalCouette(FLOAT_TYPE kappa, FLOAT_TYPE *y, int m, int n, FLOAT_TYPE *analytical, FLOAT_TYPE ux, int h){

	int j_start, i, k;
	if(m % 2 == 0)
		j_start = m/2;
	else
		j_start = (m+1) / 2;
	if(n % 2 == 0)
		i = n/2;
	else
		i = (n+1) / 2;

	if(h % 2 == 0)
		k = h / 2;
	else
		k = (h+1) / 2;

	int ms = n * m;
	for(int j = j_start; j < m; j++){
		analytical[j] = (2.0 * y[k * ms + j*n + i] / (kappa + 1.0) + (kappa - 1.0) / (kappa + 1.0)) * ux;
	}
	for(int j = 0; j < j_start; j++){
		analytical[j] = (2.0 * kappa * y[k * ms + j*n + i] / (kappa + 1.0)) * ux;
	}
}

void deformingBubbleValid(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, int n, int m, int h){
	int j,k;
	if(m % 2 == 0)
		j = m/2;
	else
		j = (m+1) / 2;
	if(h % 2 == 0)
		k = h / 2;
	else
		k = (h+1) / 2;

	FLOAT_TYPE rho;
	FLOAT_TYPE aux = 0.0;
	int ms = m * n;
	for(int i = 0; i < n; i++){
		rho = r_rho[k * ms + j * n + i] + b_rho[k * ms + j * n + i];
		if((r_rho[k * ms + j * n + i] - b_rho[k * ms + j * n + i]) / rho > 0.0){
			aux += (r_rho[k * ms + j * n + i] - b_rho[k * ms + j * n + i]) / rho;
		}
	}

	aux /= 2.0;
	FLOAT_TYPE predicted;
	if(h > 0){
		predicted = n / 2.0 / pow((4.0/3.0 * M_PI),(1.0/3.0));
	}else{
		predicted = n / (2.0 * sqrt(M_PI));
	}
	printf("Predicted radius: "FLOAT_FORMAT" Final radius: "FLOAT_FORMAT"\n", predicted, aux);
	printf("Radius error % = "FLOAT_FORMAT" \n", (abs(aux - predicted) / predicted) * 100.0);
}

void initInletVelocity(FLOAT_TYPE *u, FLOAT_TYPE *v, FLOAT_TYPE u_veloc, FLOAT_TYPE v_veloc, int n, int m){

	for(int i = 0; i < n; i++){
		u[i] = 0.0;
		v[i] = 0.0;

		u[(m-1) * n + i] = u_veloc;
		v[(m-1) * n + i] = v_veloc;
	}

}

FLOAT_TYPE getMaxYOscilating(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, int n, int m, FLOAT_TYPE *nodeY){

	int i, j_start, j;
	if(n % 2 == 0)
		i = n/2;
	else
		i = (n+1) / 2;

	if(m % 2 == 0)
		j_start = m/2;
	else
		j_start = (m+1) / 2;

	FLOAT_TYPE rho;
	for(j = j_start; j < m; j++){

		rho = r_rho[j * n + i] + b_rho[j * n + i];
		if((r_rho[j * n + i] - b_rho[j * n + i]) / rho < 0.0){
			break;
		}
	}
	FLOAT_TYPE phi1, phi2;
	phi1 = (r_rho[j * n + i] - b_rho[j * n + i]) / rho;
	rho = r_rho[(j-1) * n + i] + b_rho[(j-1) * n + i];
	phi2 = (r_rho[(j-1) * n + i] - b_rho[(j-1) * n + i]) / rho;

	FLOAT_TYPE aux_m = (phi1 - phi2) / (nodeY[j * n + i] - nodeY[(j - 1) * n + i]);
	FLOAT_TYPE aux_b =  phi1 - aux_m * nodeY[j * n + i];
	return -aux_b / aux_m;
}

FLOAT_TYPE getMinYRT(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, int n, int m, FLOAT_TYPE *nodeY){

	int i, j_start, j;
	if(n % 2 == 0)
		i = n/2;
	else
		i = (n+1) / 2;

	if(m % 2 == 0)
		j_start = m/2;
	else
		j_start = (m+1) / 2;

	FLOAT_TYPE rho;
	for(j = j_start; j >= 0; j--){

		rho = r_rho[j * n + i] + b_rho[j * n + i];
		if((r_rho[j * n + i] - b_rho[j * n + i]) / rho < 0.0){
			break;
		}
	}
	FLOAT_TYPE phi1, phi2;
	phi1 = (r_rho[j * n + i] - b_rho[j * n + i]) / rho;
	rho = r_rho[(j+1) * n + i] + b_rho[(j+1) * n + i];
	phi2 = (r_rho[(j+1) * n + i] - b_rho[(j+1) * n + i]) / rho;

	FLOAT_TYPE aux_m = (phi1 - phi2) / (nodeY[j * n + i] - nodeY[(j + 1) * n + i]);
	FLOAT_TYPE aux_b =  phi1 - aux_m * nodeY[j * n + i];
	return -aux_b / aux_m;
}

FLOAT_TYPE validateOscilating(FLOAT_TYPE *r_rho, FLOAT_TYPE *b_rho, int n, int m, FLOAT_TYPE *extremes, int size,
		FLOAT_TYPE ST_predicted, FLOAT_TYPE r_density, FLOAT_TYPE b_density){
	int j;
	if(m % 2 == 0)
		j = m/2;
	else
		j = (m+1) / 2;

	FLOAT_TYPE rho;
	FLOAT_TYPE aux = 0.0;
	for(int i = 0; i < n; i++){
		rho = r_rho[j * n + i] + b_rho[j * n + i];
		if((r_rho[j * n + i] - b_rho[j * n + i]) / rho > 0.0){
			aux += (r_rho[j * n + i] - b_rho[j * n + i]) / rho;
		}
	}

	aux /= 2.0;
	printf("\n\nMax and Min values\n\n");
	int max_iter[6];
	int min_iter[6];
	int max_counter = 0, min_counter = 0;
	bool sign = false;
	for(int i = 500; i < size-1; i++){
		if(sign){
			if((extremes[i+1] - extremes[i-1]) / 2.0 < 0.0){
				printf("Max found at %d\n", i);
				max_iter[max_counter] = i;
				max_counter++;
				sign = false;
			}
		}
		else{
			if((extremes[i+1] - extremes[i-1]) / 2.0 > 0.0){
				printf("Min found at %d\n", i);
				min_iter[min_counter] = i;
				min_counter++;
				sign = true;
			}
		}
	}

	FLOAT_TYPE theoretical_omega = sqrt(6.0 * ST_predicted / ( (r_density + b_density) *aux*aux*aux ) );
	FLOAT_TYPE numerical_omega = M_PI / (min_iter[1] - max_iter[0]);
	printf("Theoretical "FLOAT_FORMAT" vs Numerical "FLOAT_FORMAT"\n", theoretical_omega, numerical_omega);
	return (abs(theoretical_omega - numerical_omega) / theoretical_omega )* 100.0;
}

void analyticalPoiseuille(int m, int n, FLOAT_TYPE *analytical, FLOAT_TYPE r_density, FLOAT_TYPE b_density,
		FLOAT_TYPE r_visc, FLOAT_TYPE b_visc, FLOAT_TYPE g, FLOAT_TYPE *y){

	FLOAT_TYPE r_mu, b_mu;
	r_mu = r_density * r_visc;
	b_mu = b_density * b_visc;
	int j_start, i;
	if(m % 2 == 0)
		j_start = m/2;
	else
		j_start = (m+1) / 2;
	if(n % 2 == 0)
		i = n/2;
	else
		i = (n+1) / 2;

	for(int j = j_start; j < m; j++){
		analytical[j] = -g * y[j * n + i] * y[j * n + i] / (2.0 * r_mu) + g * ((3.0 * b_mu + r_mu ) * y[j * n + i] / (4.0 * r_mu * (r_mu + b_mu))) +
				g / (2.0 * r_mu)  - g * ((3.0 * b_mu + r_mu) / (4.0 * r_mu * (r_mu + b_mu)));
	}
	for(int j = 0; j < j_start; j++){
		analytical[j] = -g * y[j * n + i] * y[j * n + i] / (2.0 * b_mu) + g * ((3.0 * b_mu + r_mu) * y[j * n + i]) / (4.0 * b_mu * (r_mu + b_mu));
	}

}
