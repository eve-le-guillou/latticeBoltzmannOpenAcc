#include "CpuFunctions.h"
#include "BcMacros.h"
#include "BcMacros3D.h"
#include "CpuConstants.h"

void cpuUpdateMacro2D(int *fluid_d, FLOAT_TYPE* rho_d,
                                 FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, int *bcMask_d, FLOAT_TYPE* drag_d,
                                 FLOAT_TYPE* lift_d,
                                 FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d, int dlBoundaryId) {

    int ms = depth * length;

    FLOAT_TYPE r, u, v;

    for (int ind = 0; ind < ms; ind++) {
        if (fluid_d[ind] == 1) {
            r = u = v = 0.0;
            r = f_d[ind] + f_d[ind + ms] + f_d[ind + 2 * ms] + f_d[ind + 3 * ms]
                + f_d[ind + 4 * ms] + f_d[ind + 5 * ms] + f_d[ind + 6 * ms]
                + f_d[ind + 7 * ms] + f_d[ind + 8 * ms];
            u = f_d[ind + ms] - f_d[ind + 3 * ms] + f_d[ind + 5 * ms]
                - f_d[ind + 6 * ms] - f_d[ind + 7 * ms] + f_d[ind + 8 * ms];
            v = f_d[ind + 2 * ms] - f_d[ind + 4 * ms] + f_d[ind + 5 * ms]
                + f_d[ind + 6 * ms] - f_d[ind + 7 * ms] - f_d[ind + 8 * ms];

            rho_d[ind] = r;
            u_d[ind] = u / r;
            ///@todo code: probably should handle outlet on other sides
            v_d[ind] = ((bcMask_d[ind] & BC_OUTL_E) == BC_OUTL_E) ? 0.0 : v / r;

            //   DRAG/LIFT FORCE
            if (dlBoundaryId
                != 0&& (bcMask_d[ind] & BND_ID_ALL) == BOUND_ID(dlBoundaryId)) {
                // printf("draglift: %d\n",ind);
                drag_d[ind] = 0.33333333 * r * (20 - coordX_d[ind]) * 0.2;
                lift_d[ind] = 0.33333333 * r * (20 - coordY_d[ind]) * 0.2;
            }
        }
    }
}

void cpuUpdateMacro2DCG(FLOAT_TYPE* rho_d,
                                   FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* r_f_d, FLOAT_TYPE* b_f_d, FLOAT_TYPE *f_d, FLOAT_TYPE* r_rho_d,
                                   FLOAT_TYPE* b_rho_d, FLOAT_TYPE *p_in_d, FLOAT_TYPE *p_out_d,
                                   int *num_in_d, int *num_out_d, int *cg_direction, int test_case) {
    int ms = depth * length;

    FLOAT_TYPE r_r, b_r, u, v, r, chi;
    FLOAT_TYPE aux1, mean_nu, omega_eff;
    for (int ind = 0; ind < ms; ind++) {
        //necessary because of sum
        if(test_case == 1){
            p_in_d[ind] = 0;
            p_out_d[ind] = 0;
            num_in_d[ind] = 0;
            num_out_d[ind] = 0;
        }

        if (cg_direction[ind] == 0 || (test_case != 6) || (cg_direction[ind] == 3 || cg_direction[ind] == 4)) {

            aux1 = r_rho_d[ind] / (rho_d[ind] * r_viscosity) + b_rho_d[ind] /(rho_d[ind] * b_viscosity);
            mean_nu = 1.0/aux1;
            omega_eff = 1.0/(3.0*mean_nu+0.5);

            r_r = b_r = u = v = 0.0;
            r_r = r_f_d[ind] +
                  r_f_d[ind + ms] +
                  r_f_d[ind + 2 * ms] +
                  r_f_d[ind + 3 * ms] +
                  r_f_d[ind + 4 * ms] +
                  r_f_d[ind + 5 * ms] +
                  r_f_d[ind + 6 * ms]	+
                  r_f_d[ind + 7 * ms] +
                  r_f_d[ind + 8 * ms];
            b_r = b_f_d[ind] +
                  b_f_d[ind + ms] +
                  b_f_d[ind + 2 * ms] +
                  b_f_d[ind + 3 * ms] +
                  b_f_d[ind + 4 * ms] +
                  b_f_d[ind + 5 * ms] +
                  b_f_d[ind + 6 * ms] +
                  b_f_d[ind + 7 * ms] +
                  b_f_d[ind + 8 * ms];

            f_d[ind] = r_f_d[ind] + b_f_d[ind];
            f_d[ind + 1 * ms] = r_f_d[ind + 1 * ms] + b_f_d[ind + 1 * ms];
            f_d[ind + 2 * ms] = r_f_d[ind + 2 * ms] + b_f_d[ind + 2 * ms];
            f_d[ind + 3 * ms] = r_f_d[ind + 3 * ms] + b_f_d[ind + 3 * ms];
            f_d[ind + 4 * ms] = r_f_d[ind + 4 * ms] + b_f_d[ind + 4 * ms];
            f_d[ind + 5 * ms] = r_f_d[ind + 5 * ms] + b_f_d[ind + 5 * ms];
            f_d[ind + 6 * ms] = r_f_d[ind + 6 * ms] + b_f_d[ind + 6 * ms];
            f_d[ind + 7 * ms] = r_f_d[ind + 7 * ms] + b_f_d[ind + 7 * ms];
            f_d[ind + 8 * ms] = r_f_d[ind + 8 * ms] + b_f_d[ind + 8 * ms];

            r_rho_d[ind] = r_r;
            b_rho_d[ind] = b_r;
            r = r_r + b_r;
            rho_d[ind] = r;



            u = (r_f_d[ind + ms] + b_f_d[ind + ms]) -
                (r_f_d[ind + 3 * ms] + b_f_d[ind + 3 * ms]) +
                (r_f_d[ind + 5 * ms] + b_f_d[ind + 5 * ms]) -
                (r_f_d[ind + 6 * ms] + b_f_d[ind + 6 * ms]) -
                (r_f_d[ind + 7 * ms] + b_f_d[ind + 7 * ms]) +
                (r_f_d[ind + 8 * ms] + b_f_d[ind + 8 * ms]);

            v = (r_f_d[ind + 2 * ms] + b_f_d[ind + 2 * ms]) -
                (r_f_d[ind + 4 * ms] + b_f_d[ind + 4 * ms]) +
                (r_f_d[ind + 5 * ms] + b_f_d[ind + 5 * ms]) +
                (r_f_d[ind + 6 * ms] + b_f_d[ind + 6 * ms]) -
                (r_f_d[ind + 7 * ms] + b_f_d[ind + 7 * ms]) -
                (r_f_d[ind + 8 * ms] + b_f_d[ind + 8 * ms]);


            u_d[ind] = u / r + external_force * g / (r * omega_eff);
            v_d[ind] = v / r + (1-external_force) * g / omega_eff;

            if(test_case == 1){
                // p_in and p_out for the surface tension
                chi=(r_r-b_r)/r;

                if (chi >= control_param){
                    num_in_d[ind] = 1;
                    p_in_d[ind] = r_r;
                }
                else if (chi <= -control_param){
                    num_out_d[ind] = 1;
                    p_out_d[ind] = b_r;
                }
            }
        }
    }
}

void cpuUpdateMacro3D(int *fluid_d, FLOAT_TYPE* rho_d,
                                 FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* w_d, int* bcBoundId_d,
                                 FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* coordZ_d,
                                 FLOAT_TYPE* f_d, FLOAT_TYPE g_arg, unsigned long long *bcMask_d,int updateInltOutl)
{
    int ms = depth * length * height;

    FLOAT_TYPE r, rU, rV, rW;

    for (int ind = 0; ind < ms; ind++) {
        if (fluid_d[ind] == 1
            && (!(((bcMask_d[ind] & BC3D_OUTL_B) ==BC3D_INLT_B))||updateInltOutl)) {
            r = rU = rV = rW = 0.0;
            r = f_d[ind] +
                f_d[ind + ms] +
                f_d[ind + 2 * ms] +
                f_d[ind + 3 * ms] +
                f_d[ind + 4 * ms] +
                f_d[ind + 5 * ms] +
                f_d[ind + 6 * ms] +
                f_d[ind + 7 * ms] +
                f_d[ind + 8 * ms] +
                f_d[ind + 9 * ms] +
                f_d[ind + 10 * ms] +
                f_d[ind + 11 * ms] +
                f_d[ind + 12 * ms] +
                f_d[ind + 13 * ms] +
                f_d[ind + 14 * ms] +
                f_d[ind + 15 * ms] +
                f_d[ind + 16 * ms] +
                f_d[ind + 17 * ms] +
                f_d[ind + 18 * ms];

            rU = f_d[ind + ms] -
                 f_d[ind + 2 * ms] +
                 f_d[ind + 7 * ms] -
                 f_d[ind + 8 * ms] +
                 f_d[ind + 9 * ms] -
                 f_d[ind + 10 * ms] +
                 f_d[ind + 11 * ms] -
                 f_d[ind + 12 * ms] +
                 f_d[ind + 13 * ms] -
                 f_d[ind + 14 * ms];

            rV = f_d[ind + 3 * ms] -
                 f_d[ind + 4 * ms] +
                 f_d[ind + 7 * ms] +
                 f_d[ind + 8 * ms] -
                 f_d[ind + 9 * ms] -
                 f_d[ind + 10 * ms] +
                 f_d[ind + 15 * ms] -
                 f_d[ind + 16 * ms] +
                 f_d[ind + 17 * ms] -
                 f_d[ind + 18 * ms];

            rW = f_d[ind + 5 * ms] -
                 f_d[ind + 6 * ms] +
                 f_d[ind + 11 * ms] +
                 f_d[ind + 12 * ms] -
                 f_d[ind + 13 * ms] -
                 f_d[ind + 14 * ms] +
                 f_d[ind + 15 * ms] +
                 f_d[ind + 16 * ms] -
                 f_d[ind + 17 * ms] -
                 f_d[ind + 18 * ms];

            rho_d[ind] = r;
            u_d[ind] = rU / r + g_arg / (omega);
            v_d[ind] = rV / r;
            w_d[ind] = rW / r;
        }
    }

}

void cpuUpdateMacro3DCG(int *fluid_d, FLOAT_TYPE* rho_d,
                                   FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* w_d, int* bcBoundId_d,
                                   FLOAT_TYPE* f_d, FLOAT_TYPE g_arg, unsigned long long *bcMask_d,int updateInltOutl, FLOAT_TYPE* r_f_d, FLOAT_TYPE* b_f_d, FLOAT_TYPE* r_rho_d,
                                   FLOAT_TYPE* b_rho_d, FLOAT_TYPE *p_in_d, FLOAT_TYPE *p_out_d,
                                   int *num_in_d, int *num_out_d, int test_case)
{
   int ms = depth * length * height;

    FLOAT_TYPE r_r, b_r, r, rU, rV, rW, aux1, mean_nu, omega_eff;

    for (int ind = 0; ind < ms; ind++) {
        //necessary because of sum
        if(test_case == 1){
            p_in_d[ind] = 0;
            p_out_d[ind] = 0;
            num_in_d[ind] = 0;
            num_out_d[ind] = 0;
        }
        if (fluid_d[ind] == 1 && (!(((bcMask_d[ind] & BC3D_OUTL_B) ==BC3D_INLT_B))||updateInltOutl)) {

            aux1 = r_rho_d[ind] / (rho_d[ind] * r_viscosity) + b_rho_d[ind] /(rho_d[ind] * b_viscosity);
            mean_nu = 1.0/aux1;
            omega_eff = 1.0/(3.0*mean_nu+0.5);

            r_r = b_r = rU = rV = rW = 0.0;
            r_r = r_f_d[ind] +
                  r_f_d[ind + ms] +
                  r_f_d[ind + 2 * ms] +
                  r_f_d[ind + 3 * ms] +
                  r_f_d[ind + 4 * ms] +
                  r_f_d[ind + 5 * ms] +
                  r_f_d[ind + 6 * ms] +
                  r_f_d[ind + 7 * ms] +
                  r_f_d[ind + 8 * ms] +
                  r_f_d[ind + 9 * ms] +
                  r_f_d[ind + 10 * ms] +
                  r_f_d[ind + 11 * ms] +
                  r_f_d[ind + 12 * ms] +
                  r_f_d[ind + 13 * ms] +
                  r_f_d[ind + 14 * ms] +
                  r_f_d[ind + 15 * ms] +
                  r_f_d[ind + 16 * ms] +
                  r_f_d[ind + 17 * ms] +
                  r_f_d[ind + 18 * ms];

            b_r = b_f_d[ind] +
                  b_f_d[ind + ms] +
                  b_f_d[ind + 2 * ms] +
                  b_f_d[ind + 3 * ms] +
                  b_f_d[ind + 4 * ms] +
                  b_f_d[ind + 5 * ms] +
                  b_f_d[ind + 6 * ms] +
                  b_f_d[ind + 7 * ms] +
                  b_f_d[ind + 8 * ms] +
                  b_f_d[ind + 9 * ms] +
                  b_f_d[ind + 10 * ms] +
                  b_f_d[ind + 11 * ms] +
                  b_f_d[ind + 12 * ms] +
                  b_f_d[ind + 13 * ms] +
                  b_f_d[ind + 14 * ms] +
                  b_f_d[ind + 15 * ms] +
                  b_f_d[ind + 16 * ms] +
                  b_f_d[ind + 17 * ms] +
                  b_f_d[ind + 18 * ms];

            r_rho_d[ind] = r_r;
            b_rho_d[ind] = b_r;

            f_d[ind] = r_f_d[ind] + b_f_d[ind];
            f_d[ind + ms] = r_f_d[ind + ms] + b_f_d[ind + ms];
            f_d[ind + 2 * ms] = r_f_d[ind + 2 * ms] + b_f_d[ind + 2 * ms];
            f_d[ind + 3 * ms] = r_f_d[ind + 3 * ms] + b_f_d[ind + 3 * ms];
            f_d[ind + 4 * ms] = r_f_d[ind + 4 * ms] + b_f_d[ind + 4 * ms];
            f_d[ind + 5 * ms] = r_f_d[ind + 5 * ms] + b_f_d[ind + 5 * ms];
            f_d[ind + 6 * ms] = r_f_d[ind + 6 * ms] + b_f_d[ind + 6 * ms];
            f_d[ind + 7 * ms] = r_f_d[ind + 7 * ms] + b_f_d[ind + 7 * ms];
            f_d[ind + 8 * ms] = r_f_d[ind + 8 * ms] + b_f_d[ind + 8 * ms];
            f_d[ind + 9 * ms] = r_f_d[ind + 9 * ms] + b_f_d[ind + 9 * ms];
            f_d[ind + 10 * ms] = r_f_d[ind + 10 * ms] + b_f_d[ind + 10 * ms];
            f_d[ind + 11 * ms] = r_f_d[ind + 11 * ms] + b_f_d[ind + 11 * ms];
            f_d[ind + 12 * ms] = r_f_d[ind + 12 * ms] + b_f_d[ind + 12 * ms];
            f_d[ind + 13 * ms] = r_f_d[ind + 13 * ms] + b_f_d[ind + 13 * ms];
            f_d[ind + 14 * ms] = r_f_d[ind + 14 * ms] + b_f_d[ind + 14 * ms];
            f_d[ind + 15 * ms] = r_f_d[ind + 15 * ms] + b_f_d[ind + 15 * ms];
            f_d[ind + 16 * ms] = r_f_d[ind + 16 * ms] + b_f_d[ind + 16 * ms];
            f_d[ind + 17 * ms] = r_f_d[ind + 17 * ms] + b_f_d[ind + 17 * ms];
            f_d[ind + 18 * ms] = r_f_d[ind + 18 * ms] + b_f_d[ind + 18 * ms];

            rU = f_d[ind + ms] -
                 f_d[ind + 2 * ms] +
                 f_d[ind + 7 * ms] -
                 f_d[ind + 8 * ms] +
                 f_d[ind + 9 * ms] -
                 f_d[ind + 10 * ms] +
                 f_d[ind + 11 * ms] -
                 f_d[ind + 12 * ms] +
                 f_d[ind + 13 * ms] -
                 f_d[ind + 14 * ms];

            rV = f_d[ind + 3 * ms] -
                 f_d[ind + 4 * ms] +
                 f_d[ind + 7 * ms] +
                 f_d[ind + 8 * ms] -
                 f_d[ind + 9 * ms] -
                 f_d[ind + 10 * ms] +
                 f_d[ind + 15 * ms] -
                 f_d[ind + 16 * ms] +
                 f_d[ind + 17 * ms] -
                 f_d[ind + 18 * ms];

            rW = f_d[ind + 5 * ms] -
                 f_d[ind + 6 * ms] +
                 f_d[ind + 11 * ms] +
                 f_d[ind + 12 * ms] -
                 f_d[ind + 13 * ms] -
                 f_d[ind + 14 * ms] +
                 f_d[ind + 15 * ms] +
                 f_d[ind + 16 * ms] -
                 f_d[ind + 17 * ms] -
                 f_d[ind + 18 * ms];

            r = r_r + b_r;

            rho_d[ind] = r;
            u_d[ind] =  rU / r + external_force * g / (r * omega_eff);
            v_d[ind] = rV / r + (1-external_force) * g / omega_eff;
            w_d[ind] = rW / r;

            if(test_case == 1){
                // p_in and p_out for the surface tension
                FLOAT_TYPE chi=(r_r-b_r)/r;

                if (chi >= control_param){
                    num_in_d[ind] = 1;
                    p_in_d[ind] = r_r;
                }
                else if (chi <= -control_param){
                    num_out_d[ind] = 1;
                    p_out_d[ind] = b_r;
                }
            }
        }
    }
}
