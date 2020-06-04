/**
 * Header for variables stored in GPU constant memory
 * @file GpuConstants.h
 * @author Adam Koleszar (adam.koleszar@gmail.com) - implementation of functions used in the 2D solver
 * @author Alfonso Aguilar (a.aguilar-pontes@cranfield.ac.uk) - implementation of functions used in the 3D solver
 * @details These constants need to be extern in other files because they can only be declared in
 * one source file. If other object files need these variable, just include this header and don't
 * forget to compile them with -rdc=true
 *
 * Velocity unit vector components (#cx_d, #cy_d)
 * @verbatim
   (-1,1)   (0,1)   (1,1)
          6   2   5
            \ | /
  (-1,0) 3 -(0,0)- 1 (1,0)
            / | \
          7   4   8
   (-1,-1)  (0,-1)  (1,-1) @endverbatim
 * Lattice weights (#w_d)
 * @verbatim
    (1/36)   (1/9)   (1/36)
           6   2   5
             \ | /
    (1/9) 3 -(4/9)- 1 (1/9)
             / | \
           7   4   8
    (1/36)   (1/9)   (1/36) @endverbatim
 * Opposite lattices (#opp_d)
 * @verbatim
        (8)   (4)   (7)
           6   2   5
             \ | /
       (1) 3 -(0)- 1 (3)
             / | \
           7   4   8
        (5)   (2)   (6) @endverbatim
 */
#ifndef GPU_CONSTANTS_H
#define GPU_CONSTANTS_H

#ifndef RELEASE

extern __constant__ InletProfile inletProfile_d;     ///< inlet profile
extern __constant__ BoundaryType boundaryType_d;     ///< boundary type
extern __constant__ OutletProfile outletProfile_d;   ///< outlet profile
extern __constant__ FLOAT_TYPE rhoIn_d;              ///< input density
extern __constant__ FLOAT_TYPE uIn_d;                ///< input velocity x
extern __constant__ FLOAT_TYPE vIn_d;                ///< input velocity y
extern __constant__ int dlBoundaryId_d;              ///< boundary ID
extern __constant__ int depth_d;                     ///< number of rows Y
extern __constant__ int length_d;                    ///< number of columns X
extern __constant__ int height_d;                    ///< number of layers Z
extern __constant__ FLOAT_TYPE delta_d;              ///< grid spacing
extern __constant__ FLOAT_TYPE minInletCoordY_d;     ///< maximum inlet coordinate y
extern __constant__ FLOAT_TYPE maxInletCoordY_d;     ///< minimum inlet coordinate y
extern __constant__ FLOAT_TYPE minInletCoordZ_d;     ///< maximum inlet coordinate z
extern __constant__ FLOAT_TYPE maxInletCoordZ_d;     ///< minimum inlet coordinate z
extern __constant__ FLOAT_TYPE omega_d;            ///< collision frequency for D2Q9 \f$ \omega = \frac{1}{3\nu + 0.5} \f$
extern __constant__ FLOAT_TYPE omegaA_d;           ///< asymmetric collision frequency \f$ \omega_a = \frac{8(2-\omega)}{8-\omega} \f$
extern __constant__ FLOAT_TYPE g_d;
//#### 2D d2q9 ####//
extern __constant__ int cx2D_d[9];                   ///< velocity x unit vector components
extern __constant__ int cy2D_d[9];                   ///< velocity y unit vector components
extern __constant__ int c2D_d[9];                    ///< direction offset levels
extern __constant__ int opp2D_d[9];                  ///< opposite lattice offset
extern __constant__ FLOAT_TYPE w2D_d[9];             ///< lattice weights
extern __constant__ FLOAT_TYPE velMomMap2D_d[81];    ///< MRT constants: mapping between velocity and momentum space \f$ \mathbf{M} \f$
extern __constant__ FLOAT_TYPE momCollMtx2D_d[81];   ///< MRT constants: collision matrix in momentum space \f$ \mathbf{M}^{-1}\mathbf{S} \f$

//#### Color Gradient ####//
extern __constant__ FLOAT_TYPE beta_d;
extern __constant__ FLOAT_TYPE g_limit_d;
extern __constant__ FLOAT_TYPE r_alpha_d;
extern __constant__ FLOAT_TYPE b_alpha_d;
extern __constant__ FLOAT_TYPE bubble_radius_d;
extern __constant__ FLOAT_TYPE r_density_d;
extern __constant__ FLOAT_TYPE b_density_d;
extern __constant__ bool external_force_d;
extern __constant__ FLOAT_TYPE A_d;
extern __constant__ FLOAT_TYPE r_viscosity_d;
extern __constant__ FLOAT_TYPE b_viscosity_d;

//#### 2D Color Gradient ####//
extern __constant__ FLOAT_TYPE control_param_d;
extern __constant__ FLOAT_TYPE phi_d[9];
extern __constant__ FLOAT_TYPE teta_d[9];
extern __constant__ FLOAT_TYPE chi_d[9];
extern __constant__ FLOAT_TYPE psi_d[9];
extern __constant__ FLOAT_TYPE w_pert_d[9];
extern __constant__ FLOAT_TYPE c_norms_d[9];
extern __constant__ FLOAT_TYPE cg_w_d[9];
extern __constant__ FLOAT_TYPE hocg_w_d[25];
extern __constant__ int hocg_cx_d[25];
extern __constant__ int hocg_cy_d[25];

//#### 3D Color Gradient ####//
extern __constant__ FLOAT_TYPE c_norms3D_d[19];
extern __constant__ FLOAT_TYPE w_pert3D_d[19];
extern __constant__ FLOAT_TYPE phi3D_d[19];
extern __constant__ FLOAT_TYPE teta3D_d[19];
extern __constant__ FLOAT_TYPE chi3D_d[19];
extern __constant__ FLOAT_TYPE psi3D_d[19];
extern __constant__ FLOAT_TYPE cg_w3D_d[19];
extern __constant__ FLOAT_TYPE hocg_w3D_d[105];
extern __constant__ int hocg_cx3D_d[105];
extern __constant__ int hocg_cy3D_d[105];
extern __constant__ int hocg_cz3D_d[105];
extern __constant__ int hoc3D_d[105];

//#### 3D d3q19 ####//
extern __constant__ int cx3D_d[19];                   ///< velocity x unit vector components
extern __constant__ int cy3D_d[19];                   ///< velocity y unit vector components
extern __constant__ int cz3D_d[19];                   ///< velocity y unit vector components
extern __constant__ int c3D_d[19];                    ///< direction offset levels
extern __constant__ int opp3D_d[19];                  ///< opposite lattice offset
extern __constant__ FLOAT_TYPE w3D_d[19];             ///< lattice weights
extern __constant__ FLOAT_TYPE velMomMap3D_d[361];    ///< MRT constants: mapping between velocity and momentum space \f$ \mathbf{M} \f$
extern __constant__ FLOAT_TYPE momCollMtx3D_d[361];   ///< MRT constants: collision matrix in momentum space \f$ \mathbf{M}^{-1}\mathbf{S} \f$

#endif

#endif
