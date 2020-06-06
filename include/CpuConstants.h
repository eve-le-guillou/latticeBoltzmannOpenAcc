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

//extern __constant__ InletProfile inletProfile_d;     ///< inlet profile
//extern __constant__ BoundaryType boundaryType_d;     ///< boundary type
//extern __constant__ OutletProfile outletProfile_d;   ///< outlet profile
extern FLOAT_TYPE rhoIn;              ///< input density
extern FLOAT_TYPE uIn;                ///< input velocity x
extern FLOAT_TYPE vIn;                ///< input velocity y
//extern int dlBoundaryId_d;              ///< boundary ID
extern int depth;                     ///< number of rows Y
extern int length;                    ///< number of columns X
extern int height;                    ///< number of layers Z
extern FLOAT_TYPE delta;              ///< grid spacing
extern FLOAT_TYPE minInletCoordY;     ///< maximum inlet coordinate y
extern FLOAT_TYPE maxInletCoordY;     ///< minimum inlet coordinate y
extern FLOAT_TYPE minInletCoordZ;     ///< maximum inlet coordinate z
extern FLOAT_TYPE maxInletCoordZ;     ///< minimum inlet coordinate z
extern FLOAT_TYPE omega;            ///< collision frequency for D2Q9 \f$ \omega = \frac{1}{3\nu + 0.5} \f$
extern FLOAT_TYPE omegaA;           ///< asymmetric collision frequency \f$ \omega_a = \frac{8(2-\omega)}{8-\omega} \f$
extern FLOAT_TYPE g;
//#### 2D d2q9 ####//
extern int cx2D[9];                   ///< velocity x unit vector components
extern int cy2D[9];                   ///< velocity y unit vector components
extern int c2D[9];                    ///< direction offset levels
extern int opp2D[9];                  ///< opposite lattice offset
extern FLOAT_TYPE w2D[9];             ///< lattice weights
extern FLOAT_TYPE velMomMap2D[81];    ///< MRT constants: mapping between velocity and momentum space \f$ \mathbf{M} \f$
extern FLOAT_TYPE momCollMtx2D[81];   ///< MRT constants: collision matrix in momentum space \f$ \mathbf{M}^{-1}\mathbf{S} \f$

//#### Color Gradient ####//
extern FLOAT_TYPE beta;
extern FLOAT_TYPE g_limit;
extern FLOAT_TYPE r_alpha;
extern FLOAT_TYPE b_alpha;
extern FLOAT_TYPE bubble_radius;
extern FLOAT_TYPE r_density;
extern FLOAT_TYPE b_density;
extern bool external_force;
extern FLOAT_TYPE A;
extern FLOAT_TYPE r_viscosity;
extern FLOAT_TYPE b_viscosity;

//#### 2D Color Gradient ####//
extern FLOAT_TYPE control_param;
extern FLOAT_TYPE phi[9];
extern FLOAT_TYPE teta[9];
extern FLOAT_TYPE chi[9];
extern FLOAT_TYPE psi[9];
extern FLOAT_TYPE w_pert[9];
extern FLOAT_TYPE c_norms[9];
extern FLOAT_TYPE cg_w[9];
extern FLOAT_TYPE hocg_w[25];
extern int hocg_cx[25];
extern int hocg_cy[25];

//#### 3D Color Gradient ####//
extern FLOAT_TYPE c_norms3D[19];
extern FLOAT_TYPE w_pert3D[19];
extern FLOAT_TYPE phi3D[19];
extern FLOAT_TYPE teta3D[19];
extern FLOAT_TYPE chi3D[19];
extern FLOAT_TYPE psi3D[19];
extern FLOAT_TYPE cg_w3D[19];
extern FLOAT_TYPE hocg_w3D[105];
extern int hocg_cx3D[105];
extern int hocg_cy3D[105];
extern int hocg_cz3D[105];
extern int hoc3D[105];

//#### 3D d3q19 ####//
extern int cx3D[19];                   ///< velocity x unit vector components
extern int cy3D[19];                   ///< velocity y unit vector components
extern int cz3D[19];                   ///< velocity y unit vector components
extern int c3D[19];                    ///< direction offset levels
extern int opp3D[19];                  ///< opposite lattice offset
extern FLOAT_TYPE w3D[19];             ///< lattice weights
extern FLOAT_TYPE velMomMap3D[361];    ///< MRT constants: mapping between velocity and momentum space \f$ \mathbf{M} \f$
extern FLOAT_TYPE momCollMtx3D[361];   ///< MRT constants: collision matrix in momentum space \f$ \mathbf{M}^{-1}\mathbf{S} \f$

#endif

#endif
