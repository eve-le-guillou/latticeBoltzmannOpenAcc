/**
 * This file contains the acceptable arguments for the solver.
 * @file Arguments.h
 * @author Adam Koleszar (adam.koleszar@gmail.com) - 2D 
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - development of new features used in 3D solver
 */

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "FloatType.h"

/// Inlet profile options
typedef enum
{
    INLET=1,        ///< inlet profile
    NO_INLET,       ///< no inlet profile
    PULSATILE_INLET ///< pulsatile inlet profile @warning not implemented
} InletProfile;

/// Outlet profile options
typedef enum
{
    OUTLET=1,      ///< outlet profile
    OUTLET_SECOND, ///< open boundary second order
    OUTLET_FIRST, ///< open boundary first order
    OUTLET_HEmodel
} OutletProfile;

/// Collision models
typedef enum
{
    BGKW=1, ///< BGKW method
    TRT,    ///< TRT method
    MRT     ///< MRT method
} CollisionModel;

/// Boundary types
typedef enum
{
    CURVED=1, ///< curved boundaries
    STRAIGHT  ///< straight boundaries
} BoundaryType;

typedef enum
{
  _2D=1,
  _3D
} ProblemType;

/// Boundary Wall model
typedef enum
{
    SIMPLE=1,
    COMPLEX
} BCWallModel;

/// Output formats
typedef enum
{
    CSV=1, ///< paraview format (.csv)
    TECPLOT,     ///< tecplot format (.dat)
    PARAVIEW	///<paraview format (vti.)
} OutputFormat;

/// Result from command line argument parsing
typedef enum
{
    NORMAL=0, ///< normal execution, parameters read
    HELP,     ///< print help
    INIT,     ///< read values from init file
    TEST,     ///< run unit tests
    COMPARE,  ///< compare results
    ERROR     ///< error happened
} ArgResult;
typedef enum
{
	L2=1, ///< normal execution, parameters read
    FdRelDiff,     ///< print help
    MacroDiff

} ResidualsModel;


/// Input parameters
typedef struct arg
{
	char problemName[512];		   ///< problem files name
    FLOAT_TYPE u;                  ///< velocity x component
    FLOAT_TYPE v;                  ///< velocity y component
    FLOAT_TYPE w;                  ///< velocity z component
    FLOAT_TYPE rho;                ///< density
    FLOAT_TYPE viscosity;          ///< viscosity
    FLOAT_TYPE g;				   ///< force
    InletProfile inletProfile;     ///< inlet profile
    OutletProfile outletProfile;   ///< outlet profile
    CollisionModel collisionModel; ///< collision model
    BoundaryType boundaryType;     ///< boundary type
    BCWallModel bcwallmodel;       ///< boundary type
    OutputFormat outputFormat;     ///< output format
    int iterations;                ///< number of iterations
    FLOAT_TYPE StopCondition[4];	///< Stop condition to stop simulator
    int autosaveEvery;             ///< autosave every n iteration
    int autosaveAfter;             ///< autosave after nth iteration
    int boundaryId;                ///< boundary ID for drag/lift calculation
    char initialConditionsName[512]; ///< filename with initialconditions
    int UseInitialCondFromFile;
    int TypeOfResiduals;
    int ShowMacroDiff;
    int TypeOfProblem;
    int UpdateInltOutl;
    //New
    int multiPhase;
    FLOAT_TYPE r_density;
    FLOAT_TYPE b_density;
    FLOAT_TYPE r_alpha;
    FLOAT_TYPE b_alpha;
    FLOAT_TYPE gamma;
    FLOAT_TYPE r_viscosity;
    FLOAT_TYPE b_viscosity;
    FLOAT_TYPE A;
    FLOAT_TYPE beta;
    FLOAT_TYPE control_param;
    FLOAT_TYPE g_limit;
    int test_case;
    bool external_force;
    bool high_order;
    bool enhanced_distrib;

    //Bubble case
    FLOAT_TYPE bubble_radius;
    //COUETTE case
    FLOAT_TYPE kappa;

} Arguments;

/// Input filenames
typedef struct inputf
{
    char init[512];   ///< init file name
    char node[512];   ///< node file name
    char bc[512];     ///< bc file name
    char result[512]; ///< results directory name
    char comp[512];   ///< filename to compare
    char final[512];  ///< filename to compare to
    char InitialConditions[512];  ///< initial conditions filename
} InputFilenames;

/**
 * Parse command line arguments
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param inFn input filenames
 * @param args input parameters
 * @return task that need to be done, or error
 *
 * @test MANUAL:
 * In order to test argument handling try some of these cases
 *   - run with -h, should display help
 *   - run without arguments, should work as written in SetUpData.ini
 *   - run with -f \<ini\> and any other argument, should work as written in SetUpData.ini
 *   - run with -o \<folder\> should output results to given folder
 *   - run with -i \<n\> should do iterations n times
 *   - run with -n \<file\> with non-existent file, should print out error message
 *   - try using long argument options as well, should work as the short version
 *   - run in compare mode with one file given, should compare Results/FinalData.csv to given file
 *   - run in compare mode with two files, should compare the 2 given files
 */
ArgResult handleArguments(int argc, char **argv, InputFilenames *inFn, Arguments *args);

#endif
