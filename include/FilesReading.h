#ifndef FILESREADING_H
#define FILESREADING_H

#include "FloatType.h"
#include "Arguments.h"
/**
 * @author Istvan Tamas Jozsa (jozsait@gmail.com) - 2D version
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - development of the existing functions into 3D version
 */
 
/**
 * Get the number of lines in a file
 * @param[in] filename filename
 * @return number of lines
 */
int getNumberOfLines(const char *filename);

/**
 * Read init file (usually SetUpData.ini)
 * @param[in] filename filename
 * @param[out] args input parameters
 */
void readInitFile(const char* filename, Arguments *args);
int readInitConditionsFile(const char* filename, int NoOfNodes, int n, int m, int h, FLOAT_TYPE **u_init, FLOAT_TYPE **v_init, FLOAT_TYPE **w_init, FLOAT_TYPE **rho_init);

/**
 * Read node file (usually D3node.dat)
 * Arrays supplied are allocated inside, freeing these arrays is the caller's responsibility
 *
 * @param[in] filename filename
 * @param[out] ni,nj,nz node ID (x,y,z)
 * @param[out] nx,ny,nz node coordinate (x,y,z)
 * @param[out] nf node type (0: solid, 1:fluid)
 * @return number of lines read
 */
int readNodeFile(const char *filename, int **ni, int **nj, int **nk, FLOAT_TYPE **nx, FLOAT_TYPE **ny, FLOAT_TYPE **nz, int **nf,int problemtype);


/**
 * Read boundary conditions file (usually BCconnectors.dat)
 * Arrays supplied are allocated inside, freeing these arrays is the caller's responsibility
 *
 *
 * @param[in] filename filename
 * @param[out] ni,nj,nz node ID (x,y,z)
 * @param[out] dir direction
 * @param[out] bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param[out] bcx,bcy,bcz BC coordinate (x,y,z)
 * @param[out] id boundary ID (for drag/lift computation)
 * @return number of lines read
 */
int readConnFile(const char *filename, int **ni, int **nj, int **nk, int **dir, int **bc,
                 FLOAT_TYPE **bcx, FLOAT_TYPE **bcy, FLOAT_TYPE **bcz, int **id, int problemtype);
/**
 * Read result file (usually FinalData.csv)
 * @param[in] filename filename
 * @param[out] results coordinates (x,y), velocity (u,v), vel_mag, rho, pressure
 * @param[out] fluid node type (0: solid, 1:fluid)
 * @return number of lines read
 */
int readResultFile(const char *filename, FLOAT_TYPE ***results);

/**
 * @brief Compare result files
 * @note Use it to compare FinalData.csv with previous result.
 *
 * @param f1,f2 files to compare
 * @return 0 for identical files and 1 for different files
 */
__host__ int compareFiles(const char* f1, const char* f2);

/**
 * Gets the last value of the array
 * Use it with nodeIdx to get number rows, and with nodeIdy to get number of columns
 * @param arr input array
 * @param n array size
 * @return last element plus one
 */
int getLastValue(int *arr, int n);

/**
 * Returns the grid spacing of the mesh
 * @param ni,nj node ID (x,y)
 * @param nx node coordinate x
 * @param n number of nodes
 * @return grid spacing
 */
FLOAT_TYPE getGridSpacing(int *ni, int *nj, FLOAT_TYPE *nx, int n);

/**
 * Get number of inlet nodes
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param n number of BC nodes
 * @return number of inlet nodes
 */
int getNumInletNodes(int *bc, int *dir, int n, int problemtype);

/**
 * Get the maximum coordinate of the inlet in the Y direction
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param bcy BC coordinate y
 * @param delta grid spacing
 * @param n number of BC nodes
 * @return maximum coordinate of the inlet in the Y direction
 */
FLOAT_TYPE getMaxInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy, FLOAT_TYPE delta, int n,int problemtype);

/**
 * Get the maximum coordinate of the inlet in the Z direction
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param bcy BC coordinate z
 * @param delta grid spacing
 * @param n number of BC nodes
 * @return maximum coordinate of the inlet in the Z direction
 */
FLOAT_TYPE getMaxInletCoordZ(int *bc, int *dir, FLOAT_TYPE *bcz, FLOAT_TYPE delta, int n,int problemtype);

/**
 * Get the minimum coordinate of the inlet in the Y direction
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param bcy BC coordinate y
 * @param delta grid spacing
 * @param n number of BC nodes
 * @return minimum coordinate of the inlet in the Y direction
 */
FLOAT_TYPE getMinInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy, FLOAT_TYPE delta, int n,int problemtype);

/**
 * Get the minimum coordinate of the inlet in the Z direction
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param bcz BC coordinate z
 * @param delta grid spacing
 * @param n number of BC nodes
 * @return minimum coordinate of the inlet in the Z direction
 */
FLOAT_TYPE getMinInletCoordZ(int *bc, int *dir, FLOAT_TYPE *bcz, FLOAT_TYPE delta, int n,int problemtype);

int readArray(const char *file,FLOAT_TYPE **array);
int compareTestFiles(const char* f1, const char* f2);
#endif
