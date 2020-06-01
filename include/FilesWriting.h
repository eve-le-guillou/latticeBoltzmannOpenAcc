/**
 * Functions for file writing
 * @file FilesWriting.h
 * @author Istvan Tamas Jozsa (jozsait@gmail.com) - 2D version
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - development of the existing functions into 3D version
 */

#ifndef FILESWRITING_H
#define FILESWRITING_H

#include "FloatType.h"
#include "Arguments.h"

/**
 * @brief Print results to output file
 *
 * @param OutputFile name of the outpu file
 * @param CoordX,CoordY x,y coordinates
 * @param U,V velocity x,y coordinates
 * @param Rho density
 * @param Fluid node type (0: solid, 1: fluid)
 * @param n number of rows
 * @param m number of columns
 * @param outputFormat output file format
 */
void WriteResults3D(char* OutputFile, int* fluid_d, FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* CoordZ,
					FLOAT_TYPE* U, FLOAT_TYPE* V, FLOAT_TYPE* W, FLOAT_TYPE* Rho, int* Fluid,
					int n, int m, int h, OutputFormat outputFormat);
void WriteLidDrivenCavityMidLines3D(FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* CoordZ,
					FLOAT_TYPE* U, FLOAT_TYPE* W, int n, int m, int h, FLOAT_TYPE Ulid);
void WriteChannelCrossSection3D(FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* CoordZ,
					FLOAT_TYPE* U, FLOAT_TYPE* V, FLOAT_TYPE* W, int n, int m, int h, FLOAT_TYPE Ulid);
void eraseAlertLog();

void WriteResultsMultiPhase(char* OutputFile, int* fluid, FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY,
		FLOAT_TYPE* CoordZ,
		FLOAT_TYPE* U, FLOAT_TYPE* V, FLOAT_TYPE* W, FLOAT_TYPE* Rho,FLOAT_TYPE* r_Rho, FLOAT_TYPE* b_Rho,
		int* Fluid, int n, int m, int h, OutputFormat outputFormat);

void WriteArray(char* fileName, FLOAT_TYPE *arr, int n, int m, int h = 1, int dir = 1);

void writeCouetteSolution(char* fileName, FLOAT_TYPE *analytical, FLOAT_TYPE *computed, FLOAT_TYPE *y, int m, int n, int h = 0);

void writeOscilatingSolution(char* fileName, FLOAT_TYPE *extremes, int size);
#endif
