/**
 * File reading functions
 * @file FilesReading.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "GpuSum.h"
#include "FilesReading.h"
#include "ArrayUtils.h"
#include "math.h"
#include<string>
#include <iostream>
#include<fstream>
#include<stdexcept>
int getNumberOfLines(const char *filename) {
	FILE *f = fopen(filename, "r");
	if (!f) {
		fprintf(stderr, "Error reading file %s: %s\n", filename,
				strerror(errno));
		return 0;
	}

	int lines = 0;
	char c;
	while (!feof(f)) {
		c = fgetc(f);
		if (c == '\n' || c == '\r') {
			lines++;
		}
	}
	fclose(f);

	return lines;
}

/**
 * @brief Set enum options from init file
 *
 * @param args input parameters
 * @param ipr inlet profile
 * @param coll collision model
 * @param curved boundary type
 * @param opr outlet profile
 * @param format output format
 */
void setOptions(Arguments *args, int ipr, int coll, int curved, int opr,
		int format, int bcwalltype, int residualstype) {
	args->inletProfile = (InletProfile) ipr;
	args->collisionModel = (CollisionModel) coll;
	args->boundaryType = (BoundaryType) curved;
	args->bcwallmodel = (BCWallModel) bcwalltype;
	args->outletProfile = (OutletProfile) opr;
	args->outputFormat = (OutputFormat) format;
	args->TypeOfResiduals = (ResidualsModel) residualstype;

}

void readInitFile(const char* filename, Arguments *args) {
	int ipr, coll, curved, opr, format, bcwalltype, ResidualsType;
	std::ifstream initfile;
	initfile.open(filename);
	std::string line[23];
	for (int i = 0; i < 23; i++) {
		std::getline(initfile, line[i]);
		std::cout << line[i] << std::endl;
		if (line[i].find("= ") != std::string::npos)
			line[i] = line[i].substr(line[i].find("= ") + 2);
		else
			fprintf(stderr, "Error reading initial file\n");
	}
	initfile.close();
	if (line[0].find("\"") != std::string::npos)
		line[0] = line[0].substr(line[0].find("\"") + 1);
	else
		fprintf(stderr, "Error reading initial file\n");
	if (line[0].find("\"") != std::string::npos)
		line[0] = line[0].substr(0, line[0].find("\""));
	else
		fprintf(stderr, "Error reading initial file\n");
	sprintf(args->problemName, "%s", line[0].c_str());

	if (line[17].find("\"") != std::string::npos)
		line[17] = line[17].substr(line[17].find("\"") + 1);
	else
		fprintf(stderr, "Error reading initial file\n");
	if (line[17].find("\"") != std::string::npos)
		line[17] = line[17].substr(0, line[17].find("\""));
	else
		fprintf(stderr, "Error reading initial file\n");
	sprintf(args->initialConditionsName, "%s", line[17].c_str());
	args->u = atof(line[1].c_str());
	args->v = atof(line[2].c_str());
	args->w = atof(line[3].c_str());
	args->rho = atof(line[4].c_str());
	args->viscosity = atof(line[5].c_str());
	args->iterations = atoi(line[11].c_str());
	args->autosaveAfter = atoi(line[12].c_str());
	args->autosaveEvery = atoi(line[13].c_str());
	sscanf(line[14].c_str(), FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT, &args->StopCondition[0],
			&args->StopCondition[1], &args->StopCondition[2],
			&args->StopCondition[3]);

	args->ShowMacroDiff = atoi(line[19].c_str());
	args->g = atof(line[20].c_str());
	args->boundaryId = atof(line[21].c_str());
	if (line[15].find("CSV") != std::string::npos)
		format = 1;
	else {
		if (line[15].find("DAT") != std::string::npos)
			format = 2;
		else {
			if (line[15].find("VTI") != std::string::npos)
				format = 3;
			else {
				fprintf(stderr, "Error reading initial file\n");
			}
		}
	}
	if (line[16].find("yes") != std::string::npos)
		args->UseInitialCondFromFile = 1;
	else {
		if (line[16].find("no") != std::string::npos)
			args->UseInitialCondFromFile = 0;
		else {
			fprintf(stderr, "Error reading initial file\n");
		}
	}
	if (line[6].find("INLET_PROF") != std::string::npos)
		ipr = 1;
	else {
		if (line[6].find("INLET_CTE") != std::string::npos)
			ipr = 2;
		else {
			if (line[6].find("INLET_PUL") != std::string::npos)
				ipr = 3;
			else {
				fprintf(stderr, "Error reading initial file\n");
			}
		}
	}
	if (line[7].find("SRT") != std::string::npos) {
		coll = 1;

	} else {
		if (line[7].find("TRT") != std::string::npos) {
			coll = 2;
		}

		else {
			if (line[7].find("MRT") != std::string::npos) {
				coll = 3;
			}

			else {
				fprintf(stderr, "Error reading initial file\n");
			}
		}
	}
	if (line[8].find("STRAIGHT") != std::string::npos)
		curved = 2;
	else {
		if (line[8].find("CURVED") != std::string::npos)
			curved = 1;
		else {
			fprintf(stderr, "Error reading initial file\n");
		}
	}
	if (line[9].find("BBack") != std::string::npos) {
		printf("here\n");
		bcwalltype = 1;
	} else {
		if (line[9].find("HHmodel") != std::string::npos)
			bcwalltype = 2;
		else {
			fprintf(stderr, "Error reading initial file\n");
		}
	}
	if (line[10].find("Vout") != std::string::npos)
		opr = 1;
	else {
		if (line[10].find("1stE") != std::string::npos)
			opr = 2;
		else {
			if (line[10].find("2ndE") != std::string::npos)
				opr = 3;
			else {
				if (line[10].find("HEmodel") != std::string::npos)
					opr = 4;
				else {
					fprintf(stderr, "Error reading initial file\n");
				}
			}
		}
	}
	if (line[18].find("L2") != std::string::npos)
		ResidualsType = 1;
	else {
		if (line[18].find("MaxRelativeFdDiff") != std::string::npos)
			ResidualsType = 2;
		else {
			if (line[18].find("MaxMacroDiff") != std::string::npos)
				ResidualsType = 3;

			else {
				fprintf(stderr, "Error reading initial file\n");
			}
		}
	}
	if (line[22].find("yes") != std::string::npos)
		args->UpdateInltOutl = 1;
	else {
		if (line[22].find("no") != std::string::npos)
			args->UpdateInltOutl = 0;
		else {
			fprintf(stderr, "Error reading initial file\n");
		}
	}

	setOptions(args, ipr, coll, curved, opr, format, bcwalltype, ResidualsType);
}

int readInitConditionsFile(const char* filename, int NoOfNodes, int n, int m,
		int h,
		FLOAT_TYPE **u_init, FLOAT_TYPE **v_init, FLOAT_TYPE **w_init,
		FLOAT_TYPE **rho_init) {
	*u_init = createHostArrayFlt(m * n * h);
	*v_init = createHostArrayFlt(m * n * h);
	*w_init = createHostArrayFlt(m * n * h);
	*rho_init = createHostArrayFlt(m * n * h);

	int NoOfLines = getNumberOfLines(filename);
	FLOAT_TYPE **res1;

	if (!NoOfLines) {
		printf("empty/no initial condition file!");
		return -1;
	}
	char * line = NULL;
	size_t len = 0;
	FILE *f = fopen(filename, "r");
	getline(&line, &len, f);
	if (NoOfLines - 3 == NoOfNodes && line[0] == 'T' && line[1] == 'i') {
		//tecplot format
	} else {
		if (line[0] == '<' && line[1] == '?'
				&& NoOfLines == (18 + 3 * h * m * n)) {
			//paraview format
		} else {
			if (NoOfLines - 1 == NoOfNodes && line[0] == 'x'
					&& line[1] == ',') { //csv format
			} else {
				printf(
						"wrong file! Number of nodes in initial condition file is incorrect\n");
				return -1;
			}
		}
	}
	fclose(f);
	readResultFile(filename, &res1);
	*u_init = res1[3];
	*v_init = res1[4];
	*w_init = res1[5];
	*rho_init = res1[7];
	return 0;

}

int readNodeFile(const char *filename, int **ni, int **nj, int **nk,
		FLOAT_TYPE **nx, FLOAT_TYPE **ny, FLOAT_TYPE **nz, int **nf, int problemtype) {
	int n = getNumberOfLines(filename);
	if (!n) {
		return 0;
	}

	*ni = createHostArrayInt(n);
	*nj = createHostArrayInt(n);
	*nk = createHostArrayInt(n);
	*nx = createHostArrayFlt(n);
	*ny = createHostArrayFlt(n);
	*nz = createHostArrayFlt(n);
	*nf = createHostArrayInt(n);

	FILE *f = fopen(filename, "r");
	int i;
	for (i = 0; i < n; ++i) {
		switch (problemtype) {
		case _2D:
			fscanf(f, "%d %d "FLOAT_FORMAT" "FLOAT_FORMAT" %d", (*ni) + i,
					(*nj) + i, (*nx) + i, (*ny) + i, (*nf) + i);
			break;
		case _3D:
			fscanf(f,
					"%d\t%d\t%d\t"FLOAT_FORMAT"\t"FLOAT_FORMAT"\t"FLOAT_FORMAT"\t%d",
					(*ni) + i, (*nj) + i, (*nk) + i, (*nx) + i, (*ny) + i,
					(*nz) + i, (*nf) + i);
			break;

		}

	}
	fclose(f);
	return n;
}

int readConnFile(const char *filename, int **ni, int **nj, int **nk, int **dir,
		int **bc,
		FLOAT_TYPE **bcx, FLOAT_TYPE **bcy, FLOAT_TYPE **bcz, int **id,
		int problemtype) {
	int n = getNumberOfLines(filename);
	if (!n) {
		return 0;
	}

	*ni = createHostArrayInt(n);
	*nj = createHostArrayInt(n);
	*nk = createHostArrayInt(n);
	*dir = createHostArrayInt(n);
	*bc = createHostArrayInt(n);
	*bcx = createHostArrayFlt(n);
	*bcy = createHostArrayFlt(n);
	*bcz = createHostArrayFlt(n);
	*id = createHostArrayInt(n);

	FILE *f = fopen(filename, "r");
	int i;
	for (i = 0; i < n; ++i) {
		switch (problemtype) {
		case _2D:
			fscanf(f, "%d %d %d %d "FLOAT_FORMAT" "FLOAT_FORMAT" %d", (*ni) + i,
					(*nj) + i, (*dir) + i, (*bc) + i, (*bcx) + i, (*bcy) + i,
					(*id) + i);
			break;
		case _3D:
			fscanf(f,
					"%d\t%d\t%d\t\t%d\t%d\t"FLOAT_FORMAT"\t"FLOAT_FORMAT"\t"FLOAT_FORMAT"\t%d",
					(*ni) + i, (*nj) + i, (*nk) + i, (*dir) + i, (*bc) + i,
					(*bcx) + i, (*bcy) + i, (*bcz) + i, (*id) + i);
			break;

		}

	}
	fclose(f);
	return n;
}
int readResultFile(const char *filename, FLOAT_TYPE ***results) {
	int i, n = getNumberOfLines(filename);
	if (!n) {
		return 0;
	}

	FILE *f = fopen(filename, "r");
	*results = create2DHostArrayFlt(n, 9);
	char *line = NULL;
	size_t len;
	getline(&line, &len, f);
	if (line[0] == 'x' && line[1] == ',') { //VTI format

		n -= 1;
		for (i = 0; i < n; ++i) {
			fscanf(f,
					FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT"\n",
					results[0][0] + i, results[0][1] + i, results[0][2] + i,
					results[0][3] + i, results[0][4] + i, results[0][5] + i,
					results[0][6] + i, results[0][7] + i, results[0][8] + i);
		}
	} else {

		if (line[0] == 'T' && line[1] == 'i') { 	//tecplot format
			getline(&line, &len, f);
			getline(&line, &len, f);
			for (i = 0; i < n - 3; ++i) {

				fscanf(f,
						FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT"\n",
						results[0][0] + i, results[0][1] + i, results[0][2] + i,
						results[0][3] + i, results[0][4] + i, results[0][5] + i,
						results[0][6] + i, results[0][7] + i,
						results[0][8] + i);
			}
		} else {
			if (line[0] == '<' && line[1] == '?') { //vti format
				for (i = 0; i < 6; ++i)
					getline(&line, &len, f);
				for (i = 0; i < (n - 18) / 3; ++i) {

					fscanf(f,
							FLOAT_FORMAT" "FLOAT_FORMAT" "FLOAT_FORMAT"\n",
							results[0][3] + i, results[0][4] + i,
							results[0][5] + i);
				}
				getline(&line, &len, f);
				getline(&line, &len, f);
				for (i = 0; i < (n - 18) / 3; ++i) {

					fscanf(f,
							FLOAT_FORMAT"\n", results[0][7] + i);
				}
			}
		}
	}

	fclose(f);
	return n;
}

__host__ int compareFiles(const char* f1, const char* f2) {
	printf("Comparing results\n");
	int l1 = getNumberOfLines(f1);
	int l2 = getNumberOfLines(f2);
	if (l1 != l2) {
		printf("Line number mismatch\n");
		exit(1);
	}
	int i;

	const char *columnName[] = { "x", "y", "z", "u", "v", "w", "vel_mag", "rho",
			"pressure" };

	FLOAT_TYPE **res1;
	FLOAT_TYPE **res2;

	printf("Reading files...");
	readResultFile(f1, &res1);
	printf("...");
	readResultFile(f2, &res2);
	printf("...done\nComparing results...\n");

	dim3 gridDim(l1 / THREADS + 1);
	FLOAT_TYPE *da = createGpuArrayFlt(l1);
	FLOAT_TYPE *db = createGpuArrayFlt(l1);
	FLOAT_TYPE *dc = createGpuArrayFlt(l1);
	FLOAT_TYPE *dd = createGpuArrayFlt(l1);
	FLOAT_TYPE diff_sum[9];
	FLOAT_TYPE diff_max[9];
	bool divergence = false;
	for (i = 0; i < 9; ++i) {

		cudaMemcpy(da, res1[i], SIZEFLT(l1), cudaMemcpyHostToDevice);
		cudaMemcpy(db, res2[i], SIZEFLT(l1), cudaMemcpyHostToDevice);

		gpu_abs_sub<<<gridDim, THREADS>>>(da, db, dc, l1, &divergence);
		diff_sum[i] = gpu_sum_h(dc, dd, l1);
		gpu_abs_sub<<<gridDim, THREADS>>>(da, db, dc, l1, &divergence);
		diff_max[i] = gpu_max_h(dc, dd, l1);

	}

	printf("     array |        diff sum |        diff max |      diff/nodes\n"
			"----------------------------------------------\n");

	int b = 0;
	for (i = 0; i < 9; ++i) {
		printf("%10s | %15g | %15g | %15g\n", columnName[i], diff_sum[i],
				diff_max[i], diff_sum[i] / l1);
		b |= diff_sum[i] > 0.001;
		free(res1[i]);
		free(res2[i]);
	}

	free(res1);
	free(res2);

	return b;
}

int getLastValue(int *arr, int n) {
	return arr[n - 1] + 1;
}

FLOAT_TYPE getGridSpacing(int *ni, int *nj, FLOAT_TYPE *nx, int n) {
	int i;
	FLOAT_TYPE delta1, delta2;
	for (i = 0; i < n; ++i) {
		if (ni[i] == 0 && nj[i] == 0) {
			delta1 = nx[i];
			break;
		}
	}
	for (i = 0; i < n; ++i) {
		if (ni[i] == 1 && nj[i] == 0) {
			delta2 = nx[i];
			break;
		}
	}

	return (FLOAT_TYPE) fabs(delta2 - delta1);
}

int getNumInletNodes(int *bc, int *dir, int n, int problemtype) {
	int i;
	int nodes = 0;
	for (i = 0; i < n; ++i) {
		switch (problemtype) {
		case _2D:
			if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 4) {
				++nodes;
			}
			break;
		case _3D:
			if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 6) {
				++nodes;
			}
			break;
		}
	}
	return nodes;
}

FLOAT_TYPE getMaxInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy,
		FLOAT_TYPE delta, int n, int problemtype) {
	int i = 0;
	FLOAT_TYPE maxY = 0.0;
	while (bc[i] != 2 && i < n) //inlet
	{
		maxY = bcy[++i];
	}
	for (i = 0; i < n; ++i) {
		switch (problemtype) {
		case _2D:
			if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 4) {
				maxY = (bcy[i] > maxY) ? bcy[i] : maxY;
			}
			break;
		case _3D:
			if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 6) {
				maxY = (bcy[i] > maxY) ? bcy[i] : maxY;
			}
			break;
		}

	}

	return maxY + delta / 2;
}

FLOAT_TYPE getMaxInletCoordZ(int *bc, int *dir, FLOAT_TYPE *bcz,
		FLOAT_TYPE delta, int n, int problemtype) {
	if (problemtype == (ProblemType) _2D)
		return 0.0;
	int i = 0;
	FLOAT_TYPE maxZ;
	while (bc[i] != 2 && i < n) //inlet
	{
		maxZ = bcz[++i];
	}
	for (i = 0; i < n; ++i) {

		if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 6) {
			maxZ = (bcz[i] > maxZ) ? bcz[i] : maxZ;
		}
	}

	return maxZ + delta / 2;
}

FLOAT_TYPE getMinInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy,
		FLOAT_TYPE delta, int n, int problemtype) {
	int i = 0;
	FLOAT_TYPE minY= 0.0;
	while (bc[i] != 2 && i < n) //inlet
	{
		minY = bcy[++i];
	}

	for (i = 0; i < n; ++i) {
		switch (problemtype) {
		case _2D:
			if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 4) {
				minY = (bcy[i] < minY) ? bcy[i] : minY;
			}
			break;
		case _3D:
			if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 6) {
				minY = (bcy[i] < minY) ? bcy[i] : minY;
			}
			break;
		}

	}

	return minY - delta / 2;
}

FLOAT_TYPE getMinInletCoordZ(int *bc, int *dir, FLOAT_TYPE *bcz,
		FLOAT_TYPE delta, int n, int problemtype) {
	if (problemtype == (ProblemType) _2D)
		return 0.0;
	int i = 0;
	FLOAT_TYPE minZ;
	while (bc[i] != 2 && i < n) //inlet
	{
		minZ = bcz[++i];
	}

	for (i = 0; i < n; ++i) {
		if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 6) {
			minZ = (bcz[i] < minZ) ? bcz[i] : minZ;
		}
	}

	return minZ - delta / 2;
}

int readArray(const char *file,FLOAT_TYPE **array){

	int i, n = getNumberOfLines(file);
	if (!n) {
		return 0;
	}

	FILE *f = fopen(file, "r");
	*array = createHostArrayFlt(n,ARRAY_NONE);

	for(i = 0; i < n; i++){
		fscanf(f, FLOAT_FORMAT"\n", (array[0] + i));
	}
	fclose(f);
	return 1;
}

int compareTestFiles(const char* f1, const char* f2){
	printf("Comparing results\n");
	int l1 = getNumberOfLines(f1);
	int l2 = getNumberOfLines(f2);
	if (l1 != l2) {
		printf("Line number mismatch %d vs %d\n", l1,l2);
		exit(1);
	}

	FLOAT_TYPE *arr1;
	FLOAT_TYPE *arr2;

	if(!readArray(f1, &arr1)){
		printf("Error reading file\n");
		exit(1);
	}

	if(!readArray(f2, &arr2)){
		printf("Error reading file\n");
		exit(1);
	}

	FLOAT_TYPE max = 0.0, sum = 0.0, diff;
	int max_index = -1;
	for(int i = 0; i < l1; i++){
		diff = abs(arr1[i] - arr2[i]);
		if(diff > max){
			max = diff;
			max_index = i;
		}
		sum += diff;
	}

	printf("\nFiles finished comparing\nMax error: %f at %d\nSum: %f\n", max,max_index,sum);

	return 0;
}
