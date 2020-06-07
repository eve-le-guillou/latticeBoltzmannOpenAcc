/**
 * @author Alfonso Aguilar (a.aguilar-pontes@cranfield.ac.uk) - implementation of the physics
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - software aspects of the implementation
 */

#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <errno.h>
#include "CpuFunctions.h"       // GPU kernels
#include "ShellFunctions.h"     // For convenience
#include "FilesReading.h"       // For reading files
#include "FilesWriting.h"       // For writing files e.g. tecplot
#include "CellFunctions.h"      // For cell modifications
#include "ComputeResiduals.h"   // residuals
#include "LogWriter.h"
#include "Iterate.h"
#include "ArrayUtils.h"
//#include "Check.h"
//#include "cuda.h"
#include "CpuSum.h"
#include "Multiphase.h"


#define CUDA 1

int Iterate3D(InputFilenames *inFn, Arguments *args) {
	// Time measurement: declaration, begin
	clock_t tStart = clock();

	FILE* logFile;               // file for log
	char autosaveFilename[768];  // autosave filename
	char outputFilename[768];    // initial data will be written to this file
	char finalFilename[768];     // final data will be written to this file
	char logFilename[768];       // path of the .log file
	char residualsFilename[768]; // path of the residuals file
	char timeFilename[768];      // path of time measurement file
	bool firstIter = true;
	//bool *d_divergence;
	int AuxMacroDiff = 1;
	FLOAT_TYPE r = -1.0;
	logFilename[0] = '\0';
	residualsFilename[0] = '\0';
	timeFilename[0] = '\0';

	if (strlen(inFn->result)) {
		strcat(logFilename, inFn->result);
		strcat(residualsFilename, inFn->result);
		strcat(timeFilename, inFn->result);
	}
	strcat(logFilename, "lbmsolver.log");
	strcat(residualsFilename, "residuals.dat");
	strcat(timeFilename, "runtimes.dat");

	int autosaveIt = 1; // autosave i variable, will be incremented after every autosave
	int numNodes, numConns; // This will store the number of lines of the read files
	FLOAT_TYPE delta;          // grid spacing
	int n, m, h;                 // number of nodes in the x, y and z directions
	FLOAT_TYPE maxInletCoordY; // maximum inlet coordinate in y
	FLOAT_TYPE minInletCoordY; // minimum inlet coordinate in y
	FLOAT_TYPE maxInletCoordZ; // maximum inlet coordinate in z
	FLOAT_TYPE minInletCoordZ; // minimum inlet coordinate in z
	int numInletNodes;         // number of inlet nodes
	FLOAT_TYPE uMaxDiff = -1, vMaxDiff = -1, wMaxDiff = -1, rhoMaxDiff = -1, fMaxDiff = -1;
	int *nodeIdX, *nodeIdY, *nodeIdZ, *nodeType, *bcNodeIdX, *bcNodeIdY,
	*bcNodeIdZ, *latticeId, *bcType, *bcBoundId;
	FLOAT_TYPE *nodeX, *nodeY, *nodeZ, *bcX, *bcY, *bcZ;

	FLOAT_TYPE taskTime[9];
	int i;
	for (i = 0; i < 9; ++i) {
		taskTime[i] = 0.0;
	}

	clock_t tInstant1, tInstant2; // Time measurement points, universal
	clock_t tIterStart, tIterEnd; // Time measurement points: main loop

	// cuda time measurement variables
	/*cudaEvent_t start, stop;
	float cudatime;
	CHECK(cudaEventCreate(&start));
	CHECK(cudaEventCreate(&stop));*/

	numNodes = readNodeFile(inFn->node, &nodeIdX, &nodeIdY, &nodeIdZ, &nodeX,
			&nodeY, &nodeZ, &nodeType, args->TypeOfProblem);
	if (numNodes == 0) {
		printf("NODES NOT FOUND in file\n");
		return 2;
	}

	/*int *fluid_d = createGpuArrayInt(numNodes, ARRAY_COPY, 0, nodeType);
	FLOAT_TYPE *coordX_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeX);
	FLOAT_TYPE *coordY_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeY);
	FLOAT_TYPE *coordZ_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeZ);*/

	numConns = readConnFile(inFn->bc, &bcNodeIdX, &bcNodeIdY, &bcNodeIdZ,
			&latticeId, &bcType, &bcX, &bcY, &bcZ, &bcBoundId,
			args->TypeOfProblem);
	if (numConns == 0) {
		printf("NEIGHBOURING NOT FOUND in file\n");
		return 2;
	}

	/*int *bcNodeIdX_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX);
	int *bcNodeIdY_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdY);
	int *bcNodeIdZ_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdZ);
	int *latticeId_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, latticeId);
	int *bcType_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcType);
	int *bcBoundId_d = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcBoundId);
	FLOAT_TYPE *bcX_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcX);
	FLOAT_TYPE *bcY_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcY);
	FLOAT_TYPE *bcZ_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcZ);*/

	m = getLastValue(nodeIdY, numNodes);
	n = getLastValue(nodeIdX, numNodes);
	h = getLastValue(nodeIdZ, numNodes);

	delta = getGridSpacing(nodeIdX, nodeIdY, nodeX, numNodes);
	//  printf("checkComment: delta, %f \n",delta);//checkComment
	numInletNodes = getNumInletNodes(bcType, latticeId, numConns,
			args->TypeOfProblem);
	maxInletCoordY = getMaxInletCoordY(bcType, latticeId, bcY, delta, numConns,
			args->TypeOfProblem);
	minInletCoordY = getMinInletCoordY(bcType, latticeId, bcY, delta, numConns,
			args->TypeOfProblem);
	maxInletCoordZ = getMaxInletCoordZ(bcType, latticeId, bcZ, delta, numConns,
			args->TypeOfProblem);
	minInletCoordZ = getMinInletCoordZ(bcType, latticeId, bcZ, delta, numConns,
			args->TypeOfProblem);


	printf("Nx: n= %d \n", n); //checkComment
	printf("Ny: m= %d \n", m); //checkComment
	printf("Nz: h= %d \n", h); //checkComment

	writeInitLog(logFilename, args, delta, m, n, h, numInletNodes,
			maxInletCoordY, minInletCoordY, maxInletCoordZ, minInletCoordZ);
	logFile = fopen(logFilename, "a");
	// In case of no autosave
	sprintf(autosaveFilename, "NOWHERE!");

	initConstants3D(args, maxInletCoordY, minInletCoordY, maxInletCoordZ,
			minInletCoordZ, delta, m, n, h);

	/*dim3 tpb(THREADS, THREADS); 					     // THREADS/block
	dim3 bpg1((int) (sqrt(m * n * h) / THREADS) + 1,
			(int) (sqrt(m * n * h) / THREADS) + 1);       // blocks/grid   MxNxH
	dim3 bpg18((int) (sqrt(18 * m * n * h) / THREADS) + 1,
			(int) (sqrt(18 * m * n * h) / THREADS) + 1);  // blocks/grid 18MxNxH
	dim3 bpg19((int) (sqrt(19 * m * n * h) / THREADS) + 1,
			(int) (sqrt(19 * m * n * h) / THREADS) + 1);  // blocks/grid 19MxNxH
	dim3 bpgBC((int) (sqrt(numConns) / THREADS) + 1,
			(int) (sqrt(numConns) / THREADS) + 1); 	 // blocks/grid N_BC
	*/
	// residuals
	FLOAT_TYPE *norm = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *dragSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *liftSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *latFSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);


	fprintf(logFile, "\n:::: Initializing ::::\n");
	printf("\n:::: Initializing ::::\n");
	//CHECK(cudaEventRecord(start, 0));

	FLOAT_TYPE *u, *v, *w, *rho;

	int InitialCondLoadingErr = -1;
	if (args->UseInitialCondFromFile) {
		InitialCondLoadingErr = readInitConditionsFile(inFn->InitialConditions,
				numNodes, n, m, h, &u, &v, &w, &rho);
	} else {
		u = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		v = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		w = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		rho = createHostArrayFlt(m * n * h, ARRAY_ZERO);
	}
	FLOAT_TYPE *rho_calc;
	if (InitialCondLoadingErr)
		rho_calc = createHostArrayFlt(m * n * h, ARRAY_FILL, args->rho);
	else
		rho_calc = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, rho);
	FLOAT_TYPE *u1, *v1, *w1;
	if (args->inletProfile == NO_INLET) {
		if (InitialCondLoadingErr) {
			u1 = createHostArrayFlt(m * n * h, ARRAY_FILL, args->u);
			v1 = createHostArrayFlt(m * n * h, ARRAY_FILL, args->v);
			w1 = createHostArrayFlt(m * n * h, ARRAY_FILL, args->w);
		} else {
			u1 = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, u);
			v1 = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, v);
			w1 = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, w);
			printf("Initial conditions loaded from file\n");
		}
	}

	if (args->inletProfile == INLET) { 	 //m*h means to do in the inlet face
		printf(
				"Inlet profile is not currently available! Please initiate Inlet profile from file!\n");
		return 0;
		//		cpuInitInletProfile3D(u1, v1,
		//				w1, nodeY, coordZ_d, m * h);
	}
	FLOAT_TYPE *u_prev, *v_prev, *w_prev, *rho_prev, *f_prev;
	if (args->TypeOfResiduals == MacroDiff) {
		if(args->multiPhase)
			f_prev = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		else{
			u_prev = createHostArrayFlt(m * n * h, ARRAY_ZERO);
			v_prev = createHostArrayFlt(m * n * h, ARRAY_ZERO);
			w_prev = createHostArrayFlt(m * n * h, ARRAY_ZERO);
			rho_prev = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		}
	}

	//Multiphase Color Gradient
	FLOAT_TYPE *f, *r_rho, *b_rho, *st_error, *color_gradient, *r_f, *b_f, *r_fColl, *b_fColl;
	if(args->multiPhase){
		r_rho = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		b_rho = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		st_error = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	}

	FLOAT_TYPE aux1 = args->r_density / ((args->r_density + args->b_density) * args->r_viscosity) +
			args->b_density / ((args->r_density + args->b_density) * args->b_viscosity);
	FLOAT_TYPE mean_nu = 1.0/aux1;
	FLOAT_TYPE omega_eff = 1.0/(3.0*mean_nu+0.5);

	FLOAT_TYPE st_predicted = 4.0 * args->A / 9.0 / omega_eff;
#if !CUDA
	f = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
	color_gradient = createHostArrayFlt(m * n * h * 3, ARRAY_ZERO);
	r_f = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
	b_f = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
	r_fColl = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
	b_fColl = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
	FLOAT_TYPE w_pert[19];
	int cx[19] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0 };
	int cy[19] = { 0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1 };
	int cz[19] = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1 };
	FLOAT_TYPE phi[19];
	FLOAT_TYPE teta[19];
	FLOAT_TYPE chi[19];
	FLOAT_TYPE psi[19];
	FLOAT_TYPE cg_w[19];
	FLOAT_TYPE weight[19] =
	{ 1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1.
			/ 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1.
			/ 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1.
			/ 36. };
	if(args->multiPhase){
		int i;
		w_pert[0] = -2.0/ 9.0;
		phi[0]=0;
		teta[0] = 1;
		chi[0] = -5.0/2.0;
		psi[0] = 0;
		cg_w[0] = 0.0;
		for(i = 1; i < 7; i++){
			w_pert[i] = 1.0 / 54.0;
			phi[i] = 1.0 / 12.0;
			teta[i] = -1.0 / 12.0;
			chi[i] = -1.0 / 6.0;
			psi[i] = 1.0 / 4.0;
			cg_w[0] = 1.0 / 6.0;
		}
		for(i = 7; i < 19; i++){
			w_pert[i] = 1.0 / 27.0;
			phi[i] = 1.0 / 24.0;
			teta[i] = -1.0 / 24.0;
			chi[i] = 1.0 / 24.0;
			psi[i] = 1.0 / 8.0;
			cg_w[0] = 1.0 / 12.0;
		}

		createBubble3D(nodeX, nodeY,nodeZ, n, m, h,args->bubble_radius, r_f, b_f,r_rho,b_rho, args->r_density, args->b_density, phi, rho, f);
	}
#endif
	int *cg_directions, *cg_directions;
	FLOAT_TYPE *r_rho, *b_rho, *r_f, *b_f, *r_fColl, *b_fColl, *p_in, *p_out;
	int *num_in, *num_out;
	if(args->multiPhase){
		r_rho = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		b_rho = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		r_f = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		b_f = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		r_fColl = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		b_fColl = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		cg_directions = createHostArrayInt(n * m * h, ARRAY_ZERO);
		if(args->test_case == 1){
			p_in = createHostArrayFlt(n*m*h, ARRAY_ZERO);
			p_out = createHostArrayFlt(n*m*h, ARRAY_ZERO);
			num_in = createHostArrayInt(n*m*h, ARRAY_ZERO);
			num_out = createHostArrayInt(n*m*h, ARRAY_ZERO);
		}
	}

	FLOAT_TYPE *f = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	FLOAT_TYPE *fColl = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	FLOAT_TYPE *f1, *fprev;
	if (args->TypeOfResiduals == FdRelDiff) {
		fprev = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	}

#if CUDA
	FLOAT_TYPE p_in_mean;
	FLOAT_TYPE p_out_mean;
	FLOAT_TYPE ms = n * m * h;
	if(args->multiPhase){
		if(args->high_order)
			initHOColorGradient3D(cg_directions, n, m, h);
		else
			initColorGradient3D(cg_directions, n, m, h);
		//CHECK(cudaMemcpy(cg_directions, cg_directions, SIZEINT(m*n*h), cudaMemcpyHostToDevice));
		initCGBubble3D(nodeX,nodeY,nodeZ,r_rho, b_rho, rho_calc, r_f, b_f, f, args->test_case);
	}
#endif

	if(args->multiPhase){
#if !CUDA
		/*CHECK(cudaMemcpy(f,f,SIZEFLT(m*n*h*19),cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(r_rho,r_rho,SIZEFLT(m*n*h),cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(b_rho,b_rho,SIZEFLT(m*n*h),cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(rho_calc,rho,SIZEFLT(m*n*h),cudaMemcpyHostToDevice));*/
#endif
	}
//Arrived here
	FLOAT_TYPE *temp19a, *temp19b;
	if(args->multiPhase){
		temp19a = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
		temp19b = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	}
	else if (args->TypeOfResiduals != MacroDiff) {
		temp19a = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
		temp19b = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	}
	FLOAT_TYPE *tempA = createHostArrayFlt(m * n * h, ARRAY_ZERO);
	FLOAT_TYPE *tempB = createHostArrayFlt(m * n * h, ARRAY_ZERO);

	int *mask = createHostArrayInt(m * n * h, ARRAY_ZERO);
	unsigned long long *bcMask = createHostArrayLongLong(m * n * h, ARRAY_ZERO);
	int *bcIdx = createHostArrayInt(m * n * h, ARRAY_ZERO);

	/*FLOAT_TYPE *u, *v, *w;

	u = createHostArrayFlt(m * n * h, ARRAY_CPYD, 0, u1);
	v = createHostArrayFlt(m * n * h, ARRAY_CPYD, 0, v1);
	w = createHostArrayFlt(m * n * h, ARRAY_CPYD, 0, w1);*/


	bool *stream = createHostArrayBool(18 * m * n * h, ARRAY_FILL, 1);
	FLOAT_TYPE *q = createHostArrayFlt(18 * m * n * h, ARRAY_FILL, 0.5);

	int bcCount = initBoundaryConditions3D(bcNodeIdX, bcNodeIdY, bcNodeIdZ, q,
			bcBoundId, nodeType, bcX, bcY, bcZ, nodeX, nodeY, nodeZ, latticeId,
			stream, bcType, bcMask, bcIdx, mask, delta, m, n, h, numConns,
			args->boundaryType);

	int *bcIdxCollapsed = createHostArrayInt(bcCount, ARRAY_ZERO);
	unsigned long long *bcMaskCollapsed = createHostArrayLongLong(bcCount,
			ARRAY_ZERO);

	FLOAT_TYPE *qCollapsed;
	if (args->boundaryType == CURVED)
		qCollapsed = createHostArrayFlt(18 * bcCount, ARRAY_ZERO);

	/*dim3 bpgB((int) (sqrt(bcCount) / THREADS) + 1,
			(int) (sqrt(bcCount) / THREADS) + 1); // blocks/grid*/

	collapseBc3D(bcIdx, bcIdxCollapsed, bcMask, bcMaskCollapsed, q,
			qCollapsed, mask, m, n, h, bcCount, args->boundaryType);


/*#if CUDA
	CHECK(cudaMemcpy(u, u_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(v, v_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(w, w_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(rho, rho_calc, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	if(args->multiPhase){
		CHECK(cudaMemcpy(r_rho, r_rho, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(b_rho, b_rho, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	}
#endif
	CHECK(cudaEventRecord(stop, 0));
	CHECK(cudaEventSynchronize(stop));
	CHECK(cudaEventElapsedTime(&cudatime, start, stop));
	taskTime[T_INIT] += cudatime / 1000;*/

	fclose(logFile);
	writeNodeNumbers(logFilename, numNodes, numConns, bcCount);
	logFile = fopen(logFilename, "a");

	void *hostArrays[] = { nodeIdX, nodeIdY, nodeIdZ, nodeX, nodeY, nodeZ,
			nodeType, bcNodeIdX, bcNodeIdY, bcNodeIdZ, latticeId, bcType, bcX,
			bcY, bcZ, bcBoundId, u, v, w, rho, mask, bcMask, bcIdx, stream, q,
			norm, dragSum, liftSum, latFSum, };

	void *gpuArrays[] =
	{ rho_calc, u1, v1, w1, f, fColl, tempA, tempB,
			bcMaskCollapsed, bcIdx, bcIdxCollapsed,
			stream, qCollapsed}; //drag_d, lift_d, latF_d,


	void *mpHostArrays[] = {
			r_rho, b_rho, st_error, color_gradient, r_f,b_f, r_fColl, b_fColl, cg_directions, f
	};


	void *mpTC1GpuArrays[] = {
			p_in,	p_out, num_in, num_out
	};

	void *FDdifGpuArrays[] = {
			fprev, f1
	};

	void *nonMacroDiffGpuArrays[] = {
			temp19a, temp19b
	};

	fprintf(logFile, "\n:::: Initialization done! ::::\n");

	printf("Initialization took %f seconds\n", taskTime[T_INIT]);

	// Write Initialized data
	switch (args->outputFormat) {
	case CSV:
		sprintf(outputFilename, "%sInitialData.csv", inFn->result);
		break;
	case TECPLOT:
		sprintf(outputFilename, "%sInitialData.dat", inFn->result);
		break;
	case PARAVIEW:
		sprintf(outputFilename, "%sInitialData.vti", inFn->result);
		break;
	}

	tInstant1 = clock(); // Start measuring time
	if(args->multiPhase){
		WriteResultsMultiPhase(outputFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,r_rho,b_rho, nodeType,
				n, m, h, args->outputFormat);
	}
	else{
		WriteResults3D(outputFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,
				nodeType, n, m, h, args->outputFormat);
	}
	tInstant2 = clock();
	taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1) / CLOCKS_PER_SEC;

	printf("\nInitialized data was written to %s\n", outputFilename);

	////////////////// ITERATION ///////////////////////

	fprintf(logFile, "\n:::: Start Iterations ::::\n");
	printf("\n:::: Start Iterations ::::\n");

	printf("%d is the number of iterations \n", args->iterations);

	tIterStart = clock(); // Start measuring time of main loop
	size_t free, total;

	/*cuMemGetInfo(&free, &total);
	printf("^^^^ Free : %llu Mbytes \n",
			(unsigned long long) free / 1024 / 1024);

	printf("^^^^ Total: %llu Mbytes \n",
			(unsigned long long) total / 1024 / 1024);

	printf("^^^^ %f%% free, %f%% used\n", 100.0 * free / (double) total,
			100.0 * (total - free) / (double) total);*/
	int iter = 0;
	while (iter < args->iterations) {
		//CHECK(cudaThreadSynchronize());
		//CHECK(cudaEventRecord(start, 0)); // Start measuring time
		////////////// COLLISION ///////////////
		switch (args->collisionModel) {
		case BGKW:
			if(args->multiPhase){
#if !CUDA
				mp3DColl(n, m, h, rho, u, v, w, r_rho, b_rho, w_pert, color_gradient, args->beta,
						args->g_limit, args->A,r_fColl, b_fColl, weight, cx, cy, cz, f, args->r_viscosity,
						args->b_viscosity, args->r_alpha, args->b_alpha, chi, phi, psi, teta, cg_w);
#else
				if(!args->enhanced_distrib)
					cpuCollBgkwGC3D(nodeType, rho_calc, r_rho, b_rho, u, v, w, f, r_fColl, b_fColl, cg_directions, args->high_order);
				else
					cpuCollEnhancedBgkwGC3D(nodeType, rho_calc, r_rho, b_rho, u, v, w, f, r_fColl, b_fColl, cg_directions, args->high_order);
#endif
			}
			else{
				cpuCollBgkw3D(nodeType, rho_calc, u, v, w, f,
						fColl);
			}
			break;

		case TRT:
			printf("TRT not implemented in 3D go for MRT \n");
			//        cpuCollTrt(nodeType, rho_calc, u, v, w, f, fColl);
			break;

		case MRT:
			cpuCollMrt3D(nodeType, rho_calc, u, v, w, f,
					fColl);
			break;
		}

		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_COLL] += cudatime;

		////////////// STREAMING ///////////////
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0));*/

		if(args->multiPhase){
#if !CUDA
			streamMP3D(n, m, h, r_f, b_f, r_fColl, b_fColl, stream);
#else
			cpuStreaming3D(nodeType, stream, r_f, r_fColl);
			cpuStreaming3D(nodeType, stream, b_f, b_fColl);
#endif
		}
		else{
			cpuStreaming3D(nodeType, stream, f, fColl);
		}

		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_STRM] += cudatime;

		// make the host block until the device is finished with foo
		CHECK(cudaThreadSynchronize());

		// check for error
		cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess) {
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			exit(-1);
		}*/

		////////////// BOUNDARIES ///////////////
		//CHECK(cudaThreadSynchronize());
		//CHECK(cudaEventRecord(start, 0));

		if(args->multiPhase){
#if !CUDA
			peridicBoundaries3D(n, m, h,r_f, b_f, r_rho, b_rho);
#else
			cpuBcInlet3D(bcIdxCollapsed, bcMaskCollapsed, r_f,
					u1, v1, w1, bcCount);
			cpuBcInlet3D(bcIdxCollapsed, bcMaskCollapsed, b_f,
					u1, v1, w1, bcCount);
			switch (args->bcwallmodel) {
			case SIMPLE:
				cpuBcSimpleWall3D(bcIdxCollapsed,
						bcMaskCollapsed, r_f, r_fColl, qCollapsed, bcCount, args->boundaryType);
				cpuBcSimpleWall3D(bcIdxCollapsed,
						bcMaskCollapsed, b_f, b_fColl, qCollapsed, bcCount, args->boundaryType);

				break;
			case COMPLEX:
				cpuBcComplexWall3D(bcIdxCollapsed,
						bcMaskCollapsed, r_f, r_fColl, qCollapsed, bcCount);
				cpuBcComplexWall3D(bcIdxCollapsed,
						bcMaskCollapsed, b_f, b_fColl, qCollapsed, bcCount);

				break;
			}

			cpuBcPeriodic3D(bcIdxCollapsed, bcMaskCollapsed, r_f,
					bcCount);
			cpuBcPeriodic3D(bcIdxCollapsed, bcMaskCollapsed, b_f,
					bcCount);
#endif
		}
		else{
			cpuBcInlet3D(bcIdxCollapsed, bcMaskCollapsed, f,
					u1, v1, w1, bcCount);
			switch (args->bcwallmodel) {
			case SIMPLE:
				cpuBcSimpleWall3D(bcIdxCollapsed,
						bcMaskCollapsed, f, fColl, qCollapsed, bcCount, args->boundaryType);

				break;
			case COMPLEX:
				cpuBcComplexWall3D(bcIdxCollapsed,
						bcMaskCollapsed, f, fColl, qCollapsed, bcCount);

				break;
			}
			cpuBcOutlet3D(bcIdxCollapsed, bcMaskCollapsed, f,
					u, v, w, rho_calc, bcCount, args->outletProfile);
			cpuBcPeriodic3D(bcIdxCollapsed, bcMaskCollapsed, f,
					bcCount);
			cpuBcSymm3D(bcIdxCollapsed, bcMaskCollapsed, f,
					bcCount);
		}
		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_BNDC] += cudatime;

		// UPDATE VELOCITY AND DENSITY
		CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0));*/

		if(args->multiPhase){
#if !CUDA
			updateMacroMP3D(n, m, h, u, v, w, r_rho, b_rho, r_f, b_f, rho, args->control_param,args->r_alpha, args->b_alpha,
					args->bubble_radius,st_error, iter, 1, cx, cy, cz, f);
#else
			cpuUpdateMacro3DCG(nodeType, rho_calc, u, v, w,
					bcBoundId, f, args->g,bcMask,args->UpdateInltOutl, r_f, b_f, r_rho, b_rho, p_in, p_out, num_in, num_out, args->test_case);
			switch(args->test_case){
			case 1:
				p_in_mean = cpu_sum_h(p_in, p_in, ms) / cpu_sum_int_h(num_in, num_in, ms);
				p_out_mean = cpu_sum_h(p_out, p_out, ms) / cpu_sum_int_h(num_out, num_out, ms);
				st_error[iter] = calculateSurfaceTension3D(p_in_mean, p_out_mean,args->r_alpha, args->b_alpha, args->bubble_radius * n, st_predicted);
				break;
			default:
				break;
			}
#endif
		}
		else{
			cpuUpdateMacro3D(nodeType, rho_calc, u, v, w,
					bcBoundId, nodeX, nodeY, nodeZ, f, args->g,bcMask,args->UpdateInltOutl);
		}
		tInstant2 = clock();
		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_MACR] += cudatime;*/

		// COMPUTE RESIDUALS

		if (AuxMacroDiff * args->ShowMacroDiff == iter + 1) {
			//CHECK(cudaThreadSynchronize());
			//CHECK(cudaEventRecord(start, 0));

			if (args->TypeOfResiduals == L2) {
				if(args->multiPhase){
#if !CUDA
					//CHECK(cudaMemcpy(f,r_f,SIZEFLT(m*n*h*19),cudaMemcpyHostToDevice));
					memcpy(&f,&r_f, SIZEFLT(m*n*h*19));
					//CHECK(cudaMemcpy(fColl,r_fColl,SIZEFLT(m*n*h*19),cudaMemcpyHostToDevice));
					memcpy(&fColl,&r_fColl, SIZEFLT(m*n*h*19));
#endif
				}
				r = computeResidual3D(f, fColl, temp19a, temp19b, m, n,
						h);
			}
			else {
				if (args->TypeOfResiduals == FdRelDiff) {
					if(args->multiPhase){
#if !CUDA
						//CHECK(cudaMemcpy(f,r_f,SIZEFLT(m*n*h*19),cudaMemcpyHostToDevice));
						memcpy(&f, &r_f, SIZEFLT(m*n*h*19));
#endif
					}

					if (firstIter) {
						firstIter = false;
						f1 = createHostArrayFlt(19 * n * m * h, ARRAY_CPYD, 0,
								f);

					}
					r = computeNewResidual3D(f, fprev, f1, temp19a,
							temp19b, m, n, h);
					//CHECK(cudaFree(fprev));
					free(fprev);
					fprev = createHostArrayFlt(19 * n * m * h, ARRAY_CPYD, 0,
							f);

				} else {
					bool h_divergence = false;
					//CHECK(cudaMalloc(&d_divergence,sizeof(bool)));
					//CHECK(cudaMemcpy(d_divergence,&h_divergence,sizeof(bool),cudaMemcpyHostToDevice));
					if(args->multiPhase){
						cpu_abs_sub(f, f_prev, temp19a, n * m * h * 19, h_divergence);
						fMaxDiff = cpu_max_h(temp19a, temp19b, n * m * h * 19);
					}
					else{
						cpu_abs_sub(u, u_prev, tempA,
								n * m * h, h_divergence);
						uMaxDiff = cpu_max_h(tempA, tempB, n * m * h);
						cpu_abs_sub(v, v_prev, tempA,
								n * m * h, h_divergence);
						vMaxDiff = cpu_max_h(tempA, tempB, n * m * h);
						cpu_abs_sub(w, w_prev, tempA,
								n * m * h, h_divergence);
						wMaxDiff = cpu_max_h(tempA, tempB, n * m * h);
						cpu_abs_sub(rho_calc, rho_prev, tempA,
								n * m * h, h_divergence);
						rhoMaxDiff = cpu_max_h(tempA, tempB, n * m * h);
					}
					//CHECK(cudaMemcpy(&h_divergence,h_divergence,sizeof(bool),cudaMemcpyDeviceToHost));
					//CHECK(cudaFree(d_divergence));
					if (h_divergence) {
						fprintf(stderr, "\nDIVERGENCE!\n");
						break;
					}

					if(args->multiPhase){
						if(abs(fMaxDiff) < args->StopCondition[0]){
							printf("simulation converged!\n");
							break;
						}
					}
					else if (abs(uMaxDiff) < args->StopCondition[0] &&
							abs(vMaxDiff) < args->StopCondition[1] &&
							abs(wMaxDiff) < args->StopCondition[2] &&
							abs(rhoMaxDiff) < args->StopCondition[3]) {
						printf("simulation converged!\n");
						break;
					}

					if(args->multiPhase){
						//CHECK(cudaFree(f_prev));
						free(f_prev);
						f_prev = createHostArrayFlt(n * m * h * 19, ARRAY_CPYD, 0, f);
					}else{
						writeMacroDiffs(iter + 1, uMaxDiff, vMaxDiff, wMaxDiff,	rhoMaxDiff);
						//CHECK(cudaFree(u_prev));
						//CHECK(cudaFree(v_prev));
						//CHECK(cudaFree(w_prev));
						//CHECK(cudaFree(rho_prev));
						free(u_prev);
                        free(v_prev);
                        free(w_prev);
                        free(rho_prev);
						u_prev = createHostArrayFlt(n * m * h, ARRAY_CPYD, 0, u);
						v_prev = createHostArrayFlt(n * m * h, ARRAY_CPYD, 0, v);
						w_prev = createHostArrayFlt(n * m * h, ARRAY_CPYD, 0, w);
						rho_prev = createHostArrayFlt(n * m * h, ARRAY_CPYD, 0,
								rho_calc);
					}
				}
			}

			if (abs(r) < args->StopCondition[0]) {
				printf("simulation converged!\n");
				break;
			}
			if (r != r) {
				fprintf(stderr, "\nDIVERGENCE!\n");
				break;
			}

			/*CHECK(cudaEventRecord(stop, 0));
			CHECK(cudaEventSynchronize(stop));
			CHECK(cudaEventElapsedTime(&cudatime, start, stop));
			taskTime[T_RESI] += cudatime;*/

			AuxMacroDiff++;

		}
		if(args->multiPhase){
			//CHECK(cudaFree(f_prev));
			free(f_prev);
			f_prev = createHostArrayFlt(n * m * h * 19, ARRAY_CPYD, 0, f);
		}
		norm[iter] = r;
		if(args->multiPhase){
			printf(
					"Iterating... %d/%d (%3.1f %%) Max macro diffs: f= %.10f\r",
					iter + 1, args->iterations,
					(FLOAT_TYPE) (iter + 1) * 100
					/ (FLOAT_TYPE) (args->iterations), fMaxDiff);
		}
		else if (args->TypeOfResiduals == MacroDiff) {
			printf(
					"Iterating... %d/%d (%3.1f %%) Max macro diffs: u= %.10f v= %.10f w= %.10f rho= %.10f \r",
					iter + 1, args->iterations,
					(FLOAT_TYPE) (iter + 1) * 100
					/ (FLOAT_TYPE) (args->iterations), uMaxDiff,
					vMaxDiff, wMaxDiff, rhoMaxDiff);
		} else {
			printf("Iterating... %d/%d (%3.1f %%)  residual="FLOAT_FORMAT" \r",
					iter + 1, args->iterations,
					(FLOAT_TYPE) (iter + 1) * 100
					/ (FLOAT_TYPE) (args->iterations), r);
		}

		iter++; // update loop variable

		////////////// Autosave ///////////////

		if (iter == (args->autosaveEvery * autosaveIt)) {
			autosaveIt++;
			if (iter > args->autosaveAfter) {
				printf("autosave\n\n");
				//////////// COPY VARIABLES TO HOST ////////////////
				/*CHECK(cudaMemcpy(u, u_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
				CHECK(cudaMemcpy(v, v_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
				CHECK(cudaMemcpy(w, w_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
				CHECK(cudaMemcpy(rho, rho_calc, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));*/

				switch (args->outputFormat) {
				case CSV:
					sprintf(autosaveFilename, "%sautosave_iter%05d.csv",
							inFn->result, iter);
					break;
				case TECPLOT:
					sprintf(autosaveFilename, "%sautosave_iter%05d.dat",
							inFn->result, iter);
					break;
				case PARAVIEW:
					sprintf(autosaveFilename, "%sautosave_iter%05d.vti",
							inFn->result, iter);
					break;
				}

				tInstant1 = clock(); // Start measuring time
				WriteResults3D(autosaveFilename, nodeType, nodeX, nodeY, nodeZ,
						u, v, w, rho, nodeType, n, m, h, args->outputFormat);
				tInstant2 = clock();
				taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1)
																																																																																														/ CLOCKS_PER_SEC;
			}
		}
	}     ////////////// END OF MAIN WHILE CYCLE! ///////////////

	tIterEnd = clock(); // End measuring time of main loop
	taskTime[T_ITER] = (FLOAT_TYPE) (tIterEnd - tIterStart) / CLOCKS_PER_SEC;

	clock_t tEnd = clock();
	taskTime[T_OALL] = (FLOAT_TYPE) (tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time
	taskTime[T_COLL] /= 1000;
	taskTime[T_STRM] /= 1000;
	taskTime[T_BNDC] /= 1000;
	taskTime[T_MACR] /= 1000;
	taskTime[T_RESI] /= 1000;

	fclose(logFile);
	writeEndLog(logFilename, taskTime);
	writeTimerLog(timeFilename, taskTime);
	if (args->TypeOfResiduals != MacroDiff) {
		writeResiduals(residualsFilename, norm, dragSum, liftSum, m * n * h,
				args->iterations);
	}
	// Write final data
/*#if CUDA
	CHECK(cudaMemcpy(u, u_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(v, v_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(w, w_d, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(rho, rho_calc, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	if(args->multiPhase){
		CHECK(cudaMemcpy(r_rho, r_rho, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
		CHECK(cudaMemcpy(b_rho, b_rho, SIZEFLT(m*n*h), cudaMemcpyDeviceToHost));
	}
#endif*/

	switch (args->outputFormat) {
	case CSV:
		sprintf(finalFilename, "%sFinalData.csv", inFn->result);
		break;
	case TECPLOT:
		sprintf(finalFilename, "%sFinalData.dat", inFn->result);
		break;
	case PARAVIEW:
		sprintf(finalFilename, "%sFinalData.vti", inFn->result);
		break;
	}
	if(args->multiPhase){
		FLOAT_TYPE *analytical = createHostArrayFlt(m, ARRAY_ZERO);
		switch (args->test_case) {
		case 1:
			printf("Suface tension error: "FLOAT_FORMAT"\n", st_error[iter-1]);
			WriteArray("surface tension",st_error, args->iterations,1);
			break;
		case 2:
			deformingBubbleValid(r_rho, b_rho, n, m, h);
			break;
		case 3:
			validateCoalescenceCase(r_rho, b_rho, n, m, args->bubble_radius, h);
			break;
		case 4: //COUETTE
			analyticalCouette(args->kappa, nodeY, m, n, analytical, args->u, h);
			writeCouetteSolution("Profile_Couette", analytical, u, nodeY, m, n, h);
			printf("Couette profile written to Profile_Couette in Results/\n");
			break;
		default:
			break;
		}
		WriteResultsMultiPhase(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,r_rho,b_rho, nodeType,
				n, m, h, args->outputFormat);
	}
	else{
		WriteResults3D(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,
				nodeType, n, m, h, args->outputFormat);
	}

	WriteLidDrivenCavityMidLines3D(nodeX, nodeY, nodeZ, u, w, n, m, h, args->u);
	WriteChannelCrossSection3D(nodeX, nodeY, nodeZ, u, v, w, n, m, h, args->u);

	// Write information for user
	printf("\n\nLog was written to %s\n", logFilename);
	printf("Last autosave result can be found at %s\n", autosaveFilename);
	printf("residuals were written to %s\n", residualsFilename);
	printf("Profiling results were written to %s\n", timeFilename);
	//	compareTestFiles("./TestValues/CUDA/cg2d.txt", "./TestValues/CUDA/cg3d.txt");
	printf("Final results were written to %s\n", finalFilename);

	//cudaEventDestroy(start);
	//cudaEventDestroy(stop);
	freeAllHost(hostArrays, sizeof(hostArrays) / sizeof(hostArrays[0]));
	freeAllHost(gpuArrays, sizeof(gpuArrays) / sizeof(gpuArrays[0]));
	if(args->multiPhase){
#if CUDA
		if(args->test_case == 1)
			freeAllHost(mpTC1HostArrays, sizeof(mpTC1HostArrays) / sizeof(mpTC1HostArrays[0]));
#else
		freeAllHost(mpHostArrays, sizeof(mpHostArrays) / sizeof(mpHostArrays[0]));
#endif
	}
	if (args->TypeOfResiduals == FdRelDiff) {
		freeAllHost(FDdifHostArrays, sizeof(FDdifHostArrays) / sizeof(FDdifHostArrays[0]));
	}
	if (args->TypeOfResiduals != MacroDiff) {
		freeAllHost(nonMacroDiffHostArrays, sizeof(nonMacroDiffHostArrays) / sizeof(nonMacroDiffHostArrays[0]));
	}

	return 0;
}
