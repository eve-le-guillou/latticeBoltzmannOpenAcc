/**
 * @author Alfonso Aguilar (a.aguilar-pontes@cranfield.ac.uk) - implementation of the physics
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) - software aspects of the implementation
 */

#include "Arguments.h"
#include "GpuConstants.h"
#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <errno.h>
#include "GpuFunctions.h"       // GPU kernels
#include "ShellFunctions.h"     // For convenience
#include "FilesReading.h"       // For reading files
#include "FilesWriting.h"       // For writing files e.g. tecplot
#include "CellFunctions.h"      // For cell modifications
#include "ComputeResiduals.h"   // residuals
#include "LogWriter.h"
#include "Iterate.h"
#include "ArrayUtils.h"
#include "GpuSum.h"
#include "Multiphase.h"

#ifdef MAKE_SERIAL
#define OPENACC 0
#else
#define OPENACC 1
#include <cuda_runtime.h>
#endif


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
	bool *d_divergence;
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
	*bcNodeIdZ, *latticeId, *bcType, *bcBoundId_d;
	FLOAT_TYPE *nodeX, *nodeY, *nodeZ, *bcX, *bcY, *bcZ;

	FLOAT_TYPE taskTime[9];
	int i;
	for (i = 0; i < 9; ++i) {
		taskTime[i] = 0.0;
	}

	clock_t tInstant1, tInstant2; // Time measurement points, universal
	clock_t tIterStart, tIterEnd; // Time measurement points: main loop

	numNodes = readNodeFile(inFn->node, &nodeIdX, &nodeIdY, &nodeIdZ, &nodeX,
			&nodeY, &nodeZ, &nodeType, args->TypeOfProblem);
	if (numNodes == 0) {
		printf("NODES NOT FOUND in file\n");
		return 2;
	}

	numConns = readConnFile(inFn->bc, &bcNodeIdX, &bcNodeIdY, &bcNodeIdZ,
			&latticeId, &bcType, &bcX, &bcY, &bcZ, &bcBoundId_d,
			args->TypeOfProblem);
	if (numConns == 0) {
		printf("NEIGHBOURING NOT FOUND in file\n");
		return 2;
	}

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

	// residuals
	FLOAT_TYPE *norm = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *dragSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *liftSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *latFSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);


	fprintf(logFile, "\n:::: Initializing ::::\n");
	printf("\n:::: Initializing ::::\n");

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
	FLOAT_TYPE *rho_d;
	if (InitialCondLoadingErr)
		rho_d = createHostArrayFlt(m * n * h, ARRAY_FILL, args->rho);
	else
		rho_d = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, rho);
	FLOAT_TYPE *u1_d, *v1_d, *w1_d;
	if (args->inletProfile == NO_INLET) {
		if (InitialCondLoadingErr) {
			u1_d = createHostArrayFlt(m * n * h, ARRAY_FILL, args->u);
			v1_d = createHostArrayFlt(m * n * h, ARRAY_FILL, args->v);
			w1_d = createHostArrayFlt(m * n * h, ARRAY_FILL, args->w);
		} else {
			u1_d = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, u);
			v1_d = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, v);
			w1_d = createHostArrayFlt(m * n * h, ARRAY_COPY, 0, w);
			printf("Initial conditions loaded from file\n");
		}
	}

	if (args->inletProfile == INLET) { 	 //m*h means to do in the inlet face
		printf(
				"Inlet profile is not currently available! Please initiate Inlet profile from file!\n");
		return 0;
		//		gpuInitInletProfile3D<<<(int) (m * h / THREADS) + 1, tpb>>>(u1_d, v1_d,
		//				w1_d, nodeY, nodeZ, m * h);
	}
	FLOAT_TYPE *u_prev_d, *v_prev_d, *w_prev_d, *rho_prev_d, *f_prev_d;
	if (args->TypeOfResiduals == MacroDiff) {
		if(args->multiPhase)
			f_prev_d = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		else{
			u_prev_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
			v_prev_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
			w_prev_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
			rho_prev_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		}
	}

	//Multiphase Color Gradient
	FLOAT_TYPE *st_error;
	if(args->multiPhase){
		st_error = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	}

	FLOAT_TYPE aux1 = args->r_density / ((args->r_density + args->b_density) * args->r_viscosity) +
			args->b_density / ((args->r_density + args->b_density) * args->b_viscosity);
	FLOAT_TYPE mean_nu = 1.0/aux1;
	FLOAT_TYPE omega_eff = 1.0/(3.0*mean_nu+0.5);

	FLOAT_TYPE st_predicted = 4.0 * args->A / 9.0 / omega_eff;

	int *cg_dir_d;
	FLOAT_TYPE *r_rho_d, *b_rho_d, *r_f_d, *b_f_d, *r_fColl_d, *b_fColl_d, *p_in_d, *p_out_d;
	int *num_in_d, *num_out_d;
	if(args->multiPhase){
		r_rho_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		b_rho_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
		r_f_d = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		b_f_d = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		r_fColl_d = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		b_fColl_d = createHostArrayFlt(m * n * h * 19, ARRAY_ZERO);
		cg_dir_d = createHostArrayInt(m * n * h, ARRAY_ZERO);
		if(args->test_case == 1){
			p_in_d = createHostArrayFlt(n*m*h, ARRAY_ZERO);
			p_out_d = createHostArrayFlt(n*m*h, ARRAY_ZERO);
			num_in_d = createHostArrayInt(n*m*h, ARRAY_ZERO);
			num_out_d = createHostArrayInt(n*m*h, ARRAY_ZERO);
		}
	}

	FLOAT_TYPE *f_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	FLOAT_TYPE *fColl_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	FLOAT_TYPE *f1_d, *fprev_d;
	if (args->TypeOfResiduals == FdRelDiff) {
		fprev_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	}
	f1_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);

	FLOAT_TYPE p_in_mean;
	FLOAT_TYPE p_out_mean;
	int ms = n * m * h;
	if(args->multiPhase){
		if(args->high_order)
			initHOColorGradient3D(cg_dir_d, n, m, h);
		else
			initColorGradient3D(cg_dir_d, n, m, h);
#pragma acc enter data copyin(nodeX[0:ms], nodeY[0:ms], nodeZ[0:ms], rho_d[0:ms],r_rho_d[0:ms], b_rho_d[0:ms], r_f_d[0:ms*19], b_f_d[0:ms*19], f_d[0:ms*19])
		initCGBubble3D(nodeX,nodeY,nodeZ,r_rho_d, b_rho_d, rho_d, r_f_d, b_f_d, f_d, args->test_case);
	}

	if(args->multiPhase){
	}

	FLOAT_TYPE *temp19a_d, *temp19b_d;
	if(args->multiPhase){
		temp19a_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
		temp19b_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	}
	else if (args->TypeOfResiduals != MacroDiff) {
		temp19a_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
		temp19b_d = createHostArrayFlt(19 * m * n * h, ARRAY_ZERO);
	}
	FLOAT_TYPE *tempA_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);
	FLOAT_TYPE *tempB_d = createHostArrayFlt(m * n * h, ARRAY_ZERO);

	int *mask = createHostArrayInt(m * n * h, ARRAY_ZERO);
	unsigned long long *bcMask = createHostArrayLongLong(m * n * h, ARRAY_ZERO);
	int *bcIdx = createHostArrayInt(m * n * h, ARRAY_ZERO);

	FLOAT_TYPE *u_d, *v_d, *w_d;

	u_d = createHostArrayFlt(m * n * h, ARRAY_CPYD, 0, u1_d);
	v_d = createHostArrayFlt(m*n*h, ARRAY_ZERO);
	w_d = createHostArrayFlt(m*n*h, ARRAY_ZERO);
//	v_d = createHostArrayFlt(m * n * h,ARRAY_CPYD, 0, v1_d);
//	w_d = createHostArrayFlt(m * n * h, ARRAY_CPYD, 0, w1_d);


        bool *stream_d = createHostArrayBool(18 * m * n * h, ARRAY_FILL,1);
	FLOAT_TYPE *q = createHostArrayFlt(18 * m * n * h, ARRAY_FILL, 0.5);

	int bcCount = initBoundaryConditions3D(bcNodeIdX, bcNodeIdY, bcNodeIdZ, q,
			bcBoundId_d, nodeType, bcX, bcY, bcZ, nodeX, nodeY, nodeZ, latticeId,
			stream_d, bcType, bcMask, bcIdx, mask, delta, m, n, h, numConns,
			args->boundaryType);
	unsigned long long *bcMask_d = createHostArrayLongLong(m * n * h, ARRAY_COPY,
			0, bcMask);
	int *bcIdxCollapsed_d = createHostArrayInt(bcCount, ARRAY_ZERO);
	unsigned long long *bcMaskCollapsed_d = createHostArrayLongLong(bcCount,
			ARRAY_ZERO);

	FLOAT_TYPE *qCollapsed_d;
	if (args->boundaryType == CURVED)
		qCollapsed_d = createHostArrayFlt(18 * bcCount, ARRAY_ZERO);

	collapseBc3D(bcIdx, bcIdxCollapsed_d, bcMask, bcMaskCollapsed_d, q,
			qCollapsed_d, mask, m, n, h, bcCount, args->boundaryType);

	if(args->multiPhase){
	}
	fclose(logFile);
	writeNodeNumbers(logFilename, numNodes, numConns, bcCount);
	logFile = fopen(logFilename, "a");

	void *hostArrays[] = { nodeIdX, nodeIdY, nodeIdZ, nodeX, nodeY, nodeZ,
			nodeType, bcNodeIdX, bcNodeIdY, bcNodeIdZ, latticeId, bcType, bcX,
			bcY, bcZ, u, v, w, rho, mask, bcMask, bcIdx, q,
			norm, dragSum, liftSum, latFSum};

	void *gpuArrays[] ={ bcBoundId_d, u_d, v_d, w_d, rho_d, u1_d, v1_d, w1_d, f_d, fColl_d, tempA_d, tempB_d, bcMaskCollapsed_d, bcIdxCollapsed_d, stream_d/*, qCollapsed_d*/}; //drag_d, lift_d, latF_d,


	void *mpHostArrays[] = {st_error};

	void *mpGpuArrays[] = {r_rho_d, b_rho_d, r_f_d, b_f_d, r_fColl_d, b_fColl_d, cg_dir_d };

	void *mpTC1GpuArrays[] = {
			p_in_d,	p_out_d, num_in_d, num_out_d
	};

	void *FDdifGpuArrays[] = {
			fprev_d, f1_d
	};

	void *nonMacroDiffGpuArrays[] = {
			temp19a_d, temp19b_d
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
		WriteResultsMultiPhase(outputFilename, nodeType, nodeX, nodeY, nodeZ, u_d, v_d, w_d, rho_d, r_rho_d, b_rho_d, nodeType, n, m, h, args->outputFormat);
	}
	else{
		WriteResults3D(outputFilename, nodeType, nodeX, nodeY, nodeZ, u_d, v_d, w_d, rho_d, nodeType, n, m, h, args->outputFormat);
	}
	tInstant2 = clock();
	taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1) / CLOCKS_PER_SEC;

	printf("\nInitialized data was written to %s\n", outputFilename);

	////////////////// ITERATION ///////////////////////

	fprintf(logFile, "\n:::: Start Iterations ::::\n");
	printf("\n:::: Start Iterations ::::\n");

	printf("%d is the number of iterations \n", args->iterations);

	tIterStart = clock(); // Start measuring time of main loop

	int iter = 0;
#pragma acc data copyin(u_d[0:ms], v_d[0:ms], w_d[0:ms],  nodeType[0:numNodes], stream_d[0:18*ms], bcIdxCollapsed_d[0:bcCount], bcMaskCollapsed_d[0:bcCount], qCollapsed_d[0:bcCount*18]) create(r_fColl_d[ms*19], b_fColl_d[ms*19], p_in_d[0:ms], p_out_d[0:ms], num_in_d[0:ms], num_out_d[0:ms], d_divergence[0:1]) \
                copyin(cg_dir_d[0:ms], bcMask_d[0:ms], f_prev_d[0:19*ms], fprev_d[0:19],  f1_d[0:19*ms], temp19a_d[0:ms*19], temp19b_d[0:ms*19], u1_d[0:ms], v1_d[0:ms], w1_d[0:ms]) 
		
 //               copyin(nodeX[0:ms],nodeY[0:ms], nodeZ[0:ms],  rho_d[0:ms],r_rho_d[0:ms], b_rho_d[0:ms], r_f_d[0:ms*19], b_f_d[0:ms*19], f_d[0:ms*19])
{
	while (iter < args->iterations) {
		////////////// COLLISION ///////////////
		switch (args->collisionModel) {
		case BGKW:
			if(args->multiPhase){
				if(!args->enhanced_distrib)
					gpuCollBgkwGC3D(nodeType, rho_d, r_rho_d, b_rho_d, u_d, v_d, w_d, f_d, r_fColl_d, b_fColl_d, cg_dir_d, args->high_order);
				else
					gpuCollEnhancedBgkwGC3D(nodeType, rho_d, r_rho_d, b_rho_d, u_d, v_d, w_d, f_d, r_fColl_d, b_fColl_d, cg_dir_d, args->high_order);
			}
			else{
			//	gpuCollBgkw3D<<<bpg1, tpb>>>(nodeType, rho_d, u_d, v_d, w_d, f_d, fColl_d);
			}
			break;

		case TRT:
			printf("TRT not implemented in 3D go for MRT \n");
			//        gpuCollTrt<<<bpg1,tpb>>>(nodeType, rho_d, u_d, v_d, w_d, f_d, fColl_d);
			break;

		case MRT:
			//gpuCollMrt3D<<<bpg1, tpb>>>(nodeType, rho_d, u_d, v_d, w_d, f_d, fColl_d);
			break;
		}

		////////////// STREAMING ///////////////
		if(args->multiPhase){
			gpuStreaming3D(nodeType, stream_d, r_f_d, r_fColl_d);
			gpuStreaming3D(nodeType, stream_d, b_f_d, b_fColl_d);
		}
		else{
		//	gpuStreaming3D<<<bpg1, tpb>>>(nodeType, stream_d, f_d, fColl_d);
		}

		////////////// BOUNDARIES ///////////////

		if(args->multiPhase){
			gpuBcInlet3D(bcIdxCollapsed_d, bcMaskCollapsed_d, r_f_d,
					u1_d, v1_d, w1_d, bcCount);
			gpuBcInlet3D(bcIdxCollapsed_d, bcMaskCollapsed_d, b_f_d,
					u1_d, v1_d, w1_d, bcCount);
			switch (args->bcwallmodel) {
			case SIMPLE:
				
				gpuBcSimpleWall3D(bcIdxCollapsed_d,
						bcMaskCollapsed_d, r_f_d, r_fColl_d, qCollapsed_d, bcCount);
				gpuBcSimpleWall3D(bcIdxCollapsed_d,
						bcMaskCollapsed_d, b_f_d, b_fColl_d, qCollapsed_d, bcCount);
				
				break;
			case COMPLEX:
				/*gpuBcComplexWall3D(bcIdxCollapsed_d,
						bcMaskCollapsed_d, r_f_d, r_fColl_d, qCollapsed_d, bcCount);
				gpuBcComplexWall3D(bcIdxCollapsed_d,
						bcMaskCollapsed_d, b_f_d, b_fColl_d, qCollapsed_d, bcCount);
				*/
				break;
			}

			gpuBcPeriodic3D(bcIdxCollapsed_d, bcMaskCollapsed_d, r_f_d,
					bcCount);
			gpuBcPeriodic3D(bcIdxCollapsed_d, bcMaskCollapsed_d, b_f_d,
					bcCount);
		}
		else{
		/*	gpuBcInlet3D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					u1_d, v1_d, w1_d, bcCount);
			switch (args->bcwallmodel) {
			case bcIdxCollapsed_d,
						bcMaskCollapsed_d, f_d, fColl_d, qCollapsed_d, bcCount);

				break;
			case COMPLEX:
				gpuBcComplexWall3D<<<bpgB, tpb>>>(bcIdxCollapsed_d,
						bcMaskCollapsed_d, f_d, fColl_d, qCollapsed_d, bcCount);

				break;
			}
			gpuBcOutlet3D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					u_d, v_d, w_d, rho_d, bcCount);
			gpuBcPeriodic3D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					bcCount);
			gpuBcSymm3D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					bcCount);*/
		}
		// UPDATE VELOCITY AND DENSITY
		if(args->multiPhase){
			gpuUpdateMacro3DCG(nodeType, rho_d, u_d, v_d, w_d,
					bcBoundId_d, f_d, args->g,bcMask_d,args->UpdateInltOutl, r_f_d, b_f_d, r_rho_d, b_rho_d, p_in_d, p_out_d, num_in_d, num_out_d, args->test_case);
			switch(args->test_case){
			case 1:
				p_in_mean = gpu_sum_h(p_in_d, p_in_d, ms) / gpu_sum_int_h(num_in_d, num_in_d, ms);
				p_out_mean = gpu_sum_h(p_out_d, p_out_d, ms) / gpu_sum_int_h(num_out_d, num_out_d, ms);
				st_error[iter] = calculateSurfaceTension3D(p_in_mean, p_out_mean,args->r_alpha, args->b_alpha, args->bubble_radius * n, st_predicted);
				break;
			default:
				break;
			}
		}
		else{
			//gpuUpdateMacro3D(nodeType, rho_d, u_d, v_d, w_d,
			//		bcBoundId_d, nodeX, nodeY, nodeZ, f_d, args->g,bcMask_d,args->UpdateInltOutl);
		}
		tInstant2 = clock();
		// COMPUTE RESIDUALS

		if (AuxMacroDiff * args->ShowMacroDiff == iter + 1) {

			if (args->TypeOfResiduals == L2) {
				if(args->multiPhase){
				#pragma acc update device(f_d[0:19*ms], r_fColl_d[0:ms])
				}
				r = computeResidual3D(f_d, fColl_d, temp19a_d, temp19b_d, m, n, h);
			}
			else {
				if (args->TypeOfResiduals == FdRelDiff) {
					if(args->multiPhase){
					}

					if (firstIter) {
						firstIter = false;
#if OPENACC
                                #pragma acc host_data use_device(f_d, f1_d) 
                                {
                                cudaMemcpy(&f1_d,&f_d, ms*19,cudaMemcpyDeviceToDevice);
                                }
#else
                                delete[] f1_d;
                                f1_d = createHostArrayFlt(ms * 19, ARRAY_CPYD, 0, f_d);
#endif

					}
					r = computeNewResidual3D(f_d, fprev_d, f1_d, temp19a_d, temp19b_d, m, n, h);
					
#if OPENACC
                                #pragma acc host_data use_device(f_d, fprev_d) 
                                {
                                cudaMemcpy(&fprev_d,&f_d, ms*19,cudaMemcpyDeviceToDevice);
                                }
#else
                                delete[] fprev_d;
                                fprev_d = createHostArrayFlt(ms * 19, ARRAY_CPYD, 0, f_d);
#endif
				} else {
					bool d_divergence = false;
					if(args->multiPhase){
						fMaxDiff = gpu_abs_relSub_max(f_d, f_prev_d, n * m * h * 19);
					}
					else{/*
						gpu_abs_sub<<<bpg1, tpb>>>(u_d, u_prev_d, tempA_d,
								n * m * h, d_divergence);
						uMaxDiff = gpu_max_h(tempA_d, tempB_d, n * m * h);
						gpu_abs_sub<<<bpg1, tpb>>>(v_d, v_prev_d, tempA_d,
								n * m * h, d_divergence);
						vMaxDiff = gpu_max_h(tempA_d, tempB_d, n * m * h);
						gpu_abs_sub<<<bpg1, tpb>>>(w_d, w_prev_d, tempA_d,
								n * m * h, d_divergence);
						wMaxDiff = gpu_max_h(tempA_d, tempB_d, n * m * h);
						gpu_abs_sub<<<bpg1, tpb>>>(rho_d, rho_prev_d, tempA_d,
								n * m * h, d_divergence);
						rhoMaxDiff = gpu_max_h(tempA_d, tempB_d, n * m * h);*/
					}
					if (d_divergence) {
						fprintf(stderr, "\nDIVERGENCE!\n");
						break;
					}

					if(args->multiPhase){
						if(abs(fMaxDiff) < args->StopCondition[0]){
							printf("simulation converged!\n");
							break;
						}
					}
					else /*if (abs(uMaxDiff) < args->StopCondition[0] &&
							abs(vMaxDiff) < args->StopCondition[1] &&
							abs(wMaxDiff) < args->StopCondition[2] &&
							abs(rhoMaxDiff) < args->StopCondition[3]) {
						printf("simulation converged!\n");
						break;
					}*/

					if(args->multiPhase){
#if OPENACC
                                #pragma acc host_data use_device(f_d, f_prev_d) 
                                {
                                cudaMemcpy(&f_prev_d,&f_d, ms*19,cudaMemcpyDeviceToDevice);
                                }
#else
                                delete[] f_prev_d;
                                f_prev_d = createHostArrayFlt(ms * 19, ARRAY_CPYD, 0, f_d);
#endif
					}else{/*
						writeMacroDiffs(iter + 1, uMaxDiff, vMaxDiff, wMaxDiff,	rhoMaxDiff);
						CHECK(cudaFree(u_prev_d));
						CHECK(cudaFree(v_prev_d));
						CHECK(cudaFree(w_prev_d));
						CHECK(cudaFree(rho_prev_d));
						u_prev_d = createGpuArrayFlt(n * m * h, ARRAY_CPYD, 0, u_d);
						v_prev_d = createGpuArrayFlt(n * m * h, ARRAY_CPYD, 0, v_d);
						w_prev_d = createGpuArrayFlt(n * m * h, ARRAY_CPYD, 0, w_d);
						rho_prev_d = createGpuArrayFlt(n * m * h, ARRAY_CPYD, 0,
								rho_d);
				*/	}
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

			AuxMacroDiff++;

		}
		if(args->multiPhase){
#if OPENACC
                                #pragma acc host_data use_device(f_d, f_prev_d) 
                                {
                                cudaMemcpy(&f_prev_d,&f_d, ms*19,cudaMemcpyDeviceToDevice);
                                }
#else
                                delete[] f_prev_d;
                                f_prev_d = createHostArrayFlt(ms * 19, ARRAY_CPYD, 0, f_d);
#endif
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
				#pragma acc update host(u_d[0:ms], v_d[0:ms], w_d[0:ms], rho_d[0:ms])
				WriteResults3D(autosaveFilename, nodeType, nodeX, nodeY, nodeZ,
						u_d, v_d, w_d, rho_d, nodeType, n, m, h, args->outputFormat);
				tInstant2 = clock();
				taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1)/ CLOCKS_PER_SEC;
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
	#pragma acc update host(u_d[0:ms], v_d[0:ms], w_d[0:ms], rho_d[0:ms])
	if(args->multiPhase){
	#pragma acc update host(r_rho_d[0:ms], b_rho_d[0:ms])
	}

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
			deformingBubbleValid(r_rho_d, b_rho_d, n, m, h);
			break;
		case 3:
			validateCoalescenceCase(r_rho_d, b_rho_d, n, m, args->bubble_radius, h);
			break;
		case 4: //COUETTE
			analyticalCouette(args->kappa, nodeY, m, n, analytical, args->u, h);
			writeCouetteSolution("Profile_Couette", analytical, u_d, nodeY, m, n, h);
			printf("Couette profile written to Profile_Couette in Results/\n");
			break;
		default:
			break;
		}
		WriteResultsMultiPhase(finalFilename, nodeType, nodeX, nodeY, nodeZ, u_d, v_d, w_d, rho_d,r_rho_d,b_rho_d, nodeType,
				n, m, h, args->outputFormat);
	}
	else{
		WriteResults3D(finalFilename, nodeType, nodeX, nodeY, nodeZ, u_d, v_d, w_d, rho_d,
				nodeType, n, m, h, args->outputFormat);
	}

	WriteLidDrivenCavityMidLines3D(nodeX, nodeY, nodeZ, u_d, w_d, n, m, h, args->u);
	WriteChannelCrossSection3D(nodeX, nodeY, nodeZ, u_d, v_d, w_d, n, m, h, args->u);
	}
	// Write information for user
	printf("\n\nLog was written to %s\n", logFilename);
	printf("Last autosave result can be found at %s\n", autosaveFilename);
	printf("residuals were written to %s\n", residualsFilename);
	printf("Profiling results were written to %s\n", timeFilename);
	//	compareTestFiles("./TestValues/CUDA/cg2d.txt", "./TestValues/CUDA/cg3d.txt");
	printf("Final results were written to %s\n", finalFilename);

	freeAllHost(hostArrays, sizeof(hostArrays) / sizeof(hostArrays[0]));
	freeAllHost(gpuArrays, sizeof(gpuArrays) / sizeof(gpuArrays[0]));
	if(args->multiPhase){
		freeAllHost(mpGpuArrays, sizeof(mpGpuArrays) / sizeof(mpGpuArrays[0]));
		if(args->test_case == 1)
			freeAllHost(mpTC1GpuArrays, sizeof(mpTC1GpuArrays) / sizeof(mpTC1GpuArrays[0]));
		freeAllHost(mpHostArrays, sizeof(mpHostArrays) / sizeof(mpHostArrays[0]));
	}
	if (args->TypeOfResiduals == FdRelDiff) {
		freeAllHost(FDdifGpuArrays, sizeof(FDdifGpuArrays) / sizeof(FDdifGpuArrays[0]));
	}
	if (args->TypeOfResiduals != MacroDiff) {
		freeAllHost(nonMacroDiffGpuArrays, sizeof(nonMacroDiffGpuArrays) / sizeof(nonMacroDiffGpuArrays[0]));
	}
        if (args->boundaryType == (BoundaryType) CURVED) free(qCollapsed_d);
	return 0;
}
