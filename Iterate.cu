#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <errno.h>
#include <cstring>
#include <cstdlib>
#include "GpuFunctions.h"       // GPU kernels
#include "ShellFunctions.h"     // For convenience
#include "FilesReading.h"       // For reading files
#include "FilesWriting.h"       // For writing files e.g. tecplot
#include "CellFunctions.h"      // For cell modifications
#include "ComputeResiduals.h"   // residuals
#include "LogWriter.h"
#include "Iterate.h"
#include "ArrayUtils.h"
#include "Check.h"
#include "Multiphase.h"
#include "GpuSum.h"

#define CUDA 1

int Iterate2D(InputFilenames *inFn, Arguments *args) {
	// Time measurement: declaration, begin
	clock_t tStart = clock();

	FILE* logFile;               // file for log
	char autosaveFilename[768];  // autosave filename
	char outputFilename[768];    // initial data will be written to this file
	char finalFilename[768];     // final data will be written to this file
	char logFilename[768];       // path of the .log file
	char residualsFilename[768]; // path of the residuals file
	char timeFilename[768];      // path of time measurement file

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
	int n, m;                   // number of nodes in the x and y directions
	FLOAT_TYPE maxInletCoordY; // maximum inlet coordinate in y
	FLOAT_TYPE minInletCoordY; // minimum inlet coordinate in y
	int numInletNodes;         // number of inlet nodes

	int AuxMacroDiff = 1;

	int *nodeIdX, *nodeIdY, *nodeType, *bcNodeIdX, *bcNodeIdY, *latticeId,
	*bcType, *bcBoundId,*tempi;
	FLOAT_TYPE *nodeX, *nodeY, *bcX, *bcY,*temp;

	FLOAT_TYPE taskTime[9];
	int i;
	for (i = 0; i < 9; ++i) {
		taskTime[i] = 0.0;
	}

	clock_t tInstant1, tInstant2; // Time measurement points, universal
	clock_t tIterStart, tIterEnd; // Time measurement points: main loop

	numNodes = readNodeFile(inFn->node, &nodeIdX, &nodeIdY, &tempi, &nodeX,
			&nodeY, &temp, &nodeType,args->TypeOfProblem);

	if (numNodes == 0) {
		printf("NODES NOT FOUND in file\n");
		return 2;
	}

	numConns = readConnFile(inFn->bc, &bcNodeIdX, &bcNodeIdY, &tempi,
			&latticeId, &bcType, &bcX, &bcY, &temp, &bcBoundId,args->TypeOfProblem);

	if (numConns == 0) {
		printf("NEIGHBOURING NOT FOUND in file\n");
		return 2;
	}

	if(!args->multiPhase || args->test_case != 6){
		m = getLastValue(nodeIdY, numNodes);
		n = getLastValue(nodeIdX, numNodes);

		delta = getGridSpacing(nodeIdX, nodeIdY, nodeX, numNodes);
		numInletNodes = getNumInletNodes(bcType, latticeId, numConns,
				args->TypeOfProblem);
		maxInletCoordY = getMaxInletCoordY(bcType, latticeId, bcY, delta, numConns,args->TypeOfProblem);
		minInletCoordY = getMinInletCoordY(bcType, latticeId, bcY, delta, numConns,args->TypeOfProblem);
	}
	else{
		m = 1024;
		n = 256;
		free(nodeX);
		free(nodeY);
		nodeX = createHostArrayFlt(n*m,ARRAY_ZERO);
		nodeY = createHostArrayFlt(n*m,ARRAY_ZERO);
		FLOAT_TYPE deltaX = 1.0 / (n-1), deltaY = (4.0 / (m-1));
		for(int j = 0; j < m; j++){
			for(int i = 0; i < n; i++){
				nodeX[j * n + i] = i * deltaX;
				nodeY[j * n + i] = j * deltaY;
			}
		}
	}

	FLOAT_TYPE *nodeZ = createHostArrayFlt(m * n, ARRAY_ZERO);
	writeInitLog(logFilename, args, delta, m, n, 1, numInletNodes, maxInletCoordY, minInletCoordY, 0.0, 0.0);
	logFile = fopen(logFilename, "a");

	// In case of no autosave
	sprintf(autosaveFilename, "NOWHERE!");
	initConstants2D(args, maxInletCoordY, minInletCoordY, delta, m, n);

	// residuals
	FLOAT_TYPE *norm = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *dragSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE *liftSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);



	fprintf(logFile, "\n:::: Initializing ::::\n");
	printf("\n:::: Initializing ::::\n");
	FLOAT_TYPE *u = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *v = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *w = createHostArrayFlt(m * n, ARRAY_ZERO);

	//Multiphase
	FLOAT_TYPE *st_error = createHostArrayFlt(args->iterations, ARRAY_ZERO);
	FLOAT_TYPE aux1 = args->r_density / ((args->r_density + args->b_density) * args->r_viscosity) +
			args->b_density / ((args->r_density + args->b_density) * args->b_viscosity);
	FLOAT_TYPE mean_nu = 1.0/aux1;
	FLOAT_TYPE omega_eff = 1.0/(3.0*mean_nu+0.5);

	FLOAT_TYPE st_predicted = 4.0 * args->A / 9.0 / omega_eff;

	int *cg_directions = createHostArrayInt(n*m, ARRAY_ZERO);

	FLOAT_TYPE *color_gradient = createHostArrayFlt(m * n * 2, ARRAY_ZERO);

	FLOAT_TYPE *r_f = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_f = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *r_fColl = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_fColl = createHostArrayFlt(m * n * 9, ARRAY_ZERO);


	FLOAT_TYPE *rho = createHostArrayFlt(m * n, ARRAY_FILL, args->rho);

	FLOAT_TYPE *r_rho = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *b_rho = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *r_f_d = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_f_d = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *r_fColl_d = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	FLOAT_TYPE *b_fColl_d = createHostArrayFlt(m * n * 9, ARRAY_ZERO);
	int *cg_dir_d = createHostArrayInt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *p_in_d = createHostArrayFlt(n*m, ARRAY_ZERO);
	FLOAT_TYPE *p_out_d = createHostArrayFlt(n*m, ARRAY_ZERO);
	int *num_in_d = createHostArrayInt(n*m, ARRAY_ZERO);
	int *num_out_d = createHostArrayInt(n*m, ARRAY_ZERO);
	FLOAT_TYPE *oscilating_y = createHostArrayFlt(args->iterations, ARRAY_NONE);
	FLOAT_TYPE *u0_d, *v0_d;
	FLOAT_TYPE fMaxDiff = -1;
	bool h_divergence;
	bool *d_divergence;
	//d_divergence = malloc(sizeof(bool));

	if (args->inletProfile == NO_INLET) {
		u0_d = createHostArrayFlt(m * n, ARRAY_FILL, args->u);
		v0_d = createHostArrayFlt(m * n, ARRAY_FILL, args->v);
	} else {
		u0_d = createHostArrayFlt(m * n, ARRAY_ZERO);
		v0_d = createHostArrayFlt(m * n, ARRAY_ZERO);
	}
	if (args->inletProfile == INLET) {
		//gpuInitInletProfile2D<<<bpg1, tpb>>>(u0_d, v0_d, nodeY, m * n);
	}
	FLOAT_TYPE *drag_d = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *lift_d = createHostArrayFlt(m * n, ARRAY_ZERO);

	FLOAT_TYPE *f_d = createHostArrayFlt(9 * m * n, ARRAY_ZERO);
	FLOAT_TYPE *fColl_d = createHostArrayFlt(9 * m * n, ARRAY_ZERO);

	FLOAT_TYPE *temp9a_d = createHostArrayFlt(9 * m * n, ARRAY_ZERO);
	FLOAT_TYPE *temp9b_d = createHostArrayFlt(9 * m * n, ARRAY_ZERO);
	FLOAT_TYPE *tempA_d = createHostArrayFlt(m * n, ARRAY_ZERO);
	FLOAT_TYPE *tempB_d = createHostArrayFlt(m * n, ARRAY_ZERO);


	FLOAT_TYPE p_in_mean;
	FLOAT_TYPE p_out_mean;
	FLOAT_TYPE ms = n * m;
	if(args->multiPhase){
		if(args->high_order)
			initHOColorGradient(cg_directions, n, m);
		else
			initColorGradient(cg_directions, n, m);
		memcpy(cg_dir_d, cg_directions, SIZEINT(m*n));
		initCGBubble(nodeX,nodeY,r_rho, b_rho, rho, r_f_d, b_f_d, f_d, args->test_case);
	}

	FLOAT_TYPE *f_prev_d = createHostArrayFlt(9 * m * n, ARRAY_COPY,0,r_f_d);
	int *mask = createHostArrayInt(m * n, ARRAY_ZERO);
	int *bcMask = createHostArrayInt(m * n, ARRAY_ZERO);
	int *bcIdx = createHostArrayInt(m * n, ARRAY_ZERO);

	FLOAT_TYPE *u_d = createHostArrayFlt(m * n, ARRAY_CPYD, 0, u0_d);
	FLOAT_TYPE *v_d = createHostArrayFlt(m * n, ARRAY_CPYD, 0, v0_d);
	int *stream = createHostArrayInt(8 * m * n, ARRAY_FILL, 1);
	FLOAT_TYPE *q = createHostArrayFlt(8 * m * n, ARRAY_FILL, 0.5);

	int bcCount = initBoundaryConditions2D(bcNodeIdX, bcNodeIdY, q, bcBoundId,
			nodeType, bcX, bcY, nodeX, nodeY, latticeId, stream, bcType, bcMask,
			bcIdx, mask, delta, m, n, numConns);
	int *bcIdxCollapsed_d = createHostArrayInt(bcCount, ARRAY_ZERO);
	int *bcMaskCollapsed_d = createHostArrayInt(bcCount, ARRAY_ZERO);
	FLOAT_TYPE *qCollapsed_d = createHostArrayFlt(8 * bcCount, ARRAY_ZERO);

	dim3 bpgB((int) (bcCount / THREADS) + 1); // blocks/grid
	int *bcMask_d = createHostArrayInt(m * n, ARRAY_COPY, 0, bcMask);
	int *bcIdx_d = createHostArrayInt(m * n, ARRAY_COPY, 0, bcIdx);

	collapseBc2D(bcIdx, bcIdxCollapsed_d, bcMask, bcMaskCollapsed_d, q,
			qCollapsed_d, mask, m, n, bcCount);

	int *stream_d = createHostArrayInt(8 * m * n, ARRAY_COPY, 0, stream);



	if(args->multiPhase && args->test_case == 2) //only for couette
	{
		initInletVelocity(u, v, args->u, args->v, n, m);
		CHECK(cudaMemcpy(u_d, u, SIZEFLT(m*n), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(v_d, v, SIZEFLT(m*n), cudaMemcpyHostToDevice));
	}
	memcpy(u, u_d, SIZEFLT(m*n));
	memcpy(v, v_d, SIZEFLT(m*n));

	fclose(logFile);
	writeNodeNumbers(logFilename, numNodes, numConns, bcCount);
	logFile = fopen(logFilename, "a");


	void *hostArrays[] = { nodeIdX, nodeIdY, nodeX, nodeY, nodeType, bcNodeIdX,
			bcNodeIdY, latticeId, bcType, bcX, bcY, bcBoundId, u, v, rho, mask,
			bcMask, bcIdx, stream, q, norm, dragSum, liftSum, r_rho, b_rho,
			color_gradient,r_f,b_f,r_fColl, b_fColl, st_error, cg_directions};
	void *gpuArrays[] = { u_d, v_d,
			u0_d, v0_d, drag_d, lift_d, f_d, fColl_d, temp9a_d, temp9b_d,
			tempA_d, tempB_d, bcMask_d, bcMaskCollapsed_d, bcIdx_d,
			bcIdxCollapsed_d, stream_d,  qCollapsed_d, r_f_d, r_fColl_d, b_f_d,
			b_fColl_d, cg_dir_d, p_in_d, p_out_d, num_in_d, num_out_d};
	fprintf(logFile, "\n:::: Initialization done! ::::\n");

	printf("Initialization took %f seconds\n", taskTime[T_INIT]);

	// Write Initialized data
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
	size_t freeSpace, total;

	/*cuMemGetInfo(&freeSpace, &total);
	printf("^^^^ Free : %llu Mbytes \n",
			(unsigned long long) freeSpace / 1024 / 1024);

	printf("^^^^ Total: %llu Mbytes \n",
			(unsigned long long) total / 1024 / 1024);

	printf("^^^^ %f%% free, %f%% used\n", 100.0 * freeSpace / (double) total,
			100.0 * (total - freeSpace) / (double) total);*/

	tInstant1 = clock(); // Start measuring time
	if(args->multiPhase){
		WriteResultsMultiPhase(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,r_rho,b_rho, nodeType,
				n, m, 1, args->outputFormat);
	}
	else
		WriteResults3D(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho, nodeType,
				n, m, 1, args->outputFormat);
	tInstant2 = clock();
	taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1) / CLOCKS_PER_SEC;
	printf("\nInitialized data was written to %s\n", outputFilename);


	////////////////// ITERATION ///////////////////////

	fprintf(logFile, "\n:::: Start Iterations ::::\n");
	printf("\n:::: Start Iterations ::::\n");

	printf("%d is the number of iterations \n", args->iterations);

	tIterStart = clock(); // Start measuring time of main loop
	int iter = 0;
	while (iter < args->iterations) {
		/*CHECK(cudaThreadSynchronize());
		CHECK(cudaEventRecord(start, 0)); // Start measuring time*/
		switch (args->collisionModel) {

		case BGKW:
			if(args->multiPhase){
				//Collision
				if(!args->enhanced_distrib)
					gpuCollBgkwGC2D(rho, r_rho, b_rho, u_d, v_d, f_d, r_fColl_d, b_fColl_d, cg_dir_d, args->high_order);
				else
					gpuCollEnhancedBgkwGC2D(rho, r_rho, b_rho, u_d, v_d, f_d, r_fColl_d, b_fColl_d, cg_dir_d, args->high_order);
			}else{
				//gpuCollBgkw2D<<<bpg1, tpb>>>(nodeType, rho, u_d, v_d, f_d,
				//		fColl_d);
			}
			break;
		case TRT:
			//gpuCollTrt<<<bpg1, tpb>>>(nodeType, rho, u_d, v_d, f_d, fColl_d);
			break;

		case MRT:
			//gpuCollMrt2D<<<bpg1, tpb>>>(nodeType, rho, u_d, v_d, f_d, fColl_d);
			break;
		}

		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_COLL] += cudatime;*/

		////////////// STREAMING ///////////////
		//CHECK(cudaThreadSynchronize());
		//CHECK(cudaEventRecord(start, 0));
		if(args->multiPhase){
			//			gpuStreaming2D<<<bpg1, tpb>>>(nodeType, stream_d, r_f_d, r_fColl_d);
			//			gpuStreaming2D<<<bpg1, tpb>>>(nodeType, stream_d, b_f_d, b_fColl_d);
			gpuStreaming2DCG(nodeType, stream_d, r_f_d, r_fColl_d, b_f_d, b_fColl_d, cg_dir_d);
		}
		else{
			//gpuStreaming2D<<<bpg1, tpb>>>(nodeType, stream_d, f_d, fColl_d);
		}
		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_STRM] += cudatime;*/

		// make the host block until the device is finished with foo
		//CHECK(cudaThreadSynchronize());
		// check for error
		/*cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess) {
			// print the CUDA error message and exit
			printf("CUDA error: %s\n", cudaGetErrorString(error));
			exit(-1);
		}*/
		////////////// BOUNDARIES ///////////////
		//CHECK(cudaThreadSynchronize());
		//CHECK(cudaEventRecord(start, 0));

		if(args->multiPhase){
			gpuBcPeriodic2D(bcIdxCollapsed_d, bcMaskCollapsed_d, r_f_d, b_f_d,bcCount, cg_dir_d, args->test_case, r_rho, b_rho, rho,
					u_d, v_d);

		} else{
		/*	gpuBcInlet2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
				u0_d, v0_d, bcCount);
			gpuBcWall2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					fColl_d, qCollapsed_d, bcCount);
			gpuBcOutlet2D<<<bpgB, tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d,
					u0_d, v0_d, bcCount);*/
		}

		/*CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_BNDC] += cudatime;*/

		// UPDATE VELOCITY AND DENSITY
		//CHECK(cudaThreadSynchronize());
		//CHECK(cudaEventRecord(start, 0));
		if(args->multiPhase){
			gpuUpdateMacro2DCG(rho, u_d, v_d, r_f_d, b_f_d, f_d, r_rho, b_rho, p_in_d, p_out_d, num_in_d, num_out_d, cg_dir_d,
					args->test_case);

			//			updateSurfaceTension(r_rho,b_rho,args->control_param, st_predicted, st_error, iter,args->r_alpha, args->b_alpha, args->bubble_radius, n ,m);
			//gpu reduction is faster than serial surface tension
			switch(args->test_case){
			case 1:
				p_in_mean = gpu_sum_h(p_in_d, p_in_d, ms) / gpu_sum_int_h(num_in_d, num_in_d, ms);
				p_out_mean = gpu_sum_h(p_out_d, p_out_d, ms) / gpu_sum_int_h(num_out_d, num_out_d, ms);
				st_error[iter] = calculateSurfaceTension(p_in_mean, p_out_mean,args->r_alpha, args->b_alpha, args->bubble_radius * n, st_predicted);
				break;
			case 5:
				oscilating_y[iter] = getMaxYOscilating(r_rho, b_rho, n, m, nodeY);
				break;
			case 6:
				oscilating_y[iter] = getMinYRT(r_rho, b_rho, n, m, nodeY);
				break;
			default:
				break;
			}
		}
		else; //gpuUpdateMacro2D<<<bpg1, tpb>>>(nodeType, rho, u_d, v_d, bcMask_d,
			//	drag_d, lift_d, nodeX, nodeY, f_d);

		/*tInstant2 = clock();
		CHECK(cudaEventRecord(stop, 0));
		CHECK(cudaEventSynchronize(stop));
		CHECK(cudaEventElapsedTime(&cudatime, start, stop));
		taskTime[T_MACR] += cudatime;*/

		// COMPUTE RESIDUALS
		if (AuxMacroDiff * args->ShowMacroDiff == iter + 1) {
			//CHECK(cudaThreadSynchronize());
			//CHECK(cudaEventRecord(start, 0));
			FLOAT_TYPE r;
			if(args->multiPhase){
				//				gpu_abs_sub(f_d, f_prev_d, temp9a_d, n * m * 9, d_divergence);
				gpu_abs_relSub(f_d, f_prev_d, temp9a_d, n * m * 9, d_divergence);
				fMaxDiff = gpu_max_h(temp9a_d, temp9b_d, n * m * 9);
				//	printf("MAX diff "FLOAT_FORMAT"\n", fMaxDiff);
				memcpy(&h_divergence,d_divergence,sizeof(bool));
				if (h_divergence || fMaxDiff != fMaxDiff || !isfinite(fMaxDiff)) {
					fprintf(stderr, "\nDIVERGENCE!\n");
					break;
				}
				if(abs(fMaxDiff) < args->StopCondition[0]){
					printf("simulation converged!\n");
					break;
				}
				delete[] f_prev_d;
				f_prev_d = createHostArrayFlt(n * m * 9, ARRAY_CPYD, 0, f_d);

			}else{
				r = computeResidual2D(f_d, fColl_d, temp9a_d, temp9b_d, m,n);
				if (r != r) {
					fprintf(stderr, "\nDIVERGENCE!\n");

					writeResiduals(residualsFilename, norm, dragSum, liftSum, m * n,
							iter + 1);
				//	cudaEventDestroy(start);
				//	cudaEventDestroy(stop);

					freeAllHost(hostArrays, sizeof(hostArrays) / sizeof(hostArrays[0]));
					freeAllGpu(gpuArrays, sizeof(gpuArrays) / sizeof(gpuArrays[0]));

					return 1; // ERROR!
				}
				norm[iter] = r;
				if (args->boundaryId > 0) {
					dragSum[iter] = computeDragLift2D(bcMask_d, drag_d, tempA_d,
							tempB_d, m, n, args->boundaryId);
					liftSum[iter] = computeDragLift2D(bcMask_d, lift_d, tempA_d,
							tempB_d, m, n, args->boundaryId);
				}
			}

			/*CHECK(cudaEventRecord(stop, 0));
			CHECK(cudaEventSynchronize(stop));
			CHECK(cudaEventElapsedTime(&cudatime, start, stop));
			taskTime[T_RESI] += cudatime;*/
			AuxMacroDiff++;

		}
		printf("Iterating... %d/%d (%3.1f %%) Residual Max %.15f\r", iter + 1, args->iterations,
				(FLOAT_TYPE) (iter + 1) * 100
				/ (FLOAT_TYPE) (args->iterations), fMaxDiff);

		iter++; // update loop variable
		////////////// Autosave ///////////////

		if (iter == (args->autosaveEvery * autosaveIt)) {
			autosaveIt++;
			if (iter > args->autosaveAfter) {
				printf("autosave\n\n");
				//////////// COPY VARIABLES TO HOST ////////////////
				memcpy(u, u_d, SIZEFLT(m*n));
				memcpy(v, v_d, SIZEFLT(m*n));
				switch (args->outputFormat) {
				case CSV:
					sprintf(finalFilename, "%sFinalData.csv", inFn->result);
					break;
				case TECPLOT:
					sprintf(finalFilename, "%sFinalData.dat", inFn->result);
					break;
				case PARAVIEW:
					sprintf(finalFilename, "%s/Video/FinalData_%d.vti", inFn->result, autosaveIt);
					break;
				}
				tInstant1 = clock(); // Start measuring time
				WriteResults3D(finalFilename, nodeType,nodeX, nodeY, nodeZ, u, v, w, rho,
						nodeType, n, m, 1, args->outputFormat);
				tInstant2 = clock();
				taskTime[T_WRIT] += (FLOAT_TYPE) (tInstant2 - tInstant1);
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
	writeResiduals(residualsFilename, norm, dragSum, liftSum, m * n,
			args->iterations);

	//WRITE VARIABLES TO HOST
	memcpy(u, u_d, SIZEFLT(m*n));
	memcpy(v, v_d, SIZEFLT(m*n));
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
			printf("Pressure inside "FLOAT_FORMAT" Pressure outside "FLOAT_FORMAT" ST_predicted "FLOAT_FORMAT"\n", p_in_mean, p_out_mean, st_predicted);
			WriteArray("surface tension",st_error, args->iterations,1);
			break;
		case 2:
			if(args->g == 0.0){
				analyticalCouette(args->kappa, nodeY, m, n, analytical, args->u);
				writeCouetteSolution("Profile_Couette", analytical, u, nodeY, m, n);
				printf("Couette profile written to Profile_Couette in Results/\n");
			}
			else{
				analyticalPoiseuille(m, n, analytical, args->r_density, args->b_density, args->r_viscosity, args->b_viscosity, args->g, nodeY);
				writeCouetteSolution("Profile_Poiseuille", analytical, u, nodeY, m, n);
				printf("Poiseuille profile written to Profile_Poiseuille in Results/\n");
			}
			break;
		case 3:
			deformingBubbleValid(r_rho, b_rho, n, m);
			break;
		case 4:
			//validate coalescence case
			validateCoalescenceCase(r_rho, b_rho, n, m, args->bubble_radius);
			break;
		case 5:
			writeOscilatingSolution("Oscilating_frequency", oscilating_y, args->iterations);
			printf("Oscilation frequency written to Oscilating_frequency in Results/\n");
			printf("Error % : "FLOAT_FORMAT"\n", validateOscilating(r_rho, b_rho, n, m, oscilating_y, args->iterations,st_predicted, args->r_density, args->b_density));
			break;
		case 6:
			writeOscilatingSolution("RT_Interface", oscilating_y, args->iterations);
			break;
		default:
			printf("Suface tension error: "FLOAT_FORMAT"\n", st_error[iter-1]);
			break;
		}
		WriteResultsMultiPhase(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho,r_rho,b_rho, nodeType,
				n, m, 1, args->outputFormat);

		free(analytical);
	}
	else
		WriteResults3D(finalFilename, nodeType, nodeX, nodeY, nodeZ, u, v, w, rho, nodeType,
				n, m, 1, args->outputFormat);

	// Write information for user
	printf("\n\nLog was written to %s\n", logFilename);
	printf("Last autosave result can be found at %s\n", autosaveFilename);
	printf("residuals were written to %s\n", residualsFilename);
	printf("Profiling results were written to %s\n", timeFilename);
	printf("Final results were written to %s\n", finalFilename);

	//	compareTestFiles("./TestValues/CUDA/rpert.txt", "./TestValues/CUDA/rpert_gpu.txt");
	//cudaEventDestroy(start);
	//cudaEventDestroy(stop);
	freeAllHost(hostArrays, sizeof(hostArrays) / sizeof(hostArrays[0]));
	freeAllHost(gpuArrays, sizeof(gpuArrays) / sizeof(gpuArrays[0]));
	return 0;
}
