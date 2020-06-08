/**
 * Main file
 * @file main.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include <string.h>              // String operations
#include <stdio.h>

#include "Iterate.h"             // Iteration takes place
#include "ShellFunctions.h"      // For convenience
#include "FilesReading.h"        // For reading files
#include "FilesWriting.h"        // For writing files e.g. tecplot
#include "CellFunctions.h"       // For cell modifications
#include "ComputeResiduals.h"
#include "Arguments.h"
#include "errno.h"
#include "LogWriter.h"
/**
 * @brief Main function
 * @details
 *  - handles input parameters
 *  - creates directory for results
 *  - runs the solver
 *
 * @param argc number of command line parameters
 * @param argv command line parameters
 *
 * @return 0: success n: error
 */
int main(int argc, char* argv[]) {

	eraseAlertLog();
	Arguments args;
	args.u = 0.;
	args.v = 0.;
	args.w = 0.;
	args.rho = 0.;
	args.g=0.;
	args.viscosity = 0.;
	args.inletProfile = NO_INLET;
	args.outletProfile = OUTLET_SECOND;
	args.collisionModel = BGKW;
	args.boundaryType = STRAIGHT;
	args.outputFormat = PARAVIEW;
	args.iterations = 5000;
	args.autosaveEvery = 0;
	args.autosaveAfter = 0;
	args.boundaryId = 0;
	args.bcwallmodel = SIMPLE;
	args.StopCondition[0] = 0.0;
	args.StopCondition[1] = 0.0;
	args.StopCondition[2] = 0.0;
	args.StopCondition[3] = 0.0;
	args.UseInitialCondFromFile = 0;
	args.TypeOfResiduals = MacroDiff;
	args.ShowMacroDiff = 0;
	args.UpdateInltOutl=1;

	args.multiPhase = 1;
	if(args.multiPhase){
		//Read from file
		args.r_density = 1.0;
		args.gamma = 3.0;
		args.kappa = 1.0;
		args.b_alpha = 0.6;
		args.r_viscosity = 0.005;
		//args.b_viscosity = args.r_viscosity / args.kappa * args.gamma;
		args.b_viscosity = 0.005;
		args.beta = 0.99;
		args.A = 0.0;
		args.control_param = 0.9;
		args.g_limit = 0;
		args.bubble_radius = 9;
		args.test_case = 6;
		//Not file
		args.b_density = args.r_density / args.gamma;
		args.r_alpha = (1.0 - ((1.0 - args.b_alpha) / args.gamma));
		args.external_force = 0; //0 is gravity, 1 for pressure difference
		args.high_order = 0; // order of color gradient
		args.enhanced_distrib = 1; // 1 to use enhanced
	}

	InputFilenames inFn;
	strcpy(inFn.init, "SetUpData.ini");
	strcpy(inFn.result, "Results");
	strcpy(inFn.final, "Results/FinalData.csv");

	if (argc > 1) {
		switch (handleArguments(argc, argv, &inFn, &args)) {
		case INIT:
			readInitFile(inFn.init, &args);
			break;
		case HELP:
			return 0;
		case COMPARE:
			return compareFiles(inFn.comp, inFn.final);
		case ERROR:
			return 1;
		default:
			printf("NOTHING TO BE DONE in main.cu switch at line 55\n");
			break;
		}
	} else {
		readInitFile(inFn.init, &args);
	}
	snprintf(inFn.node, sizeof(inFn.node), "Mesh/%s_Node.dat",
			args.problemName);
	printf("%s\n", inFn.node);
	snprintf(inFn.bc, sizeof(inFn.bc), "Mesh/%s_BC.dat", args.problemName);
	printf("%s\n", inFn.bc);
	snprintf(inFn.InitialConditions, sizeof(inFn.InitialConditions), "%s",
			args.initialConditionsName);
	printf("%s\n", inFn.InitialConditions);

	CreateDirectory(inFn.result);

	//// 2D or 3D ?? ////
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	fp = fopen(inFn.node, "r");
	if (!fp) {
		fprintf(stderr, "Error reading file %s: %s\n", inFn.node,
				strerror(errno));
		return 0;
	}
	getline(&line, &len, fp);
	fclose(fp);

	int NoOfCol, i, j, k, m;
	FLOAT_TYPE f, g, h;
	NoOfCol = sscanf(line,
			"%d\t%d\t%d\t"FLOAT_FORMAT"\t"FLOAT_FORMAT"\t"FLOAT_FORMAT"\t%d",
			&i, &j, &k, &f, &g, &h, &m);

	FILE * f1;
	f1 = fopen("Results/MacroDiffs.txt", "w");
	fprintf(f1, "iteration\tMax u diff\tMax v diff\tMax w diff\tMax rho diff\n");
	fclose(f1);

	if (NoOfCol == 7) {
		args.TypeOfProblem=(ProblemType)_3D;
		printf("3D problem\n");
		Iterate3D(&inFn, &args);
	} else {
		args.TypeOfProblem=(ProblemType)_2D;
		printf("2D problem\n");
		Iterate2D(&inFn, &args);
	}

	//cudaDeviceReset(); //needed for cuda-gdb and cuda-memcheck

	return 0;
}

