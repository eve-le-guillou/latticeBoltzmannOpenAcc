#include <stdio.h>   // for calloc();
#include "FilesWriting.h"
#include <math.h>

void WriteResultsMultiPhase(char* OutputFile, int* fluid, FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY,
		FLOAT_TYPE* CoordZ,
		FLOAT_TYPE* U, FLOAT_TYPE* V, FLOAT_TYPE* W, FLOAT_TYPE* Rho,FLOAT_TYPE* r_Rho, FLOAT_TYPE* b_Rho,
		int* Fluid, int n, int m, int h, OutputFormat outputFormat) {
	int i;                      // Loop variables
	FILE * fp1;                 // file pointer to output file
	fp1 = fopen(OutputFile, "w"); // open file
	switch (outputFormat) {

	case CSV:

		fprintf(fp1, "x,y,z,u,v,w,vel_mag,rho,press\n");
		for(i=0; i<m*n*h; ++i)
		{
			fprintf(fp1, "%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
					CoordX[i], // x
					CoordY[i], // y
					CoordZ[i], // z
					U[i],      // u
					V[i],      // v
					W[i],      // w
					sqrt(pow(U[i],2)+pow(V[i],2)+pow(W[i],2)), // u magnitude
					Rho[i],    // density
					Rho[i]/3);
		}
		fclose(fp1);
		break;
	case TECPLOT:
		fprintf(fp1, "Title = \"LBM results\"\n");
		fprintf(fp1,
				"Variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"vel mag\",\"rho\",\"press\"\n");
		fprintf(fp1, "Zone i=%d, j=%d, k=%d, f=point\n", n, m, h);

		for (i = 0; i < m * n * h; ++i) {
			//if (fluid[i] == 1){
			fprintf(fp1,
					"%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
					CoordX[i], // x
					CoordY[i], // y
					CoordZ[i], // z
					U[i],      // u
					V[i],      // v
					W[i],      // w
					sqrt(pow(U[i], 2) + pow(V[i], 2) + pow(W[i], 2)), // u magnitude
					Rho[i],    // density
					Rho[i] / 3);  // pressure
			//	}
		}


		fclose(fp1);
		break;
	case PARAVIEW:

		int x, y, z;
		fprintf(fp1, "<?xml version=\"1.0\"?>\n");
		fprintf(fp1, "<!-- 3D LBM Solver -->\n");
		fprintf(fp1,
				"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
		fprintf(fp1,
				"  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",
				n - 1, m - 1, h - 1);
		fprintf(fp1, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n", n - 1, m - 1,
				h - 1);
		fprintf(fp1, "    <PointData Scalars=\"scalars\">\n");

		//write U
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f %.10f %.10f\n", U[x + y * n + z * n * m],
							V[x + y * n + z * n * m], W[x + y * n + z * n * m]);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");
		//write density
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f\n", Rho[x + y * n + z * n * m]);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");

		//write red density
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Red Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f\n", r_Rho[x + y * n + z * n * m]);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");

		//write blue density
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Blue Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f\n", b_Rho[x + y * n + z * n * m]);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");

		//write pressure
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f\n", Rho[x + y * n + z * n * m]/3.0);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");




		fprintf(fp1, "    </PointData>\n");
		fprintf(fp1, "    <CellData>\n");
		fprintf(fp1, "    </CellData>\n");
		fprintf(fp1, "  </Piece>\n");
		fprintf(fp1, "  </ImageData>\n");
		fprintf(fp1, "</VTKFile>\n");
		fclose(fp1);
		break;
	}

}

void WriteResults3D(char* OutputFile, int* fluid, FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY,
		FLOAT_TYPE* CoordZ,
		FLOAT_TYPE* U, FLOAT_TYPE* V, FLOAT_TYPE* W, FLOAT_TYPE* Rho,
		int* Fluid, int n, int m, int h, OutputFormat outputFormat) {
	int i;                      // Loop variables
	FILE * fp1;                 // file pointer to output file
	fp1 = fopen(OutputFile, "w"); // open file
	switch (outputFormat) {

	case CSV:

		fprintf(fp1, "x,y,z,u,v,w,vel_mag,rho,press\n");
		for(i=0; i<m*n*h; ++i)
		{
			fprintf(fp1, "%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
					CoordX[i], // x
					CoordY[i], // y
					CoordZ[i], // z
					U[i],      // u
					V[i],      // v
					W[i],      // w
					sqrt(pow(U[i],2)+pow(V[i],2)+pow(W[i],2)), // u magnitude
					Rho[i],    // density
					Rho[i]/3);
		}
		fclose(fp1);
		break;
	case TECPLOT:
		fprintf(fp1, "Title = \"LBM results\"\n");
		fprintf(fp1,
				"Variables = \"x\",\"y\",\"z\",\"u\",\"v\",\"w\",\"vel mag\",\"rho\",\"press\"\n");
		fprintf(fp1, "Zone i=%d, j=%d, k=%d, f=point\n", n, m, h);

		for (i = 0; i < m * n * h; ++i) {
			//if (fluid[i] == 1){
			fprintf(fp1,
					"%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
					CoordX[i], // x
					CoordY[i], // y
					CoordZ[i], // z
					U[i],      // u
					V[i],      // v
					W[i],      // w
					sqrt(pow(U[i], 2) + pow(V[i], 2) + pow(W[i], 2)), // u magnitude
					Rho[i],    // density
					Rho[i] / 3);  // pressure
			//	}
		}


		fclose(fp1);
		break;
	case PARAVIEW:
		int x, y, z;
		fprintf(fp1, "<?xml version=\"1.0\"?>\n");
		fprintf(fp1, "<!-- 3D LBM Solver -->\n");
		fprintf(fp1,
				"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
		fprintf(fp1,
				"  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n",
				n - 1, m - 1, h - 1);
		fprintf(fp1, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n", n - 1, m - 1,
				h - 1);
		fprintf(fp1, "    <PointData Scalars=\"scalars\">\n");
		//write U
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f %.10f %.10f\n", U[x + y * n + z * n * m],
							V[x + y * n + z * n * m], W[x + y * n + z * n * m]);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");
		//write density
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f\n", Rho[x + y * n + z * n * m]);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");
		//write pressure
		fprintf(fp1,
				"      <DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		for (z = 0; z < h; z++) {
			for (y = 0; y < m; y++) {
				for (x = 0; x < n; x++) {
					fprintf(fp1, "%.10f\n", Rho[x + y * n + z * n * m]/3.0);
				}
			}
		}
		fprintf(fp1, "      </DataArray>\n");


		fprintf(fp1, "    </PointData>\n");
		fprintf(fp1, "    <CellData>\n");
		fprintf(fp1, "    </CellData>\n");
		fprintf(fp1, "  </Piece>\n");
		fprintf(fp1, "  </ImageData>\n");
		fprintf(fp1, "</VTKFile>\n");
		fclose(fp1);
		break;
	}

}

void WriteLidDrivenCavityMidLines3D(FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* CoordZ,
		FLOAT_TYPE* U, FLOAT_TYPE* W, int n, int m, int h, FLOAT_TYPE Ulid){

	FILE* linesFile;               // file for lines
	char linesFilename[768];       // path of the .lines file
	linesFilename[0] = '\0';
	strcat(linesFilename, "./Results/midLines.txt");
	linesFile = fopen(linesFilename, "w"); // open file
	fprintf(linesFile, "CASE: Ulid = %f Nx = %d\n", Ulid, n);
	fprintf(linesFile, "U/Ulid\tZ\tX\tW/Ulid\n");

	int x, y, z;
	FLOAT_TYPE tempUZ[n][2], tempXW[n][2];
	int tempUZ_idx, tempXW_idx, tempIdx;

	tempUZ_idx=0;
	tempXW_idx=0;

	for (z = 0; z < h; z++) {
		for (y = 0; y < m; y++) {
			for (x = 0; x < n; x++) {
				if(y==(int(m/2))){
					if(x==(int(n/2))){
						tempUZ[tempUZ_idx][0]=U[x + y * n + z * n * m]/Ulid;
						tempUZ[tempUZ_idx][1]=CoordZ[x + y * n + z * n * m];
						tempUZ_idx++;
					}
					if(z==(int(h/2))){
						tempXW[tempXW_idx][0]=CoordX[x + y * n + z * n * m];
						tempXW[tempXW_idx][1]=W[x + y * n + z * n * m]/Ulid;
						tempXW_idx++;
					}
				}
			}
		}
	}
	for(tempIdx = 0; tempIdx<tempUZ_idx;tempIdx++){
		fprintf(linesFile, "%.8f\t%.8f\t%.8f\t%.8f\n",tempUZ[tempIdx][0],tempUZ[tempIdx][1],tempXW[tempIdx][0],tempXW[tempIdx][1]);
	}
	fclose(linesFile);

	printf("\n\nLinesWritten to %s\n",linesFilename);
}

void WriteChannelCrossSection3D(FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* CoordZ,
		FLOAT_TYPE* U, FLOAT_TYPE* V, FLOAT_TYPE* W, int n, int m, int h, FLOAT_TYPE Ulid){

	FILE* sectionFile;               // file for lines
	char sectionFileName[768];       // path of the .lines file
	sectionFileName[0] = '\0';
	strcat(sectionFileName, "./Results/channelPlane.txt");
	sectionFile = fopen(sectionFileName, "w"); // open file
	fprintf(sectionFile, "CASE: Uin = %f Ny = %d\n", Ulid, m);
	fprintf(sectionFile, "X\tY\tZ\tU\tV\tW\n");

	int x, y, z;

	for (z = 0; z < h; z++) {
		for (y = 0; y < m; y++) {
			for (x = 0; x < n; x++) {
				if(x==(int((3./4.)*n))){
					fprintf(sectionFile,
							"%.10f %.10f %.10f %.10f %.10f %.10f\n",
							CoordX[x + y * n + z * n * m], // x
							CoordY[x + y * n + z * n * m], // y
							CoordZ[x + y * n + z * n * m], // z
							U[x + y * n + z * n * m],      // u
							V[x + y * n + z * n * m],      // v
							W[x + y * n + z * n * m]); 		// w
				}
			}
		}
	}
	fclose(sectionFile);

	printf("CrossSection to %s\n",sectionFileName);
}

void WriteArray(char* fileName, FLOAT_TYPE *arr, int n, int m, int h, int dir){

	char outputFile[300];
	sprintf(outputFile, "./Results/%s.txt", fileName);
	FILE * fp1;                 // file pointer to output file
	fp1 = fopen(outputFile, "w"); // open file
	int i, j, k;
	//	for(i = 0; i < n*m*h;i++)
	//		fprintf(fp1, "%.10f\n", arr[i]);
	for (k = 0; k < h; k++) {
		for (j = 0; j < m; j++) {
			for(i=0; i<n; i++){
				for(int l = 0; l < dir; l++){
					fprintf(fp1, "%.10f\n", arr[(i+j*n+k*n*m) * dir + l]);
				}
			}
		}
	}
	fclose(fp1);
}

void writeCouetteSolution(char* fileName, FLOAT_TYPE *analytical, FLOAT_TYPE *computed, FLOAT_TYPE *y, int m, int n, int h){
	char outputFile[300];
	sprintf(outputFile, "./Results/%s.txt", fileName);
	FILE * fp1;                 // file pointer to output file
	fp1 = fopen(outputFile, "w"); // open file
	fprintf(fp1, "y numerical analytical\n");
	int i, k;
	if(n % 2 == 0)
		i = n/2;
	else
		i = (n+1) / 2;

	if(h % 2 == 0)
		k = h / 2;
	else
		k = (h+1) / 2;

	int ms = n * m;
	for (int j = 0; j < m; j++) {
		fprintf(fp1, "%.4f %.10f %.10f\n", y[k * ms + j * n + i], computed[k * ms + j * n + i], analytical[j]);
	}

	fclose(fp1);
}

void writeOscilatingSolution(char* fileName, FLOAT_TYPE *extremes, int size){
	char outputFile[300];
	sprintf(outputFile, "./Results/%s.txt", fileName);
	FILE * fp1;                 // file pointer to output file
	fp1 = fopen(outputFile, "w"); // open file
	for (int i = 0; i < size; i++) {
		fprintf(fp1, "%d %.10f\n", i, extremes[i]);
	}

	fclose(fp1);
}
