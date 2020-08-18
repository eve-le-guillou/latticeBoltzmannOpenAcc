# OpenACC LBM SOLVER
## Cranfield University 2020
This software is a Lattice Boltzmann solver for simple geometries.

# Build
To build the program use the provided Makefile. By default it will create object files and build the code with double precision (float)

    $ make

To build with double precision

    $ make FLOAT_TYPE=USE_DOUBLE

To create a serial build:

    $ make serial

To create Doxygen documentation (HTML, LaTeX)

    $ make doc

To create Doxygen documentation and pdf

    $ made latexdoc

To create a parallel debug build that executes on Crescent's front node

    $ make debug

To clean the project

    $ make clean

#Submit to Crescent

To have the software executed on a GPU node of Crescent, use the gpu.sub file.

To have the software profiled on a GPU node of Crescent, use the profGpu.sub file.

# Dependencies
## Argtable 3
Argtable is an open source ANSI C library that parses GNU-style command-line options. It simplifies command-line parsing by defining a declarative-style API that can be used to specify what a program's command-line syntax should look like. Argtable3 is a single-file library. All that is needed to integrate it to the project is copying argtable3.c in the main directory and including argtable3.h in the source code. As PGC++ is a compiler for C++ source code, argtable3.c was modified to argtable3.cpp and sligthly differs from the original.

# Running the program 

## Compare files
To compare results in the default place with previous ones just run

    $ ./lbmsolver compare <path to result file>

To compare files at path different than the default

    $ ./lbmsolver compare <path to result file1> <path to result file2>

## Input files
The program needs at least 3 files to run if no parameters passed to the executable
 * SetUpData.ini
 * Mesh/<Lattice_name>_BC.dat
 * Mesh/<Lattice_name>_Node.dat

## Input parameters
The following was not maintained and may not work.

To get the full list of parameters

    $ ./lbmsolver -h

    Usage: ./lbmsolver  [-ht] [-f <file>] [-n <file>] [-b <file>] [-o <file>] [-u <u>] [-v <v>] [-r <rho>] [-s <nu>] [--inlet=[yes|no|pulsatile]] [-c [BGKW|TRT|MRT]] [--curved] [-l [yes|second|first]] [-i <N>] [--every=<N>] [--after=<N>] [--format=[paraview|tecplot]] [-d <id>]
    Usage: ./lbmsolver  compare <file> [<file>]

| short | long                          | description                                                 |
|-------|-------------------------------|-------------------------------------------------------------|
| -h    | --help                        | Print help options                                          |
| -f    | --initfile=\<file\>           | Initialisation from file (default: SetUpData.ini)           |
| -t    | --test                        | Run unit tests                                              |
| -n    | --node=\<file\>               | Node file (default: Mesh/D2node.dat)                        |
| -b    | --bc=\<file\>                 | Boundary conditions file (default: Mesh/BCconnectors.dat)   |
| -o    | --output=\<file\>             | Output directory (default: Results)                         |
| -u    | --uavg=\<u\>                  | Mean U (x velocity)                                         |
| -v    | --vavg=\<v\>                  | Mean V (y velocity)                                         |
| -r    | --rho=\<rho\>                 | Density                                                     |
| -s    | --viscosity=\<nu\>            | Viscosity                                                   |
|       | --inlet=[yes\|no\|pulsatile]  | Inlet profile (default: no)                                 |
| -c    | --collision=[BGKW\|TRT\|MRT]  | Collision model (default: BGKW)                             |
|       | --curved                      | Curved boundaries                                           |
| -l    | --outlet=[yes\|second\|first] | Outlet profile (default: second)                            |
| -i    | --iter=\<N\>                  | Number of iterations (default: 1000)                        |
|       | --every=\<N\>                 | Autosave after every \<N\> iterations (default: 0)          |
|       | --after=\<N\>                 | Start autosaving after the \<N\>th iteration (default: 1000)|
|       | --format=[paraview\|tecplot]  | Output format (default: paraview)                           |
| -d    | --draglift=\<id\>             | Calculate drag/lift on \<id\> boundary (default: 0)         |

If an init file is passed with -f all other parameters are discarded.

# Author

Eve Le Guillou
