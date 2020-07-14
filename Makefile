FLOAT_TYPE=USE_DOUBLE

CFLAGS=-Minfo=accel -acc -ta=tesla:cc35 -fast
DEFINES=-D$(FLOAT_TYPE)
INCLUDES=-Iinclude
OS := $(shell uname)
LDFLAGS=-acc -ta=tesla:cc35
CC=pgc++
LIBDIR=-Mcuda=llvm 

CPP_FILES=argtable3.cpp main.cpp GpuSum.cpp GpuInit.cpp Iterate.cpp Iterate3D.cpp CellFunctions.cpp ComputeResiduals.cpp FilesReading.cpp FilesWriting.cpp \
         ShellFunctions.cpp GpuBoundaries.cpp GpuCollision.cpp GpuStream.cpp LogWriter.cpp \
         ArrayUtils.cpp Arguments.cpp GpuUpdateMacro.cpp Multiphase.cpp
HEADRF=FilesReading.h argtable3.h FilesWriting.h GpuConstants.h Iterate.h CellFunctions.h ComputeResiduals.h ShellFunctions.h \
       GpuFunctions.h LogWriter.h Arguments.h ArrayUtils.h GpuSum.h Multiphase.h
HEADERS=$(patsubst %,include/%,$(HEADRF))
ifeq ($(OS),Windows_NT)
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CPP_FILES:.cpp=.obj))
else
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CPP_FILES:.cpp=.o))
endif
OBJ_DIR=obj
DOC_DIR=docs

EXEC=lbmsolver

all: $(EXEC)
$(EXEC): $(OBJ_DIR) $(OBJ_FILES); \
	$(CC) $(LDFLAGS) $(LIBDIR) -o $@ $(OBJ_FILES)

debug: CFLAGS= -g -Minfo=accel -acc -ta=tesla:cc60
debug: LDFLAGS=-acc -ta=tesla:cc60
#debug: DEFINES+=-DMAKE_SERIAL
debug: all

serial: CFLAGS=-fast
serial: LDFLAGS=
serial: DEFINES+=-DMAKE_SERIAL
serial: all

$(OBJ_DIR):; \
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.cpp; \
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.obj: %.cpp; \
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(DOC_DIR):; \
	mkdir -p $(DOC_DIR)

doc: include/*.h Doxyfile $(DOC_DIR) README.md; \
	doxygen Doxyfile

latexdoc: doc; \
	make -C $(DOC_DIR)/latex

.PHONY: clean

clean:; \
	rm -rf $(OBJ_DIR) $(EXEC) $(RLSE) $(ITER_FILE) $(DOC_DIR)

