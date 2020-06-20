FLOAT_TYPE=USE_DOUBLE

CFLAGS=-fast -Minfo=accel -acc -ta=tesla:cc35 #-ta=tesla:cc60
DEFINES=-D$(FLOAT_TYPE)
INCLUDES=-Iinclude
OS := $(shell uname)
LDFLAGS=-acc -ta=tesla:cc35
CC=pgc++
LIBDIR=#-Mcuda=llvm #-Llibs
LIBS=#-largtable2 #-lcuda

CU_FILES=argtable3.cpp main.cpp GpuSum.cpp GpuInit.cpp Iterate.cpp Iterate3D.cpp CellFunctions.cpp ComputeResiduals.cpp FilesReading.cpp FilesWriting.cpp \
         ShellFunctions.cpp GpuBoundaries.cpp GpuCollision.cpp GpuStream.cpp LogWriter.cpp \
         ArrayUtils.cpp Arguments.cpp GpuUpdateMacro.cpp Multiphase.cpp
ITER_FILES=Iterate.cpp GpuInit.cpp ComputeResiduals.cpp GpuBoundaries.cpp GpuCollision.cpp GpuStream.cpp \
           ArrayUtils.cpp Arguments.cpp GpuSum.cpp GpuUpdateMacro.cpp
ITER_FILE=IterateCombined.cpp
RLSE_FILES=main.cpp $(ITER_FILE) CellFunctions.cpp FilesReading.cpp FilesWriting.cpp \
           ShellFunctions.cpp LogWriter.cpp
HEADRF=FilesReading.h argtable3.h FilesWriting.h GpuConstants.h Iterate.h CellFunctions.h ComputeResiduals.h ShellFunctions.h \
       GpuFunctions.h LogWriter.h Arguments.h ArrayUtils.h GpuSum.h Multiphase.h
HEADERS=$(patsubst %,include/%,$(HEADRF))
ifeq ($(OS),Windows_NT)
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cpp=.obj))
else
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cpp=.o))
endif
OBJ_DIR=obj
DOC_DIR=docs

EXEC=lbmsolver
RLSE=lbmrel

all: $(EXEC)
$(EXEC): $(OBJ_DIR) $(OBJ_FILES); \
	$(CC) $(LDFLAGS) $(LIBDIR) $(LIBS) -o $@ $(OBJ_FILES)

debug: CFLAGS+= -g -G -lineinfo
debug: DEFINES+= -DDEBUG
debug: all

release: CFLAGS=-arch=sm_35 -O3
release: DEFINES+= -DRELEASE
release: $(RLSE)

$(RLSE): $(RLSE_FILES); \
	$(CC) $(CFLAGS) $(LIBDIR) $(LIBS) $(DEFINES) $(INCLUDES) -o $@ $(RLSE_FILES)

$(ITER_FILE): $(ITER_FILES); \
	./combine.sh $(ITER_FILE) $(ITER_FILES)

$(OBJ_DIR):; \
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.cpp; \
	$(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.obj: %.cpp; \
	$(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(DOC_DIR):; \
	mkdir -p $(DOC_DIR)

doc: include/*.h Doxyfile $(DOC_DIR) README.md; \
	doxygen Doxyfile

latexdoc: doc; \
	make -C $(DOC_DIR)/latex

.PHONY: clean

clean:; \
	rm -rf $(OBJ_DIR) $(EXEC) $(RLSE) $(ITER_FILE) $(DOC_DIR)

