FLOAT_TYPE=USE_DOUBLE

CFLAGS=-Llibs -largtable2 -DUSE_DOUBLE -Iinclude -Wall
#-ta=tesla:rdc -Mcuda -Minfo=accel -largtable2 -DUSE_DOUBLE -Iinclude -Llibs
DEFINES=-D$(FLOAT_TYPE)
INCLUDES=-Iinclude
OS := $(shell uname)
ifeq ($(OS),Darwin)
LDFLAGS=-arch=sm_35 -L/usr/local/cuda/lib -lcuda
CC=/Developer/NVIDIA/CUDA-7.5/bin/nvcc
LIBDIR=-LlibsMac
LIBS=-largtable2
else
LDFLAGS=#-ta=tesla
CC=g++
LIBDIR=-Llibs
LIBS=-largtable2
endif


CU_FILES=main.cpp Iterate.cpp CellFunctions.cpp FilesReading.cpp FilesWriting.cpp \
         ShellFunctions.cpp LogWriter.cpp GpuInit.cpp\
         ArrayUtils.cpp Arguments.cpp Multiphase.cpp
ITER_FILES=Iterate.cpp ComputeResiduals.cu\
           ArrayUtils.cpp Arguments.cpp
ITER_FILE=IterateCombined.cu
RLSE_FILES=main.cpp $(ITER_FILE) CellFunctions.cpp FilesReading.cpp FilesWriting.cpp \
           ShellFunctions.cpp LogWriter.cpp GpuInit.cpp
HEADRF=FilesReading.h FilesWriting.h Iterate.h CellFunctions.h ShellFunctions.h \
       LogWriter.h Arguments.h ArrayUtils.h Multiphase.h GpuFunctions.h GpuConstants.h argtable2.h
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

release: CFLAGS=-ta=tesla
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

