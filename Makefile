FLOAT_TYPE=USE_DOUBLE

CFLAGS=-arch=sm_35 -rdc=true
PGIFLAGS=-x cu -ccbin pgc++ -Xcompiler " -ta=tesla:rdc -Mcuda -Minfo=accel -largtable2 -DUSE_DOUBLE -Iinclude -Llibs"
DEFINES=-D$(FLOAT_TYPE)
INCLUDES=-Iinclude
OS := $(shell uname)
ifeq ($(OS),Darwin)
LDFLAGS=-arch=sm_35 -L/usr/local/cuda/lib -lcuda
CC=/Developer/NVIDIA/CUDA-7.5/bin/nvcc
LIBDIR=-LlibsMac
LIBS=-largtable2
else
LDFLAGS=-arch=sm_35
CC=nvcc
LIBDIR=-Llibs
LIBS=-largtable2 -lcuda
endif

CPP_FILES=main.cpp Iterate.cpp CellFunctions.cpp ComputeResiduals.cpp FilesReading.cpp FilesWriting.cpp \
         ShellFunctions.cpp CpuInit.cpp CpuBoundaries.cpp CpuCollision.cpp CpuStream.cpp LogWriter.cpp \
         ArrayUtils.cpp Arguments.cpp CpuSum.cpp CpuUpdateMacro.cpp Multiphase.cpp

CU_FILES=Iterate3D.cu Check.cu 
ITER_FILES=Iterate.cpp CpuInit.cpp ComputeResiduals.cpp CpuBoundaries.cpp CpuCollision.cpp CpuStream.cpp \
           ArrayUtils.cpp Arguments.cpp CpuSum.cpp Check.cu CpuUpdateMacro.cpp
ITER_FILE=IterateCombined.cu
RLSE_FILES=main.cpp $(ITER_FILE) CellFunctions.cpp FilesReading.cpp FilesWriting.cpp \
           ShellFunctions.cpp LogWriter.cpp
HEADRF=FilesReading.h FilesWriting.h Iterate.h CellFunctions.h ComputeResiduals.h ShellFunctions.h \
       CpuFunctions.h LogWriter.h Arguments.h ArrayUtils.h CpuSum.h Multiphase.h
HEADERS=$(patsubst %,include/%,$(HEADRF))
ifeq ($(OS),Windows_NT)
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cu=.obj))
else
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cu=.o))+$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cpp=.o))
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

$(OBJ_DIR)/%.o: %.cu; \
	$(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: %.cpp; \
        $(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.obj: %.cpp; \
        $(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@


$(OBJ_DIR)/%.obj: %.cu; \
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

