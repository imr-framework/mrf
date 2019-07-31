# Define installation location for CUDA and compilation flags compatible
# with the CUDA include files.
CUDAHOME    = /usr/local/cuda
INCLUDEDIR  = -I$(CUDAHOME)/include
#INCLUDELIB  = -L$(CUDAHOME)/lib -lcufft -Wl,-rpath,$(CUDAHOME)/lib
INCLUDELIB  = -L$(CUDAHOME)/lib -lcudart -lcufft -Wl,-rpath,$(CUDAHOME)/lib
CFLAGS      = -fPIC -D_GNU_SOURCE -pthread -fexceptions
COPTIMFLAGS = -O3 -funroll-loops -msse2

# Define installation location for MATLAB.
export MATLAB = /usr/local/matlab
#export MATLAB = /Applications/MATLAB_R2007b
MEX           = $(MATLAB)/bin/mex
MEXEXT        = .$(shell $(MATLAB)/bin/mexext)

# nvmex is a modified mex script that knows how to handle CUDA .cu files.
#NVMEX = ./nvmex

# matlab-cuda install location
matcuda = /usr/local/matcuda

NVMEX = $(matcuda)/nvmex -f $(matcuda)/nvopts.sh

# List the mex files to be built.  The .mex extension will be replaced with the
# appropriate extension for this installation of MATLAB, e.g. .mexglx or
# .mexa64.
MEXFILES = test.mex

all: $(MEXFILES:.mex=$(MEXEXT))

clean:
	rm -f $(MEXFILES:.mex=$(MEXEXT))

.SUFFIXES: .cu .cu_o .mexglx .mexa64 .mexmaci

.cu.mexa64:
	$(NVMEX) $< $(INCLUDEDIR) $(INCLUDELIB)
