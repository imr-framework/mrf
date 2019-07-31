# Define installation location for CUDA and compilation flags compatible
# with the CUDA include files.

CUDAHOME = /usr/local/cuda

#dir_mex_src = ../../
dir_mex_src = /n/ironwood/y/fessler/l/back/src/matlab/alg/mex/src/
dir_def = $(dir_mex_src)/def/
dir_fdk = $(dir_mex_src)/fdk/fdk-orig/

INCLUDEDIR = -I$(CUDAHOME)/include \
		-I$(dir_fdk) -I$(dir_def) \
		-DMmex -Dfdk_gpu -DUse_fdk_mex
#INCLUDELIB = -L$(CUDAHOME)/lib -lcufft -Wl,-rpath,$(CUDAHOME)/lib
INCLUDELIB = -L$(CUDAHOME)/lib -lcudart -lcufft -Wl,-rpath,$(CUDAHOME)/lib
#CFLAGS = -fPIC -D_GNU_SOURCE -pthread -fexceptions
COPTIMFLAGS = -O3 -funroll-loops -msse2

# Define installation location for MATLAB.
export MATLAB = /usr/local/matlab
MEX = $(MATLAB)/bin/mex
MEXEXT = .$(shell $(MATLAB)/bin/mexext)

# nvmex is a modified mex script that knows how to handle CUDA .cu files.
#NVMEX = ./nvmex

# matlab-cuda install location
matcuda = /usr/local/matcuda

NVMEX = $(matcuda)/nvmex -f $(matcuda)/nvopts.sh

# List the mex files to be built.  The .mex extension will be replaced with the
# appropriate extension for this installation of MATLAB, e.g. .mexglx or
# .mexa64.
MEXFILES = test.mex

goal1 = $(MEXFILES:.mex=$(MEXEXT))
goal = fdk_mex$(MEXEXT)
all: $(goal1) $(goal)

test.mex:	test.cu

fdk_h := $(dir_fdk)/*.h $(dir_def)/*.h

fdk_c = $(dir_fdk)/cbct,mask2.c \
	$(dir_fdk)/fdk,st.c \
	$(dir_def)/mexarg.c \
	$(dir_def)/jf,thread1.c

testlist = fdk-mex-cuda.cu fdk-ts.cu fdk-ts-t.cu fdk-ts-h.cu $(fdk_c)

test.cu:	*.h $(fdk_h) $(testlist)
		cat $(testlist) > $@

$(goal):	$(goal1)
		cp -p $? $@
clean:
	rm -f $(MEXFILES:.mex=$(MEXEXT))

.SUFFIXES: .cu .cu_o .mexglx .mexa64 .mexmaci

.cu.mexa64:
	$(NVMEX) $< $(INCLUDEDIR) $(INCLUDELIB)
