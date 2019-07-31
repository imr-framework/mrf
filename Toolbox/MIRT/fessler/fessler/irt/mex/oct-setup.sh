# oct-setup.sh
# extra flags needed for compiling mex files for octave

version = 4.2.2

octlibdir = /opt/local/lib/octave/${version}/
octincdir = /opt/local/include/octave-${version}/octave

# http://stackoverflow.com/questions/7806418/using-setenv-in-makefile
export XTRA_CFLAGS=-std=c99 -UCountAlloc -DMmex -DUse_simd -DUse_thread -O3
