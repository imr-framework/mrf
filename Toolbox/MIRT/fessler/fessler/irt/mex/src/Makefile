# mex/src/Makefile
# master Makefile
#

all:
	cd misc;	$(MAKE)

ifeq ($(Arch), mi64)
	@echo 'now compiling 32-bit version'
	$(MAKE) Arch=i386
endif
