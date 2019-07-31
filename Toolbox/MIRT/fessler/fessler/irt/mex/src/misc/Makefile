# mex/src/misc/Makefile

ifndef src_top
	src_top := $(HOME)/l/src
endif

ifndef Makedef
	Makedef := $(src_top)/make/Makedef.$(Arch)
endif

include $(Makedef)

def = ../def
incs=	-I$(def) \
	$(MatlabInc) -DMATLAB_HOST -DMmex

CFLAGS= $(DEBUG) $(BCFLAGS) $(incs)

doth=	$(def)/*.h $(MatlabH)

goal1=	$(bin_mex)/delaysum1_mex.mex$(MatlabSuf)
goal2=	$(bin_mex)/ir_shrink1_mex.mex$(MatlabSuf)
goal3=	$(bin_mex)/ir_tridiag_inv_mex.mex$(MatlabSuf)

#
# trick: to avoid use of the IO libraries, we disable CountAlloc here
# so that the standard calloc is used.
#
com= $(MatlabCom) -UCountAlloc -v

all:	$(goal1) $(goal2) $(goal3)

g1c = $(def)/mexarg.c $(def)/jf,thread1.c delaysum1_mex.c
$(goal1):	$(doth) $(g1c) $(libs)
		$(com) -output $@ $(incs) -DNeed_delaysum1_mex_gateway \
			$(g1c) $(libs)
		chmod 755 $@

g2c = $(def)/mexarg.c $(def)/jf,thread1.c ir_shrink1.c ir_shrink1_mex.c
$(goal2):	$(doth) $(g2c) $(libs) ir_shrink1.h
		$(com) -output $@ $(incs) -DNeed_ir_shrink1_mex_gateway \
			$(g2c) $(libs)
		chmod 755 $@

g3c = $(def)/mexarg.c $(def)/jf,thread1.c ir_tridiag_inv_mex.c
$(goal3):	$(doth) $(g3c) $(libs)
		$(com) -output $@ $(incs) -DNeed_ir_tridiag_inv_mex_gateway \
			$(g3c) $(libs)
		chmod 755 $@

reset:
	rm -i $(goal1) $(goal2) $(goal3)
