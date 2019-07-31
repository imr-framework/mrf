#
# src/util/mex/Makefile
#

include ../List-util
include ./List,util,mex
#include $(IO_SRC)/List-io
#include $(misc_src)/List-misc
misc_dir = ../misc

util_mex_src=$(src_util)/mex
util_mex_lib=$(LIB0_DIR)/lib-util-mex.a

incs=	$(util_inc) \
	-I$(misc_dir) \
	$(MatlabInc) -DMmex -DMATLAB_HOST
CFLAGS= $(DEBUG) $(BCFLAGS) $(incs)
libs=	$(util_mex_lib) # $(misc_lib) $(IO_LIBS)

doth=	$(env_h) $(MatlabH) $(util_mex_h) def,mexarg.h
#doth=	$(env_h) $(MatlabH) $(IO_H) $(util_mex_h) def,mexarg.h
goal1=	$(bin_mex)/rotmex.mex$(MatlabSuf)
goal2=	$(bin_mex)/dtft_mex.mex$(MatlabSuf)
goal3=	$(bin_mex)/penalty_mex$(DEBUG).mex$(MatlabSuf)
goal4=	$(bin_mex)/jf_mex$(DEBUG).mex$(MatlabSuf)

#
# trick: to avoid use of the IO libraries, we disable CountAlloc here
# so that the standard calloc is used.
#
com= $(MatlabCom) -UCountAlloc -v

all:	$(goal1) $(goal2) $(goal3) $(goal4) $(util_mex_lib)

$(util_mex_lib):	$(util_mex_o)
		$(DO_LIB)

$(util_mex_o):	$(doth)

rh = $(misc_dir)/rotate.h $(misc_dir)/defs-misc.h
rc = rotmex.c $(misc_dir)/rotate.c $(misc_dir)/rotate0.c
$(goal1):	$(doth) $(rc) $(libs)
		$(com) -output $@ $(incs) $(rc) $(libs)
		chmod 755 $@

ph = $(misc_dir)/defs-misc.h
pc = penalty_mex.c $(misc_dir)/penalty,diff.c # penalty,2d.c
$(goal3):	$(doth) $(ph) $(pc) $(libs)
		$(com) -output $@ $(incs) $(pc) $(libs)
		chmod 755 $@

dtft_h = def,dtft.h
dtft_c = dtft,mex.c dtft,adj.c dtft,for.c
$(goal2):	$(doth) $(dtft_h) $(dtft_c) $(libs)
		$(com) -output $@ $(incs) -DNeed_dtft_mex_gateway $(dtft_c) $(libs)
		chmod 755 $@

ptab_h = ptab,def.h
ptab_c = ptab,par2,proj.c ptab,par2,back.c ptab,mex.c

fdk_h = def,fdk.h
fdk_c = fdk,st.c fdk,ts.c fdk,mex.c

jf_h = $(ptab_h) $(dtft_h) $(fdk_h)
jf_c = $(ptab_c) $(dtft_c) $(fdk_c) jf_mex.c
$(goal4):	$(doth) $(jf_h) $(jf_c) $(libs)
		$(com) -output $@ $(incs) $(jf_c) $(libs)
		chmod 755 $@

reset:
	rm -i $(goal1) $(goal2) $(goal3) $(goal4) $(util_mex_lib)
