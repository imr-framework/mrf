# Makefile
# for compiling fdk_mex

# common compile flags
aflags = -std=c99 \
	-ffast-math \
	-O3 -fexpensive-optimizations \
        -Wall -Wshadow -Wpointer-arith -W -Werror -Wmissing-prototypes \
        -DUse_thread \
	-DUse_fdk_mex

# OS specific compile flags
ifeq ($(Arch), mi64) # mac 64-bit
	oflags = -arch x86_64 -DUse_ncore_sysctl # -DUse_aff_mac1
	libs =
	mextop = ~fessler/l/matlab/extern
	mexinc = -I$(mextop)/include
	mexsuf = maci64
	mexgnu = -DMATLAB_HOST -DMmex
	mexcom = $(mextop)/../bin/mex CFLAGS='$(aflags) $(oflags) $(mexgnu)'

else # linux 64 bit
	oflags = -DNeed_uint
	libs = -lpthread -lm
	# user may need to edit the following line to set up path to matlab:
	mextop = /usr/local/matlab/extern
	mexinc = -I$(mextop)/include
	mexsuf = a64
	mexgnu = -D_GNU_SOURCE -UNeed_uint -fPIC -DMATLAB_HOST -DMmex
	mexcom = $(mextop)/../bin/mex CFLAGS='$(aflags) $(oflags) $(mexgnu)'

endif

CFLAGS = $(aflags) $(oflags)

goal = fdk_mex.mex$(mexsuf)


###all: get fdk
#	@echo 'choose "get" (jf only) or "fdk" to make $(goal)'

#
# jf: copy files from master source directories to these local copies
#
###home = $(HOME)
###jf_cbc = $(home)/l/src/util/cbct/
###jf_umx = $(home)/l/src/util/mex/
###
###l_cbc = \
	cbct,def.h \
	cbct,mask2.c
###
###l_umx = \
	def,fdk.h \
	fdk,mex.c \
	fdk,st.c \
	fdk,ts.c \
	fdk,ts,t.c
###
###get: $(l_cbc) $(l_umx)
###
###$(l_cbc): %: $(jf_cbc)/%
###	ls -l $< $@
###	cp -pi $< $@
###
###$(l_umx): %: $(jf_umx)/%
###	ls -l $< $@
###	cp -pi $< $@


# compile it all together

def = def
incs = -I. -I$(def)

allc = ./*.c \
		$(def)/mexarg.c \
		$(def)/jf,thread1.c \
		$(def)/jf,time.c
allh = ./*.h \
		$(def)/defs-env.h \
		$(def)/def,mexarg.h \
		$(def)/jf,mex,def.h \
		$(def)/jf,thread1.h \
		$(def)/jf,time.h

fdk:	$(goal)

$(goal):	$(allh) $(allc)
	$(mexcom) -o $@ $(allc) $(incs) $(libs)

reset:
	rm -i $(goal)
