# Makefile
# copy general purpose .h and .c files here for distribution
# anything in this directory is just a copy from elsewhere (except Makefile)

all:
	@echo 'choose "get"'

#
# jf: copy files from master source directories to these local copies
#
home = $(HOME)
jf_def = $(home)/l/src/defs/
jf_thr = $(home)/l/src/util/thread/
jf_umx = $(home)/l/src/util/mex/

l_def =	\
		defs-env.h
l_thr = \
		jf,thread1.c \
		jf,thread1.h \
		jf,time.c \
		jf,time.h
l_umx = \
		jf,mex,def.h \
		mexarg.c \
		def,mexarg.h

get:	$(l_def) $(l_thr) $(l_umx)

$(l_def): %: $(jf_def)/%
	ls -l $< $@
	cp -pi $< $@
	chmod 644 $@

$(l_thr): %: $(jf_thr)/%
	ls -l $< $@
	cp -pi $< $@

$(l_umx): %: $(jf_umx)/%
	ls -l $< $@
	cp -pi $< $@
