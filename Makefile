#
# Copyright (c) 2022 Roelof Bart Toonen
# License: MIT license (spdx.org MIT)
# 
# This Makefile creates a directory 'NldsPred', the root directory for the R
# package. It will copy the required files to that directory and necessary
# subdirectories and creates the package from there.
#

ORIGINROOTDIR = $(CURDIR)
ORIGINSRCDIR = $(ORIGINROOTDIR)/src
ORIGINLIBDIR = $(ORIGINROOTDIR)/nldspredl

VPATH = $(ORIGINROOTDIR):$(ORIGINSRCDIR):$(ORIGINLIBDIR)/src:$(ORIGINLIBDIR)/include

TARGETROOTDIR = $(CURDIR)/NldsPred
TARGETSRCDIR  = $(TARGETROOTDIR)/src
TARGETLIBDIR  = $(TARGETSRCDIR)/libnlpre
TARGETRDIR    = $(TARGETROOTDIR)/R

RLIBREQ = bundle.c bundle.h bundle_log.h datalimits.h dqrdc.f \
  dsvdc.f dtls.f embed.h heap.c heap.h housh.f kdt.c kdt.h fn.c \
  fn.h fn_log.h fn_exp.c fn_exp.h fn_exp_log.h fn_tls.c \
  fn_tls.h fn_tls_log.h log.c log.h log_meta.h logtotbl.c logtotbl.h \
  mkembed.c mkembed.h mkembed_log.h point.c point.h sets.c sets.h \
  sets_log.h stat.c stat.h stat_log.h tr2.f traverse.c traverse.h \
  traverse_log.h traverse_stat.h tsfile.c tsdat.h tsfile.h tstoembdef.c \
  tstoembdef.h

RROOTREQ = DESCRIPTION NAMESPACE

RSRCREQ = NldsPredGetTable.cpp NldsPredRun.cpp Makevars 

.PHONY: clean build install

build: cpall
	cd $(ORIGINROOTDIR); \
		Rscript --vanilla -e "Rcpp::compileAttributes(pkgdir=\"$(TARGETROOTDIR)\",verbose=TRUE)"; \
		R CMD build NldsPred

install: build
	cd $(ORIGINROOTDIR); \
		sudo R CMD INSTALL NldsPred

clean: cleantarget

# The code below is used to create directories in the targetdir if they do not
# yet exist and to copy the required files (searched for in VPATH directories)
# to the correct target directory.

.PHONY: rootfiles srcfiles libfiles rfiles $(TARGETROOTDIR) $(TARGETSRCDIR) \
	$(TARGETLIBDIR) $(TARGETRDIR) cprootfiles cpsrcfiles cplibfiles \
	cprfiles cpall cleantarget

# cpall: rootfiles srcfiles libfiles rfiles
cpall: rootfiles srcfiles libfiles

rootfiles: $(TARGETROOTDIR) cprootfiles

$(TARGETROOTDIR):
	+@[ -d $@ ] || mkdir -p $@

cprootfiles: $(RROOTREQ)
	cp -f -t $(TARGETROOTDIR)/ $^

srcfiles: $(TARGETSRCDIR) cpsrcfiles

$(TARGETSRCDIR):
	+@[ -d $@ ] || mkdir -p $@

cpsrcfiles: $(RSRCREQ)
	cp -f -t $(TARGETSRCDIR)/ $^

libfiles: $(TARGETLIBDIR) cplibfiles

$(TARGETLIBDIR):
	+@[ -d $@ ] || mkdir -p $@

cplibfiles: $(RLIBREQ)
	cp -f -t $(TARGETLIBDIR)/ $^

rfiles: $(TARGETRDIR) cprfiles

$(TARGETRDIR):
	+@[ -d $@ ] || mkdir -p $@

cprfiles: $(RRREQ)
	cp -f -t $(TARGETRDIR)/ $^

cleantarget:
	rm -rf $(TARGETROOTDIR)
