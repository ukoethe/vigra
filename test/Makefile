TOPTESTDIR=.

QUICKSUBDIRS =  \
        syntax \
        error \
        utilities \
        pixeltypes \
        colorspaces \
        imgproc \
        functorexpression \
        convolution \
        morphology \
        fourier_fftw3 \
        math \
        tensorimaging \
        simpleanalysis


SUBDIRS = $(QUICKSUBDIRS) \
        image \
        multiarray \
        multiconvolution \
        impex 

all:
	@rm -f testSuccess $(TOPTESTDIR)/testSuccess ; \
         for i in $(SUBDIRS); do \
	        cd $$i; \
	        make || exit 1; \
	        cd ..; \
	 done ; \
	 touch testSuccess $(TOPTESTDIR)/testSuccess

quick:
	@for i in $(QUICKSUBDIRS); do \
	        cd $$i; \
	        make || exit 1; \
	        cd ..; \
	 done

builddir:
	@if [ "x$(BUILDDIR)" = "x" ]; then \
            echo "Please call 'make builddir BUILDDIR=<builddir>'" ; \
            exit 1; \
         fi
	@if [ -d "$(BUILDDIR)" ]; then \
            echo "$(BUILDDIR) already exists" ; \
            exit 0; \
         fi
	@if [ ! -d "../current/$(BUILDDIR)" ]; then \
            echo "There must be a corresponding build dir in ../current" ; \
            exit 1; \
         fi
	@mkdir $(BUILDDIR)
	@testdir=`pwd`; \
         cat Makefile | sed "s@TOPTESTDIR=.@TOPTESTDIR=$$testdir@" > $(BUILDDIR)/Makefile
	@echo "BUILDDIR=$(BUILDDIR)" > $(BUILDDIR)/Makefile.include
	@testdir=`pwd`; \
         cat Makefile.include | sed "s@TOPTESTDIR=..@TOPTESTDIR=$$testdir@" >> $(BUILDDIR)/Makefile.include
	@for i in $(SUBDIRS); do \
	        mkdir $(BUILDDIR)/$$i; \
	        if [ "x$$i" != "xsyntax" ]; then \
                    cat $$i/Makefile | sed "s@curdir=.@curdir=$$\(TOPTESTDIR\)/$$i@" > $(BUILDDIR)/$$i/Makefile; \
                else \
                    cp $$i/Makefile $(BUILDDIR)/$$i/Makefile; \
                fi; \
	 done
        

vc:
	@rm -f testSuccess ; \
	 for i in $(SUBDIRS); do \
	        cd $$i; \
	        make $$i.exe || exit 1; \
	        cd ..; \
	 done ; \
	 touch testSuccess

clean:
	rm -f *.opt testSuccess	
	@for i in $(SUBDIRS); do \
	        cd $$i; \
	        make clean; \
	        cd ..; \
	 done
