include ../config/Makefile.include

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
        fourier \
        math \
        tensorimaging \
        simpleanalysis


SUBDIRS = $(QUICKSUBDIRS) \
        image \
        multiarray \
        multiconvolution \
        impex 

all:
	@rm -f testSuccess $(vigra_testdir)/testSuccess ; \
         for i in $(SUBDIRS); do \
	        cd $$i; \
	        make || exit 1; \
	        cd ..; \
	 done ; \
	 touch testSuccess $(vigra_testdir)/testSuccess

quick:
	@for i in $(QUICKSUBDIRS); do \
	        cd $$i; \
	        make || exit 1; \
	        cd ..; \
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
