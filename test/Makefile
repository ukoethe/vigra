QUICKSUBDIRS =  \
        error \
        pixeltypes \
        colorspaces \
        imgproc \
        functorexpression \
        convolution \
        morphology \
        fourier \
        simpleanalysis


SUBDIRS = $(QUICKSUBDIRS) \
        image \
        multiarray \
        impex 

all:
	@rm -f testSuccess ; \
         for i in $(SUBDIRS); do \
	        cd $$i; \
	        make || exit 1; \
	        cd ..; \
	 done ; \
	 touch testSuccess

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
