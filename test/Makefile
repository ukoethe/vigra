SUBDIRS = \
        error \
        image \
        impex \
        imgproc \
        functorexpression \
        convolution \
        morphology \
        simpleanalysis

all::
	@rm -f testSuccess ; \
	 for i in $(SUBDIRS); do \
	        cd $$i; \
	        make || exit 1; \
	        cd ..; \
	 done ; \
	 touch testSuccess;

clean::
	rm -f *.ncb *.opt testSuccess	
	@for i in $(SUBDIRS); do \
	        cd $$i; \
	        make clean; \
	        cd ..; \
	done
      
