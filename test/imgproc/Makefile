curdir=.

TESTIMAGES = lenna128.xv lenna128rgb.xv splineimageview2.xv splineimageview3.xv splineimageview5.xv lenna42lin.xv lenna288neu.xv lenna42neu.xv lenna288rgbneu.xv lenna42rgbneu.xv lenna367FIR.xv lenna42FIR.xv lenna367IIR.xv lenna42IIR.xv lenna42linrgb.xv lennargb42FIR.xv lennargb42IIR.xv

test: test.o $(TESTIMAGES)
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

imgproc.exe: test.d
	devenv "..\test.sln" /build Release /project imgproc
	./imgproc.exe || { rm ./imgproc.exe ; exit 1; }

$(TESTIMAGES):
	ln -s $(curdir)/$@ .

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
