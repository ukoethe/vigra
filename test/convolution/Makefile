curdir=.

TESTIMAGES = lenna128.xv lenna_simple_sharpening_orig.xv lenna_gaussian_sharpening_orig.xv lenna128sepgrad.tif lennahessxx.tif lennastxx.tif lenna128recgrad.xv lenna128nonlinear.tif resampling.xv lennahessyy.tif lennastyy.tif lennahessxy.tif lennastxy.tif


test: test.o $(TESTIMAGES)
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

convolution.exe: test.d
	devenv "..\test.sln" /build Release /project convolution
	./convolution.exe || { rm ./convolution.exe ; exit 1; }
        
$(TESTIMAGES):
	ln -s $(curdir)/$@ .

include ../Makefile.include

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
