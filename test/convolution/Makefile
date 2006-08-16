curdir=.

TESTIMAGES = lenna128.xv lenna_simple_sharpening_orig.xv lenna_gaussian_sharpening_orig.xv lenna128sepgrad.xv lennahessxx.xv lennastxx.xv lenna128recgrad.xv lenna128nonlinear.xv resampling.xv lennahessyy.xv lennastyy.xv lennahessxy.xv lennastxy.xv


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
