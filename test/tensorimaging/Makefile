curdir=.

TESTIMAGES = l2.xv riesz00.xv riesz10.xv riesz01.xv riesz20.xv riesz11.xv riesz02.xv boundaryTensor.xv l2_boundary1.xv l2_boundary.xv l2_hourglass.xv

test: test.o $(TESTIMAGES)
	$(CXX) -o test test.o $(IMPEX_LIB)
	./test || { rm ./test ; exit 1; }

tensorimaging.exe: test.d
	devenv "..\test.sln" /build Release /project tensorimaging
	./tensorimaging.exe || { rm ./tensorimaging.exe ; exit 1; }

$(TESTIMAGES):
	ln -s $(curdir)/$@ .

include ../Makefile.include

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif
