#EXTRACXXFLAGS =  -DNEWIMPEX
#EXTRACXXINCLUDES =  -I../../experiments/newimpex/include
CXXDEBUG = -g # or -O

#LIBS = -L../../experiments/newimpex/src/impex  -L/software/libpng-1.2.1/SunOS-5.8/lib -lpng $(VIGRAIMPEX_LIB)

LIBS = -L/software/libpng-1.2.1/SunOS-5.8/lib -lpng $(VIGRAIMPEX_LIB)
                 		   
test: test.o ../../current/src/impex/.libs/libvigraimpex.a
	$(CXX) -o test test.o ../../current/src/impex/.libs/libvigraimpex.a $(IMPEX_DEPLIBS)
	./test || { rm ./test ; exit 1; }

vc: impex.exe
	./impex.exe || { rm ./impex.exe ; exit 1; }

impex.exe
	devenv "..\test.sln" /build Release /project impex

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
