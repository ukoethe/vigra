curdir=.

test: test.o
	$(CXX) -o test test.o ||  { rm ./test.cxx ; exit 1; }
	./test || { rm ./test ; exit 1; }

error.exe: test.d
	devenv "..\test.sln" /build Release /project error
	./error.exe || { rm ./error.exe ; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include

CPPFLAGS := $(LIBCPPFLAGS)

HEADERS=$(shell ls ../../current/include/vigra/*.hxx)

test.cxx: $(HEADERS)
	@echo creating test.cxx
	@\rm -f test.cxx && touch test.cxx
	@for i in $(HEADERS); do \
            j=`expr $$i : ".*/\(.*\)"` ; \
            if [ "$$j" != "fftw.hxx" -a "$$j" != "rfftw.hxx" ] ; then \
                echo "#include <vigra/$$j>" >> test.cxx; \
            fi; \
        done
	@echo "int main() {}" >> test.cxx

            
            
