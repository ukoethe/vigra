test: test.o
	$(CXX) -o test test.o
	./test || { rm ./test ; exit 1; }

vc: functorexpression.exe
	./functorexpression.exe || { rm ./functorexpression.exe ; exit 1; }

functorexpression.exe:
	devenv "..\test.sln" /build Release /project functorexpression

ifneq "$(MAKECMDGOALS)" "clean"
include test.d
endif

include ../Makefile.include
