CC = cc -n32
CXX = CC -n32 -OPT:Olimit=0 -I$$VIGRA_PATH/include/IRIX
CXXLIBFLAGS = -ptused
AR = ar cq
