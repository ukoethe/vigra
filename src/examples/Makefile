# uncomment if CXX does not expand to the compiler you want to use
# CXX = g++

# update path according to where 'vigra-config' is installed
CXXFLAGS  := -O $(shell ../../config/vigra-config --cppflags)
IMPEX_LIB := $(shell ../../config/vigra-config --impex-lib)
# if your 'make'-program does not support $(shell ...), write
#    CXXFLAGS  = -O `../../config/vigra-config --cxxflags`
#    IMPEX_LIB = `../../config/vigra-config --impex-lib`

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $*.cxx

TARGETS = \
    convert \
    subimage \
    invert \
    invert_explicitly \
    resize \
    smooth \
    profile \
    pyramid \
    edge \
    watershed \
    voronoi

all: $(TARGETS)

convert:  convert.o
	$(CXX) -o convert convert.o $(IMPEX_LIB)

subimage:  subimage.o
	$(CXX) -o subimage subimage.o $(IMPEX_LIB)

invert:  invert.o
	$(CXX) -o invert invert.o $(IMPEX_LIB)

invert_explicitly:  invert_explicitly.o
	$(CXX) -o invert_explicitly invert_explicitly.o $(IMPEX_LIB)

resize:  resize.o
	$(CXX) -o resize resize.o $(IMPEX_LIB)

smooth:  smooth.o
	$(CXX) -o smooth smooth.o $(IMPEX_LIB)

profile:  profile.o
	$(CXX) -o profile profile.o $(IMPEX_LIB)

pyramid:  pyramid.o
	$(CXX) -o pyramid pyramid.o $(IMPEX_LIB)

edge:  edge.o
	$(CXX) -o edge edge.o $(IMPEX_LIB)

watershed:  watershed.o
	$(CXX) -o watershed watershed.o $(IMPEX_LIB)

voronoi:  voronoi.o
	$(CXX) -o voronoi voronoi.o $(IMPEX_LIB)

clean:
	rm -f *.o $(TARGETS) *.exe core *.plg *ilk *.lib *.pdb *.exp \
           res.gif voronoi.gif distances.gif profile.gif
	rm -rf Release Debug ii_files

