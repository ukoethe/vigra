include ../../config/Makefile.include

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
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o convert convert.o $(IMPEX_LIB)

subimage:  subimage.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o subimage subimage.o $(IMPEX_LIB)

invert:  invert.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o invert invert.o $(IMPEX_LIB)

invert_explicitly:  invert_explicitly.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o invert_explicitly invert_explicitly.o $(IMPEX_LIB)

resize:  resize.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o resize resize.o $(IMPEX_LIB)

smooth:  smooth.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o smooth smooth.o $(IMPEX_LIB)

profile:  profile.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o profile profile.o $(IMPEX_LIB)

pyramid:  pyramid.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o pyramid pyramid.o $(IMPEX_LIB)

edge:  edge.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o edge edge.o $(IMPEX_LIB)

watershed:  watershed.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o watershed watershed.o $(IMPEX_LIB)

voronoi:  voronoi.o
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o voronoi voronoi.o $(IMPEX_LIB)

clean:
	rm -f *.o $(TARGETS) *.exe core *.plg *ilk *.lib *.pdb *.exp \
           res.gif voronoi.gif distances.gif profile.gif
	rm -rf Release Debug ii_files

