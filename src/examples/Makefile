include $(VIGRA_PATH)/config/Makefile.include

CXXFLAGS = -g -O -I../include 
LIB = -lvigraimpex -ljpeg -ltiff -lm
                 		   
FILES = \
    convert \
    subimage \
    invert \
    invert_explicitly \
    smooth \
    profile \
    pyramid \
    edge \
    watershed \
    voronoi
    
all: $(FILES)

convert:  convert.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o convert convert.o $(CXXLIBDIR) $(LIB)

subimage:  subimage.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o subimage subimage.o $(CXXLIBDIR) $(LIB)

invert:  invert.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o invert invert.o $(CXXLIBDIR) $(LIB)

invert_explicitly:  invert_explicitly.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o invert_explicitly invert_explicitly.o $(CXXLIBDIR) $(LIB)

smooth:  smooth.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o smooth smooth.o $(CXXLIBDIR) $(LIB)

profile:  profile.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o profile profile.o $(CXXLIBDIR) $(LIB)

pyramid:  pyramid.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o pyramid pyramid.o $(CXXLIBDIR) $(LIB)

edge:  edge.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o edge edge.o $(CXXLIBDIR) $(LIB)

watershed:  watershed.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o watershed watershed.o $(CXXLIBDIR) $(LIB)

voronoi:  voronoi.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o voronoi voronoi.o $(CXXLIBDIR) $(LIB)

clean:
	rm -f *.o $(FILES) *.exe core *.plg *ilk *.lib *.pdb *.exp \
           res.gif voronoi.gif distances.gif profile.gif
	rm -rf Release Debug ii_files

