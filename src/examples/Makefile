include ../../config/Makefile.include

                 		   
FILES = \
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
    
all: $(FILES)

convert:  convert.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o convert convert.o $(LIBS) 

subimage:  subimage.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o subimage subimage.o $(LIBS) 

invert:  invert.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o invert invert.o $(LIBS) 

invert_explicitly:  invert_explicitly.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o invert_explicitly invert_explicitly.o $(LIBS) 

resize:  resize.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o resize resize.o $(LIBS) 

smooth:  smooth.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o smooth smooth.o $(LIBS) 

profile:  profile.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o profile profile.o $(LIBS) 

pyramid:  pyramid.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o pyramid pyramid.o $(LIBS) 

edge:  edge.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o edge edge.o $(LIBS) 

watershed:  watershed.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o watershed watershed.o $(LIBS) 

voronoi:  voronoi.o 
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o voronoi voronoi.o $(LIBS) 

clean:
	rm -f *.o $(FILES) *.exe core *.plg *ilk *.lib *.pdb *.exp \
           res.gif voronoi.gif distances.gif profile.gif
	rm -rf Release Debug ii_files

