CXXFLAGS = -O -I../include `vigra-config --cxxflags` -I/software/boost-1.25.0/include

#OBJ = one_band_basic_image_test_suite.o 
#OBJ = gray_image_test_suite.o 
#OBJ = select_band_image_test_suite.o
#OBJ = variable_bands_image_test_suite.o
#OBJ = rgb_basic_image_test_suite.o
#OBJ = rgb_imagehierarchy_test_suite.o
#OBJ = single_band_image_test_suite.o
#OBJ = vector_2_basic_image_test_suite.o
#OBJ = vector_2_imagehierarchy_test_suite.o
#OBJ = vector_3_basic_image_test_suite.o
#OBJ = vector_3_imagehierarchy_test_suite.o
#OBJ = vector_4_basic_image_test_suite.o
#OBJ = vector_4_imagehierarchy_test_suite.o

OBJ = tests_start_suite.o \
one_band_basic_image_test_suite.o \
gray_image_test_suite.o \
vector_2_basic_image_test_suite.o \
vector_2_imagehierarchy_test_suite.o \
vector_3_basic_image_test_suite.o \
vector_3_imagehierarchy_test_suite.o \
vector_4_basic_image_test_suite.o \
vector_4_imagehierarchy_test_suite.o \
variable_bands_image_test_suite.o \
rgb_imagehierarchy_test_suite.o \
rgb_basic_image_test_suite.o \
select_band_image_test_suite.o \
single_band_image_test_suite.o

test: $(OBJ)
	$(CXX) $(CXXFLAGS) -o test $(OBJ) `vigra-config --libs`
	./test || { rm ./test ; exit 1; }

# Fuer Tests von IRGBImage, FRGBImage, DRGBImage, BRGBImage den Kommentar loeschen.
#MAKROS = -DRGB_IMAGE 

# Fuer Tests von FImage, DImage, SImage, BImage, IImage den Kommentar loeschen.
#MAKROS = -DnRGB_IMAGE 

# Fuer Tests von FVector2Image, FVector3Image, FVector4Image, DVector2Image, DVector3Image, DVector4Image den Kommentar loeschen.
#MAKROS = -DVImage 

# Fuer Tests von VariableBandsImage den Kommentar loeschen.
#MAKROS = -DVARIABLE_BANDS_IMAGE -DIMAGEHIERARCHY

# Fuer Tests von SelectBandImage den Kommentar loeschen.
#MAKROS = -DSELECT_BAND_IMAGE -DIMAGEHIERARCHY

# Fuer Tests von GrayImage den Kommentar loeschen.
#MAKROS = -DGRAY_IMAGE -DIMAGEHIERARCHY

# Fuer Tests von SingleBandImage den Kommentar loeschen.
#MAKROS = -DSINGLE_BAND_IMAGE -DIMAGEHIERARCHY

# Fuer Tests von FixedBandsImage den Kommentar loeschen.
#MAKROS = -DFIXED_V_IMAGE -DIMAGEHIERARCHY

# Fuer Tests von FixedRGBImage den Kommentar loeschen.
#MAKROS = -DFIXED_RGB_IMAGE -DIMAGEHIERARCHY

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $*.cxx

%.d: %.cxx
	@echo updating header dependencies for $*.cxx
	@rm -f $*.d
	@(echo "$*.d \c" ; $(CXX) -M $(CXXFLAGS) -c $*.cxx) > $*.d \
        || { rm -f $*.d; exit 1; }

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.o, %.d, $(OBJ))
endif
