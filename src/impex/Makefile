include $(VIGRA_PATH)/config/Makefile.include

CCFLAGS = -O $(HAS_JPEG) $(HAS_TIFF) $(CCLIBFLAGS)

CXXFLAGS = -O $(HAS_JPEG) $(HAS_TIFF) $(CXXLIBFLAGS)
                             
# comment out if libjpeg is not installed
# otherwise update the path to jpeg header files		   
HAS_JPEG = -DHasJPEG -I../../../jpeg

# comment out if libtiff is not installed
# otherwise update the path to tiff header files		   
HAS_TIFF = -DHasTIFF -I../../../tiff

#SHELL = tcsh -f
TARGET = libvigraimpex.a

OBJS = \
    block_read.o \
    cr_image.o \
    freeimage.o \
    fullpath.o \
    imagesize.o \
    machorder.o \
    machtype.o \
    order.o \
    readheader.o \
    readimage.o \
    utils.o \
    writeimage.o \
    impex.o \
    bmp.o \
    gif.o \
    jpeg.o \
    sun.o \
    tiff.o \
    viff.o \
    pnm.o \
    colors.o \
    error.o \
    image.o \
    quantize.o \
    utility.o 


install: all
	@if [ ! -d $(LIBDIR) ]; then mkdir $(LIBDIR) ; fi
	cp $(TARGET) $(LIBDIR)

all: $(TARGET)

$(TARGET): $(OBJS)
	rm -f $(TARGET)
	ar cq $(TARGET) $(OBJS)
	touch $(TARGET)

clean:
	rm -f *.o core a.out $(TARGET)
	rm -f *.plg	
	rm -rf Debug Release ii_files

	
