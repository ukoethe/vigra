include ../../config/Makefile.include

TARGET = libvigraimpex.la

OBJS = \
    block_read.lo \
    cr_image.lo \
    freeimage.lo \
    fullpath.lo \
    imagesize.lo \
    machorder.lo \
    machtype.lo \
    order.lo \
    readheader.lo \
    readimage.lo \
    utils.lo \
    writeimage.lo \
    impex.lo \
    bmp.lo \
    gif.lo \
    jpeg.lo \
    sun.lo \
    tiff.lo \
    viff.lo \
    pnm.lo \
    colors.lo \
    error.lo \
    image.lo \
    quantize.lo \
    utility.lo 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LINKCXXLIB) -version-info $(MAJOR_VERSION):$(MINOR_VERSION)$(MICRO_VERSION):0 -o $(TARGET) $(OBJS)

install: all
	$(INSTALL) -d $(libdir)
	$(LIBINSTALL) $(TARGET) $(libdir)

clean:
	rm -f *.o *.lo core a.out $(TARGET)
	rm -f *.plg	
	rm -rf Debug Release ii_files .libs
