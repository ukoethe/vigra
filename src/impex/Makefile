include ../../config/Makefile.include

TARGET = libvigraimpex.la

OBJS = \
    gif.lo \
    bmp.lo \
    byteorder.lo \
    codecmanager.lo \
    imageinfo.lo \
    jpeg.lo \
    png.lo \
    pnm.lo \
    sun.lo \
    tiff.lo \
    viff.lo \
    void_vector.lo 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LINKCXXLIB) -version-info 2:$(MAJOR_VERSION)$(MINOR_VERSION)$(MICRO_VERSION):0 -o $(TARGET) $(OBJS)

install: all
	$(INSTALL) -d $(libdir)
	$(LIBINSTALL) $(TARGET) $(libdir)

vcinstall: vigraimpex.lib
	cp vigraimpex.lib /cygdrive/c/Users/Koethe/lib/msvc.net

vigraimpex.lib: $(patsubst %.lo, %.ld, $(OBJS))  # MSVC
	devenv "..\vigra.sln" /build Release /project impex

clean:
	rm -f *.o *.lo *.ld core a.out $(TARGET)
	rm -f *.plg	
	rm -rf Debug Release ii_files .libs

ifeq  "$(DEVELOPMENT_VERSION)" "yes"
ifneq "$(MAKECMDGOALS)" "clean"
  # update and include header dependencies
  include $(patsubst %.lo, %.ld, $(OBJS))
endif
endif
