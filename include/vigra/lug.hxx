/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 
 
#ifndef IMP_LUG_HXX
#define IMP_LUG_HXX

// lug includes 

// vigra includes 
#include "vigra/utilities.hxx"
#include "vigra/error.hxx"

#include <lug.h>

typedef bitmap_hdr LugImage;

extern "C"
{
    int allocatebitmap(bitmap_hdr*,int,int,int,int);
    int isagrayscaled(bitmap_hdr*);
    int write_lug_file(char*,bitmap_hdr*);
    int freebitmap(bitmap_hdr*);
    int read_lug_file(char*,bitmap_hdr*);
    int write_rgb_file(char*,bitmap_hdr*);
    int read_rgb_file(char*,bitmap_hdr*,int,int);
}

struct LugRGBEntry
{
        unsigned char red;
        unsigned char green;
        unsigned char blue;
};

inline 
byte *fill_bw_pallete(byte * buffer)
{
  int i;
  byte *ptr;

  /* I'll use a pointer */
  ptr = buffer;
  for (i = 0; i < 256; i++) {
    *ptr++= (byte) i;   /* R */
    *ptr++= (byte) i;   /* G */
    *ptr++= (byte) i;   /* B */
  }

  return buffer;
}

/** \page LUGImpex Import/export via the LUG library

    Supports GIF, TIFF, JPEG, PS, and other formats.

    <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
              \ref LUGFunctions "LUG Functions"
        <DD> Read/write/delete LUG images

        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
            \ref ConvertLugImages
        <DD> <em>Convert LUG to/from any images supporting ImageIterators</em>
    </DL>

    LUG is a simple gateway for many common image
    data formats, including GIF, TIFF, JPEG and PostScript. The source code for LUG
    (Libraria de Utilidades Graficas - Graphics Utilities Library) is
    available at 
    <a href="http://www.uniovi.es/~rivero/LUG/">http://www.uniovi.es/~rivero/LUG/</a>.
    LUG supports some formats (e.g. GIF, PostScript) on its own and uses 
    other public domain libraries for the rest (TIFF, JPEG). Unlike VIFF, 
    LUG can only handle 'byte' (scalar, RGB and mapped RGB) pixel types. 
*/

/********************************************************/
/*                                                      */
/*             read/write/delete LugImages              */
/*                                                      */
/********************************************************/

/** \defgroup LUGFunctions Read/write/delete LUG Images
    These functions are convenience functions to simplify LUG's usage 
    in the context of VIGRA. You can
    always use the original LUG functions for the same or related purposes 
    (see LUG documentation for details). 
*/
//@{
    /** \brief Create LugImage and read contents from given file.
        
        <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
        
    */
inline
LugImage *
readLugImage(const char* filename)
{
        LugImage * image = new LugImage;
        if (!image)
                return 0;

        read_lug_file(const_cast<char*>(filename), image);
        return image;
}

    /** \brief Write LugImage to given file. 
    
        The output file format is determined 
        by the file name extension. For example, 'test.gif' will be a GIF file.
        You can override this convention by using alternative write functions.
        See the LUG manual for details. 
    
        <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
        
    */
inline
void
writeLugImage(const char* filename, LugImage* image)
{
        if (!image)
                return;

        write_lug_file(const_cast<char*>(filename), image);
}

    /** \brief Create LugImage if given size and color depth.
    
        Depth is either 8 or 24.
    
        <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
        
    */
inline
LugImage *
createLugImage(int width, int height, int depth)
{
        LugImage * img = new LugImage;
        allocatebitmap(img, width, height, depth, 0);
        return img;
}

inline
LugImage *
createLugImage(int width, int height, int depth, int colors)
{
        LugImage * img = new LugImage;
        allocatebitmap(img, width, height, depth, colors);
        return img;
}

inline
LugImage *
createLugImage()
{
        LugImage * img = new LugImage;
        return img;
}

    /** \brief Delete LugImage and free memory.
    
        <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
        
    */
inline
void
freeLugImage(LugImage * image)
{
        freebitmap(image);
        delete image;
}

    /** \brief Test if LugImage is a grayscale image.
    
        <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
        
    */
int
isagrayscaled(LugImage *);

//@}

/** \defgroup ConvertLugImages Convert LUG Images
    to/from any images supporting ImageIterators
*/
//@{
/********************************************************/
/*                                                      */
/*                     importLugImage                   */
/*                                                      */
/********************************************************/

/** \brief Convert given lug into image specified by iterator range.

    Accessors are used to write the data.    
    This function calls \ref lugToScalarImage() or \ref lugToRGBImage(), depending on 
    the accessor's value_type.

    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    template <class ImageIterator, class Accessor>
    void
    importLugImage(LugImage * lug, ImageIterator iter, Accessor a)
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    template <class ImageIterator, class Accessor>
    void
    importLugImage(LugImage * lug, pair<ImageIterator, Accessor> dest)
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
    
    \code
    LugImage * lug = readimage("scalarimage.gif");
    
    BImage img(lug->xsize, lug->ysize);
    
    importLugImage(lug, destImage(img));
    
    freeLugImage(lug);
    \endcode
    
    <b> Required Interface:</b>
    
    see \ref lugToScalarImage() and \ref lugToRGBImage()
    
    <b> Preconditions:</b>
    
    see \ref lugToScalarImage() and \ref lugToRGBImage()
    
*/
template <class ImageIterator, class Accessor>
inline void
importLugImage(LugImage * lug, ImageIterator iter, Accessor a)
{
    typedef typename 
        NumericTraits<typename Accessor::value_type>::isScalar
        isScalar;
    importLugImage(lug, iter, a, isScalar());
}

template <class ImageIterator, class Accessor>
inline void
importLugImage(LugImage * lug, pair<ImageIterator, Accessor> dest)
{
    importLugImage(lug, dest.first, dest.second);
}

template <class ImageIterator, class Accessor>
inline void
importLugImage(LugImage * lug, ImageIterator iter, Accessor a, VigraTrueType)
{
    lugToScalarImage(lug, iter, a);
}

template <class ImageIterator, class Accessor>
inline void
importLugImage(LugImage * lug, ImageIterator iter, Accessor a, VigraFalseType)
{
    lugToRGBImage(lug, iter, a);
}

/********************************************************/
/*                                                      */
/*                   lugToScalarImage                   */
/*                                                      */
/********************************************************/

/** \brief Convert single-band LugImage to scalar image.

    This function uses accessors to write the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    template <class ImageIterator, class Accessor>
    void
    lugToScalarImage(LugImage * img, ImageIterator iter, Accessor a)
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    template <class ImageIterator, class Accessor>
    void
    lugToScalarImage(LugImage * img, pair<ImageIterator, Accessor> dest)
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
    
    \code
    LugImage * lug = readLugImage("scalarimage.tif");
    
    BImage img(lug->xsize, lug->ysize);
    
    lugToScalarImage(lug, destImage(img));
    
    freeLugImage(lug);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    unsigned char value;
    
    Accessor accessor;
                       
    accessor.set(value, upperleft);
    \endcode
    
    <b> Preconditions:</b>
    
    \code
    isagrayscaled(lug)
    \endcode
    
*/
template <class ImageIterator, class Accessor>
void
lugToScalarImage(LugImage * img, ImageIterator iter, Accessor a)
{
    int w = img->xsize;
    int h = img->ysize;

    vigra_precondition(img,
                "lugToScalarImage(): " 
                      "NULL pointer to input data.");
   
    vigra_precondition(isagrayscaled(img),
                "lugToScalarImage(): " 
                      "Source image is not scalar.");

    // now we have a color map with 256 entries and r=g=b
    ImageIterator yd(iter);
        
    unsigned char * pindex = img->r;
    struct LugRGBEntry * cmap = (struct LugRGBEntry *) img->cmap;
    for(int y = 0; y < h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(int x = 0; x < w; ++x, ++pindex, ++xd.x)
        {
            a.set(cmap[*pindex].red, xd);
        }
    }
}

template <class ImageIterator, class Accessor>
inline
void
lugToScalarImage(LugImage * img, pair<ImageIterator, Accessor> dest)
{
    lugToScalarImage(img, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                     lugToRGBImage                    */
/*                                                      */
/********************************************************/

/** \brief Convert LugImage to RGB image.

    This function uses accessors to write the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    template <class ImageIterator, class RGBAccessor>
    void
    lugToRGBImage(LugImage * img, ImageIterator iter, RGBAccessor a)
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    template <class ImageIterator, class RGBAccessor>
    void
    lugToRGBImage(LugImage * img, pair<ImageIterator, RGBAccessor> dest)
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
    
    \code
    LugImage * lug = readLugImage("rgbimage.gif");
    
    BRGBImage img(lug->xsize, lug->ysize);
    
    lugToRGBImage(lug, destImage(img));
    
    freeLugImage(lug);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    unsigned char value;
    
    RGBAccessor accessor;
                       
    accessor.setRed(value, upperleft);
    accessor.setGreen(value, upperleft);
    accessor.setBlue(value, upperleft);
    \endcode
    
*/
template <class ImageIterator, class Accessor>
void
lugToRGBImage(LugImage * img, ImageIterator upperleft, Accessor a)
{
    vigra_precondition(img,
            "lugToRGBImage(LugImage *, RGBImageIterator): " 
            "NULL pointer to input data.");

    int x = 0;
    int y = 0;
    int w = img->xsize;
    int h = img->ysize;
        
    if(img->depth > 8)
    {
        // real multiband image
        vigra_precondition( (img->r && img->g && img->b),
                "lugToRGBImage(): " 
                "NULL pointer to pixel data.");

        ImageIterator yd(upperleft);

        unsigned char * pred = img->r;
        unsigned char * pgreen = img->g;
        unsigned char * pblue = img->b;
        for(y = 0; y < h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(x = 0; x < w; ++x, ++pred, ++pgreen, ++pblue, ++xd.x)
            {
                    a.setRed(*pred, xd);
                    a.setGreen(*pgreen, xd);
                    a.setBlue(*pblue, xd);
            }
        }
    }
    else
    {
        // scalar data with color map
        vigra_precondition(img->r,
                "lugToRGBImage(): " 
                "NULL pointer to pixel data.");
        // data only in r-buffer
        unsigned char * ps = img->r;
        
        struct LugRGBEntry * cmap = (struct LugRGBEntry *) img->cmap;
    
        ImageIterator yd(upperleft);
    
        for(y = 0; y < h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(x = 0; x < w; ++x, ++ps, ++xd.x)
            {
                    a.setRed(cmap[*ps].red, xd);
                    a.setGreen(cmap[*ps].green, xd);
                    a.setBlue(cmap[*ps].blue, xd);
            }
        }
    }
}

template <class ImageIterator, class Accessor>
inline
void
lugToRGBImage(LugImage * img, pair<ImageIterator, Accessor> dest)
{
    lugToRGBImage(img, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                     createLugImage                   */
/*                                                      */
/********************************************************/

/** \brief Create a LUG image from the given iterator range.

    It is automatically determined whether a scalar or RGB image must 
    be created. This function uses accessors to read the data. 
    Note, however, that LUG images can only store 'unsigned char' pixel
    values, so all scalar types are converted to 'unsigned char' during conversion, 
    while all RGB types are converted to RGBValue<unsigned char>. Use the 
    <a href="VIFFrelatedfunctions.html">VIFF</a> image format, if this is not 
    acceptable.
    
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    template <class ImageIterator, class Accessor>
    inline LugImage *
    createLugImage(ImageIterator upperleft, ImageIterator lowerright, 
                    Accessor a)
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    template <class ImageIterator, class Accessor>
    inline LugImage *
    createLugImage(triple<ImageIterator, ImageIterator, Accessor> src)
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"
    
    \code
    BImage img(width, height);
    
    ...
    
    LugImage * lug = createLugImage(srcImageRange(img));
    
    // output file format GIF is determined from the file name extension '.gif'
    writeLugImage("output.gif", lug); 
    
    freeimage(lug);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    Accessor accessor;
                           
    accessor(upperleft);   // result written into xvimage
    \endcode
    
*/
template <class ImageIterator, class Accessor>
inline LugImage *
createLugImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
{
    typedef typename 
           NumericTraits<typename Accessor::value_type>::isScalar 
           isScalar;
    return createLugImage(upperleft, lowerright, a, isScalar());
}

template <class ImageIterator, class Accessor>
inline LugImage *
createLugImage(triple<ImageIterator, ImageIterator, Accessor> src)
{
    return createLugImage(src.first, src.second, src.third);
}

template <class ImageIterator, class Accessor>
inline LugImage *
createLugImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, VigraFalseType)
{
    return createRGBLugImage(upperleft, lowerright, a);
}

template <class ImageIterator, class Accessor>
inline LugImage *
createLugImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, VigraTrueType)
{
    return createScalarLugImage(upperleft, lowerright, a);
}

/********************************************************/
/*                                                      */
/*                createScalarLugImage                  */
/*                                                      */
/********************************************************/

/** \brief Create a gray-scaled LugImage from the given scalar image.

    The size of the LugImage is determined by the input image.
    This function uses accessors to read the data. All pixel data types are 
    converted to 'unsigned char' during conversion.  Use the 
    <a href="VIFFrelatedfunctions.html">VIFF</a> image format, if this is not 
    acceptable.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    template <class ImageIterator, class Accessor>
    LugImage *
    createScalarLugImage(ImageIterator upperleft, ImageIterator lowerright, Accessor a)
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    template <class ImageIterator, class Accessor>
    LugImage *
    createScalarLugImage(triple<ImageIterator, ImageIterator, Accessor> src)
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
    
    \code
    BImage img(width, height);
    ...
    
    LugImage * lug = createScalarLugImage(srcImageRange(img));
    
    writeLugImage("scalarimage.tif", lug);
    
    freeLugImage(lug);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    unsigned char value;
    
    Accessor accessor;
                       
    value = accessor(upperleft);
    \endcode
    
*/
template <class ImageIterator, class Accessor>
LugImage *
createScalarLugImage(ImageIterator upperleft, ImageIterator lowerright, Accessor a)
{

    int x = 0;
    int y = 0;
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    

    LugImage * img = createLugImage(w, h, 8, 256);
    if (!img)
            return 0;
    
    // create the gray scaled palette
    fill_bw_pallete(img->cmap);

    unsigned char * pindex = img->r;
    ImageIterator yd(upperleft);

    for(y = 0; y < h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(x = 0; x < w; ++x, ++pindex, ++xd.x)
        {
            *pindex = a(xd);
        }
    }
    return img;
}

template <class ImageIterator, class Accessor>
inline
LugImage *
createScalarLugImage(triple<ImageIterator, ImageIterator, Accessor> src)
{
    return createScalarLugImage(src.first, src.second, src.third);
}

/********************************************************/
/*                                                      */
/*                  createRGBLugImage                   */
/*                                                      */
/********************************************************/

/** \brief Create a RGB LugImage from the given RGB image.

    The size of the LugImage is determined by the input image.
    This function uses accessors to read the data. All pixel data types are 
    converted to 'RGBValue<unsigned char>' during conversion. Use the 
    <a href="VIFFrelatedfunctions.html">VIFF</a> image format, if this is not 
    acceptable.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    template <class ImageIterator, class RGBAccessor>
    LugImage *
    createRGBLugImage(ImageIterator upperleft, ImageIterator lowerright, RGBAccessor a)
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    template <class ImageIterator, class RGBAccessor>
    LugImage *
    createRGBLugImage(triple<ImageIterator, ImageIterator, RGBAccessor> src)
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="lug_8hxx-source.html">vigra/lug.hxx</a>"
    
    \code
    BRGBImage img(width, height);
    ...
    
    LugImage * lug = createRGBLugImage(srcImageRange(img));
    
    writeLugImage("rgbimage.gif", lug);
    
    freeLugImage(lug);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    unsigned char value;
    
    RGBAccessor accessor;
                       
    value = accessor.red(upperleft);
    value = accessor.green(upperleft);
    value = accessor.blue(upperleft);
    \endcode
    
*/
template <class ImageIterator, class Accessor>
LugImage *
createRGBLugImage(ImageIterator upperleft, ImageIterator lowerright, Accessor a)
{

    int x = 0;
    int y = 0;
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;

    LugImage * img = createLugImage(w, h, 24, 0);
    if (!img)
            return 0;

    unsigned char * pred = img->r;
    unsigned char * pgreen = img->g;
    unsigned char * pblue = img->b;
    ImageIterator yd(upperleft);

    for(y = 0; y < h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(x = 0; x < w; ++x, ++pred, ++pgreen, ++pblue, ++xd.x)
        {
            *pred   = a.red(xd);
            *pgreen = a.green(xd);
            *pblue  = a.blue(xd);
        }
    }
    return img;
}

template <class ImageIterator, class Accessor>
inline
LugImage *
createRGBLugImage(triple<ImageIterator, ImageIterator, Accessor> src)
{
    return createRGBLugImage(src.first, src.second, src.third);
}

//@}

#endif // IMP_LUG_HXX


