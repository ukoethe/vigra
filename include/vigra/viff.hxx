/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_VIFF_HXX
#define VIGRA_VIFF_HXX

#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/rgbvalue.hxx"

/** \page VIFFImpex Import/export of the Khoros VIFF format

    Normally, you need not call the VIFF functions directly. They are
    available much more conveniently via \ref importImage() and \ref exportImage() 
    
    <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
              \ref VIFFFunctions "VIFF Functions"
        <DD> <em>Read/write/delete VIFF images</em>

        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
            \ref ConvertVIFF
        <DD> <em>Convert VIFF to/from any images supporting ImageIterators</em>
    </DL>

    The VIFF image format was originally defined by the KHOROS public domain 
    image processing environment (KHOROS is no longer public, though). VIFF images are 
    very versatile - we can store
    many different pixel types (byte, integer, float, double, complex etc.) and
    arbitrary many spectral channels (also called 'image bands'). In particular,
    an image with one channel is a gray scale image, 3 channels or 1 channel 
    plus color map represent RGB images. See the KHOROS documentation at 
    <a href="http://www.khoral.com/">http://www.khoral.com/</a> for details.
*/

#include "vigra/viff.h"

namespace vigra {

/** \defgroup ConvertVIFF Convert VIFF Images

    VIFF images files can store byte, int, float, double etc. pixel types
*/
//@{
/********************************************************/
/*                                                      */
/*                     importViffImage                  */
/*                                                      */
/********************************************************/

/** \brief Convert given ViffImage into image specified by iterator range.

    Accessors are used to write the data.    
    This function calls \ref viffToScalarImage() or \ref viffToRGBImage(), depending on 
    the accessor's value_type.

    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importViffImage(ViffImage * viff, ImageIterator iter, Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importViffImage(ViffImage * viff, pair<ImageIterator, Accessor> dest)
    }
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    ViffImage * viff = readViffImage("scalarimage.xv");
    
    vigra::BImage img(viff->row_size, viff->col_size);
    
    vigra::importViffImage(viff, destImage(img));
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    see \ref viffToScalarImage() and \ref viffToRGBImage()
    
    <b> Preconditions:</b>
    
    see \ref viffToScalarImage() and \ref viffToRGBImage()
    
*/
template <class ImageIterator, class Accessor>
inline void
importViffImage(ViffImage * viff, ImageIterator iter, Accessor a)
{
    typedef typename 
        NumericTraits<typename Accessor::value_type>::isScalar
        isScalar;
    importViffImage(viff, iter, a, isScalar());
}

template <class ImageIterator, class Accessor>
inline void
importViffImage(ViffImage * viff, pair<ImageIterator, Accessor> dest)
{
    importViffImage(viff, dest.first, dest.second);
}

template <class ImageIterator, class Accessor>
inline void
importViffImage(ViffImage * viff, ImageIterator iter, Accessor a, VigraTrueType)
{
    viffToScalarImage(viff, iter, a);
}

template <class ImageIterator, class Accessor>
inline void
importViffImage(ViffImage * viff, ImageIterator iter, Accessor a, VigraFalseType)
{
    viffToRGBImage(viff, iter, a);
}

/********************************************************/
/*                                                      */
/*                    viffToScalarImage                 */
/*                                                      */
/********************************************************/

/** \brief Convert single-band ViffImage to scalar image.

    This function uses accessors to write the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        viffToScalarImage(ViffImage * viff, ImageIterator iter, Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        viffToScalarImage(ViffImage * viff, pair<ImageIterator, Accessor> dest)
    }
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    ViffImage * viff = readViffImage("scalarimage.xv");
    
    vigra::BImage img(viff->row_size, viff->col_size);
    
    vigra::viffToScalarImage(viff, destImage(img));
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    <unsigned char, short, long, float, double> value;
    
    Accessor accessor;
               
    accessor.set(value, upperleft);
    \endcode
    
    <b> Preconditions:</b>
    
    \code
    viff->num_data_bands == 1
    viff->data_storage_type != VFF_TYP_COMPLEX
    viff->data_storage_type != VFF_TYP_DCOMPLEX
    viff->data_storage_type != VFF_TYP_BIT
    viff->map_scheme == VFF_MS_NONE
    viff->location_type == VFF_LOC_IMPLICIT
    viff->data_encode_scheme == VFF_DES_RAW      // no compression
    \endcode
    
*/
template <class ImageIterator, class Accessor>
void
viffToScalarImage(ViffImage * viff, ImageIterator iter, Accessor a)
{
    vigra_precondition(viff != 0, 
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "NULL pointer to input data.");
    
    vigra_precondition(viff->num_data_bands == 1,
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Image is multiband - not scalar.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_COMPLEX,
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Image is VFF_TYP_COMPLEX - not scalar.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_DCOMPLEX,
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Image is VFF_TYP_DCOMPLEX - not scalar.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_BIT, 
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Unable to convert VFF_TYP_BIT.");
    
    vigra_precondition(viff->map_scheme == VFF_MS_NONE,
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Can only convert images without map (VFF_MS_NONE).");
    
    vigra_precondition(viff->location_type == VFF_LOC_IMPLICIT,
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Can only convert images with implicit location (VFF_LOC_IMPLICIT).");

    vigra_precondition(viff->data_encode_scheme == VFF_DES_RAW,
             "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
             "Can only convert uncompressed images (VFF_DES_RAW).");

    int w = viff->row_size;
    int h = viff->col_size;
    
    ImageIterator yd(iter);
    
    switch (viff->data_storage_type)
    {
      case VFF_TYP_1_BYTE:
      {
        unsigned char * ps = (unsigned char *)viff->imagedata;
   
        for(int y=0; y<h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++ps, ++xd.x)
            {
            a.set(*ps, xd);
            }
        }
        break;
      }
      case VFF_TYP_2_BYTE:
      {
        short * ps = (short *)viff->imagedata;
   
        for(int y=0; y<h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++ps, ++xd.x)
            {
            a.set(*ps, xd);
            }
        }
        break;
      }
      case VFF_TYP_4_BYTE:
      {
        long * ps = (long *)viff->imagedata;
   
        for(int y=0; y<h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++ps, ++xd.x)
            {
            a.set(*ps, xd);
            }
        }
        break;
      }
      case VFF_TYP_FLOAT:
      {
        float * ps = (float *)viff->imagedata;
   
        for(int y=0; y<h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++ps, ++xd.x)
            {
            a.set(*ps, xd);
            }
        }
        break;
      }
      case VFF_TYP_DOUBLE:
      {
        double * ps = (double *)viff->imagedata;
   
        for(int y=0; y<h; ++y, ++yd.y)
        {
            ImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++ps, ++xd.x)
            {
            a.set(*ps, xd);
            }
        }
        break;
      }
      default:
      {
        // should not happen
        vigra_precondition(0, 
          "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
          "Unknown data storage type.");
      }
    };  
}

template <class ImageIterator, class Accessor>
void
viffToScalarImage(ViffImage * viff, pair<ImageIterator, Accessor> dest)
{
    viffToScalarImage(viff, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                  viffToMultibandImage                */
/*                                                      */
/********************************************************/

/** \brief Convert multi-band ViffImage to multi-band image.

    This function uses \ref VectorComponentAccessor to write the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class VectorComponentAccessor>
        void
        viffToMultibandImage(ViffImage * viff, ImageIterator iter, VectorComponentAccessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class VectorComponentAccessor>
        void
        viffToMultibandImage(ViffImage * viff, pair<ImageIterator, VectorComponentAccessor> dest)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    using namespace vigra;
    
    ViffImage * viff = readViffImage("rgbimage.xv");
    
    BRGBImage img(viff->row_size, viff->col_size);
    
    viffToMultibandImage(viff, 
        destIter(img.upperLeft(), 
         VectorComponentAccessor<BRGBImage::Iterator, BRGBImage::PixelType>());
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    <unsigned char, short, long, float, double> value;
    
    VectorComponentAccessor accessor;
               
    accessor.setCurrentIndex(idx);   // initial setting of index 
                     // is ignored
               
    accessor.set(value, upperleft);
    \endcode
    
    <b> Preconditions:</b>
    
    \code
    viff->data_storage_type != VFF_TYP_COMPLEX
    viff->data_storage_type != VFF_TYP_DCOMPLEX
    viff->data_storage_type != VFF_TYP_BIT
    viff->map_scheme == VFF_MS_NONE
    viff->location_type == VFF_LOC_IMPLICIT
    viff->data_encode_scheme == VFF_DES_RAW
    
    viff->number_of_data_bands <= VECTOR_SIZE  // unchecked !
    \endcode
    
*/
template <class ImageIterator, class VectorComponentAccessor>
void
viffToMultibandImage(ViffImage * viff, ImageIterator iter, VectorComponentAccessor a)
{
    vigra_precondition(viff != 0,
              "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
          "NULL pointer to input data.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_COMPLEX,
             "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
             "Image is VFF_TYP_COMPLEX - not multiband.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_DCOMPLEX,
             "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
             "Image is VFF_TYP_DCOMPLEX - not multiband.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_BIT,
             "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
             "Unable to convert VFF_TYP_BIT.");
    
    vigra_precondition(viff->map_scheme == VFF_MS_NONE,
             "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
             "Can only convert images without map (VFF_MS_NONE).");
    
    vigra_precondition(viff->location_type == VFF_LOC_IMPLICIT,
             "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
             "Can only convert images with implicit location (VFF_LOC_IMPLICIT).");

    vigra_precondition(viff->data_encode_scheme == VFF_DES_RAW,
             "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
             "Can only convert uncompressed images (VFF_DES_RAW).");

    int w = viff->row_size;
    int h = viff->col_size;
    int nb = viff->num_data_bands;
    
    switch (viff->data_storage_type)
    {
      case VFF_TYP_1_BYTE:
      {
    unsigned char * ps = (unsigned char *)viff->imagedata;
    
        for(int b=0; b<nb; ++b)
        {
        ImageIterator yd(iter);
        a.setCurrentIndex(b);
        for(int y=0; y<h; ++y, ++yd.y)
        {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            a.set(*ps, xd);
        }
        }
    }
    break;
      }
      case VFF_TYP_2_BYTE:
      {
        short * ps = (short *)viff->imagedata;
   
        for(int b=0; b<nb; ++b)
        {
        ImageIterator yd(iter);
        a.setCurrentIndex(b);
        for(int y=0; y<h; ++y, ++yd.y)
        {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            a.set(*ps, xd);
        }
        }
    }
    break;
      }
      case VFF_TYP_4_BYTE:
      {
        long * ps = (long *)viff->imagedata;
   
        for(int b=0; b<nb; ++b)
        {
        ImageIterator yd(iter);
        a.setCurrentIndex(b);
        for(int y=0; y<h; ++y, ++yd.y)
        {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            a.set(*ps, xd);
        }
        }
    }
    break;
      }
      case VFF_TYP_FLOAT:
      {
        float * ps = (float *)viff->imagedata;
   
        for(int b=0; b<nb; ++b)
        {
        ImageIterator yd(iter);
        a.setCurrentIndex(b);
        for(int y=0; y<h; ++y, ++yd.y)
        {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            a.set(*ps, xd);
        }
        }
    }
    break;
      }
      case VFF_TYP_DOUBLE:
      {
        double * ps = (double *)viff->imagedata;
   
        for(int b=0; b<nb; ++b)
        {
        ImageIterator yd(iter);
        a.setCurrentIndex(b);
        for(int y=0; y<h; ++y, ++yd.y)
        {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            a.set(*ps, xd);
        }
        }
    }
    break;
      }
      default:
      {
        // should not happen
        vigra_precondition(0,
          "viffToMultibandImage(ViffImage *, VectorImageIterator): " 
          "Unknown data storage type.");
      }
    };  
}

template <class ImageIterator, class VectorComponentAccessor>
void
viffToMultibandImage(ViffImage * viff, pair<ImageIterator, VectorComponentAccessor> dest)
{
    viffToMultibandImage(viff, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                      viffToRGBImage                  */
/*                                                      */
/********************************************************/

/** \brief Convert RGB (3-band or color-mapped) or single-band ViffImage 
    to RGB image.
    
    If the source is a single band image, the result is a RGB image
    with all gray.
    This function uses \ref RGBAccessor to write the data.
    A RGBImageIterator is an iterator which is associated with a
    RGBAccessor.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        viffToRGBImage(ViffImage * viff, RGBImageIterator iter, RGBAccessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        viffToRGBImage(ViffImage * viff, pair<RGBImageIterator, RGBAccessor> dest)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    ViffImage * viff = readViffImage("rgbimage.xv");
    
    vigra::BRGBImage img(viff->row_size, viff->col_size);
    
    vigra::viffToRGBImage(viff, destImage(img));
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    <unsigned char, short, long, float, double> rvalue, gvalue, bvalue;
    
    RGBAccessor accessor;
                           
    accessor.setRed(rvalue, upperleft);
    accessor.setGreen(gvalue, upperleft);
    accessor.setBlue(bvalue, upperleft);
    \endcode
    
    <b> Preconditions:</b>
    
    \code
    viff->data_storage_type != VFF_TYP_COMPLEX
    viff->data_storage_type != VFF_TYP_DCOMPLEX
    viff->data_storage_type != VFF_TYP_BIT
    viff->location_type == VFF_LOC_IMPLICIT
    viff->data_encode_scheme == VFF_DES_RAW
    
    ((viff->map_scheme == VFF_MS_NONE) && 
                         ((viff->num_data_bands == 1) ||
              (viff->num_data_bands == 3)))
    ||
    ((viff->map_scheme == VFF_MS_ONEPERBAND)         &&
                         (viff->num_data_bands == 1) && 
             (viff->map_row_size == 3)   &&
             (viff->data_storage_type == VFF_TYP_1_BYTE) &&
                     (viff->map_storage_type == VFF_MAPTYP_1_BYTE))
    \endcode
    
*/
template <class RGBImageIterator, class RGBAccessor>
void
viffToRGBImage(ViffImage * viff, RGBImageIterator iter, RGBAccessor a)
{
    vigra_precondition(viff != 0,
              "viffToRGBImage(ViffImage *, RGBImageIterator): " 
          "NULL pointer to input data.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_COMPLEX,
             "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Image is VFF_TYP_COMPLEX - not RGB.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_DCOMPLEX,
             "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Image is VFF_TYP_DCOMPLEX - not RGB.");
    
    vigra_precondition(viff->data_storage_type != VFF_TYP_BIT,
             "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Unable to convert VFF_TYP_BIT.");
    
    vigra_precondition(viff->location_type == VFF_LOC_IMPLICIT,
             "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Can only convert images with implicit location (VFF_LOC_IMPLICIT).");

    vigra_precondition(viff->data_encode_scheme == VFF_DES_RAW,
             "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Can only convert uncompressed images (VFF_DES_RAW).");

    if(viff->map_scheme == VFF_MS_NONE)
    {
        int w = viff->row_size;
        int h = viff->col_size;
        int bandoffset;

        if(viff->num_data_bands == 3)
        {    
           bandoffset = w*h;
        }
        else if(viff->num_data_bands == 1)
        {
            bandoffset = 0;
        }
        else 
          vigra_precondition(0,
                 "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Wrong number of data bands (must be 3 or 1).");

        RGBImageIterator yd(iter);

        switch (viff->data_storage_type)
        {
          case VFF_TYP_1_BYTE:
          {
            unsigned char * pr = (unsigned char *)viff->imagedata;
            unsigned char * pg = pr + bandoffset;
            unsigned char * pb = pg + bandoffset;

            for(int y=0; y<h; ++y, ++yd.y)
            {
            RGBImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
            {
                a.setRed(*pr, xd);
                a.setGreen(*pg, xd);
                a.setBlue(*pb, xd);
            }
            }
            break;
          }
          case VFF_TYP_2_BYTE:
          {
            short * pr = (short *)viff->imagedata;
            short * pg = pr + bandoffset;
            short * pb = pg + bandoffset;

            for(int y=0; y<h; ++y, ++yd.y)
            {
            RGBImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
            {
                a.setRed(*pr, xd);
                a.setGreen(*pg, xd);
                a.setBlue(*pb, xd);
            }
            }
            break;
          }
          case VFF_TYP_4_BYTE:
          {
            long * pr = (long *)viff->imagedata;
            long * pg = pr + bandoffset;
            long * pb = pg + bandoffset;

            for(int y=0; y<h; ++y, ++yd.y)
            {
            RGBImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
            {
                a.setRed(*pr, xd);
                a.setGreen(*pg, xd);
                a.setBlue(*pb, xd);
            }
            }
            break;
          }
          case VFF_TYP_FLOAT:
          {
            float * pr = (float *)viff->imagedata;
            float * pg = pr + bandoffset;
            float * pb = pg + bandoffset;

            for(int y=0; y<h; ++y, ++yd.y)
            {
            RGBImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
            {
                a.setRed(*pr, xd);
                a.setGreen(*pg, xd);
                a.setBlue(*pb, xd);
            }
            }
            break;
          }
          case VFF_TYP_DOUBLE:
          {
            double * pr = (double *)viff->imagedata;
            double * pg = pr + bandoffset;
            double * pb = pg + bandoffset;

            for(int y=0; y<h; ++y, ++yd.y)
            {
            RGBImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
            {
                a.setRed(*pr, xd);
                a.setGreen(*pg, xd);
                a.setBlue(*pb, xd);
            }
            }
            break;
          }
          default:
          {
            // should not happen
            vigra_precondition(0,
                  "viffToScalarImage(ViffImage *, ScalarImageIterator): " 
              "Unknown data storage type.");
          }
        }
    }
    else if(viff->map_scheme == VFF_MS_ONEPERBAND)
    {
        vigra_precondition((viff->num_data_bands == 1) && (viff->map_row_size == 3),
                 "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Mapped image must have 1 data band and 3 map columns.");
             
        vigra_precondition((viff->data_storage_type == VFF_TYP_1_BYTE) &&
                     (viff->map_storage_type == VFF_MAPTYP_1_BYTE),
                 "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Can only convert mapped images with VFF_TYP_1_BYTE"
             " and VFF_MAPTYP_1_BYTE. Sorry.");
        
        unsigned char * red = (unsigned char *)viff->maps;
        unsigned char * green = red + viff->map_col_size;
        unsigned char * blue = green + viff->map_col_size;

        int w = viff->row_size;
        int h = viff->col_size;
    
        unsigned char * ps = (unsigned char *)viff->imagedata;

        RGBImageIterator yd(iter);
        for(int y=0; y<h; ++y, ++yd.y)
        {
            RGBImageIterator xd(yd);
            for(int x=0; x < w; ++x, ++ps, ++xd.x)
            {
                a.setRed(red[*ps], xd);
                a.setGreen(green[*ps], xd);
                a.setBlue(blue[*ps], xd);
            }
        }
    }
    else       
      vigra_precondition(0,
             "viffToRGBImage(ViffImage *, RGBImageIterator): " 
             "Unable to convert this kind of map.");
}

template <class RGBImageIterator, class RGBAccessor>
inline void
viffToRGBImage(ViffImage * viff, pair<RGBImageIterator, RGBAccessor> dest)
{
    viffToRGBImage(viff, dest.first, dest.second);
}

template <class T>
struct CreateViffImage;

/********************************************************/
/*                                                      */
/*                     createViffImage                  */
/*                                                      */
/********************************************************/

/** \brief Create a ViffImage from the given iterator range.

    Type and size of the ViffImage are determined by the input image. 
    Currently, the function can create scalar images of type 
    unsigned char, short, int, float, and double, and RGB images
    of type unsigned char, int, and float.
    This function uses accessors to read the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline ViffImage *
        createViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                        Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline ViffImage *
        createViffImage(triple<ImageIterator, ImageIterator, Accessor> src)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage img(width, height);
    
    ...
    
    ViffImage * viff = vigra::createViffImage(srcImageRange(img));
    
    writeViffImage("output.xv", viff);
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    Accessor accessor;
                           
    accessor(upperleft);   // result written into ViffImage
    \endcode
*/
template <class ImageIterator, class Accessor>
inline ViffImage *
createViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
{
    return CreateViffImage<typename Accessor::value_type>::
        exec(upperleft, lowerright, a);
}

template <class ImageIterator, class Accessor>
inline ViffImage *
createViffImage(triple<ImageIterator, ImageIterator, Accessor> src)
{
    return createViffImage(src.first, src.second, src.third);
}

/********************************************************/
/*                                                      */
/*                createScalarViffImage                 */
/*                                                      */
/********************************************************/

/** \brief Create a single-band ViffImage from the given scalar image.

    Type and size of the ViffImage are determined by the input image 
    (may be one of unsigned char, short, int, float, or double).
    This function uses accessors to read the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline ViffImage *
        createScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                  Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline ViffImage *
        createScalarViffImage(triple<ImageIterator, ImageIterator, Accessor> src)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage img(width, height);
    
    ...
    
    ViffImage * viff = vigra::createScalarViffImage(srcImageRange(img));
    
    writeViffImage("output.xv", viff);
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    Accessor accessor;
                           
    accessor(upperleft);   // result written into ViffImage
    \endcode
    
*/
template <class ImageIterator, class Accessor>
inline ViffImage *
createScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
{
    return CreateViffImage<typename Accessor::value_type>::
        exec(upperleft, lowerright, a);
}

template <class ImageIterator, class Accessor>
inline ViffImage *
createScalarViffImage(triple<ImageIterator, ImageIterator, Accessor> src)
{
    return createScalarViffImage(src.first, src.second, src.third);
}

template <class ImageIterator, class Accessor>
ViffImage *
createBScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createsimpleimage(h, w, VFF_TYP_1_BYTE);
    vigra_postcondition(viff != 0, "createScalarViffImage(): Unable to allocate memory");
    
    unsigned char * ps = (unsigned char *)viff->imagedata;

    ImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            *ps = a(xd);
        }
    }
    
    return viff;
}

template <class ImageIterator, class Accessor>
ViffImage *
createShortScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                   Accessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createsimpleimage(h, w, VFF_TYP_2_BYTE);
    vigra_postcondition(viff != 0, "createScalarViffImage(): Unable to allocate memory");
    
    short * ps = (short *)viff->imagedata;

    ImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            *ps = a(xd);
        }
    }
    
    return viff;
}

template <class ImageIterator, class Accessor>
ViffImage *
createIScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                   Accessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createsimpleimage(h, w, VFF_TYP_4_BYTE);
    vigra_postcondition(viff != 0, "createScalarViffImage(): Unable to allocate memory");
    
    int * ps = (int *)viff->imagedata;

    ImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            *ps = a(xd);
        }
    }
    
    return viff;
}

template <class ImageIterator, class Accessor>
ViffImage *
createFScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                   Accessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createsimpleimage(h, w, VFF_TYP_FLOAT);
    
    float * ps = (float *)viff->imagedata;

    ImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            *ps = a(xd);
        }
    }
    
    return viff;
}

template <class ImageIterator, class Accessor>
ViffImage *
createDScalarViffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                   Accessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createsimpleimage(h, w, VFF_TYP_DOUBLE);
    vigra_postcondition(viff != 0, "createScalarViffImage(): Unable to allocate memory");
        
    double * ps = (double *)viff->imagedata;

    ImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        ImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++ps, ++xd.x)
        {
            *ps = a(xd);
        }
    }
    
    return viff;
}

template <>
struct CreateViffImage<unsigned char>
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createBScalarViffImage(upperleft, lowerright, a);
    }
};

template <>
struct CreateViffImage<short>
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createShortScalarViffImage(upperleft, lowerright, a);
    }
};

template <>
struct CreateViffImage<int>
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createIScalarViffImage(upperleft, lowerright, a);
    }
};

template <>
struct CreateViffImage<float>
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createFScalarViffImage(upperleft, lowerright, a);
    }
};

template <>
struct CreateViffImage<double>
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createDScalarViffImage(upperleft, lowerright, a);
    }
};

/********************************************************/
/*                                                      */
/*                  createRGBViffImage                  */
/*                                                      */
/********************************************************/

/** \brief Create a 3-band ViffImage from the given RGB image.

    Type and size of the ViffImage are determined by the input image 
    (may be one of unsigned char, int, or float).
    This function uses \ref RGBAccessor to read the data. A
    RGBImageIterator is an iterator that is associated with a
    RGBAccessor.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        inline ViffImage *
        createRGBViffImage(RGBImageIterator upperleft, RGBImageIterator lowerright,
                   RGBAccessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        inline ViffImage *
        createRGBViffImage(triple<RGBImageIterator, RGBImageIterator, RGBAccessor> src)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="viff_8hxx-source.html">vigra/viff.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BRGBImage img(width, height);
    
    ...
    
    ViffImage * viff = vigra::createRGBViffImage(srcImageRange(img));
    
    writeViffImage("output.xv", viff);
    
    freeViffImage(viff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    RGBAccessor accessor;
                           
    accessor.red(upperleft);     // result written into ViffImage
    accessor.green(upperleft);   // result written into ViffImage
    accessor.blue(upperleft);    // result written into ViffImage
    \endcode
    
*/
template <class RGBImageIterator, class RGBAccessor>
inline ViffImage *
createRGBViffImage(RGBImageIterator upperleft, RGBImageIterator lowerright,
                   RGBAccessor a)
{
    return CreateViffImage<typename RGBAccessor::value_type>::
        exec(upperleft, lowerright, a);
//    return createRGBViffImage(upperleft, lowerright, a, a(upperleft));
}

template <class RGBImageIterator, class RGBAccessor>
inline ViffImage *
createRGBViffImage(triple<RGBImageIterator, RGBImageIterator, RGBAccessor> src)
{
    return createRGBViffImage(src.first, src.second, src.third);
}

template <class RGBImageIterator, class RGBAccessor>
ViffImage *
createBRGBViffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createmultibandimage(h, w, VFF_TYP_1_BYTE, 3);
    vigra_postcondition(viff != 0, "createRGBViffImage(): Unable to allocate memory");
    
    viff->color_space_model = VFF_CM_genericRGB;
    
    unsigned char * pr = (unsigned char *)viff->imagedata;
    unsigned char * pg = pr + w*h;
    unsigned char * pb = pg + w*h;

    RGBImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        RGBImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
        {
            *pr = a.red(xd);
            *pg = a.green(xd);
            *pb = a.blue(xd);
        }
    }
    
    return viff;
}

template <class RGBImageIterator, class RGBAccessor>
ViffImage *
createIRGBViffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                       RGBAccessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createmultibandimage(h, w, VFF_TYP_4_BYTE, 3);
    vigra_postcondition(viff != 0, "createRGBViffImage(): Unable to allocate memory");

    viff->color_space_model = VFF_CM_genericRGB;
    
    int * pr = (int *)viff->imagedata;
    int * pg = pr + w*h;
    int * pb = pg + w*h;

    RGBImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        RGBImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
        {
            *pr = a.red(xd);
            *pg = a.green(xd);
            *pb = a.blue(xd);
        }
    }
    
    return viff;
}

template <class RGBImageIterator, class RGBAccessor>
ViffImage *
createFRGBViffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    ViffImage * viff = createmultibandimage(h, w, VFF_TYP_FLOAT, 3);
    vigra_postcondition(viff != 0, "createRGBViffImage(): Unable to allocate memory");
    
    viff->color_space_model = VFF_CM_genericRGB;
    
    float * pr = (float *)viff->imagedata;
    float * pg = pr + w*h;
    float * pb = pg + w*h;

    RGBImageIterator yd(upperleft);
    for(int y=0; y<h; ++y, ++yd.y)
    {
        RGBImageIterator xd(yd);
        for(int x=0; x < w; ++x, ++pr, ++pg, ++pb, ++xd.x)
        {
            *pr = a.red(xd);
            *pg = a.green(xd);
            *pb = a.blue(xd);
        }
    }
    
    return viff;
}

template <>
struct CreateViffImage<RGBValue<unsigned char> >
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createBRGBViffImage(upperleft, lowerright, a);
    }
};

template <>
struct CreateViffImage<RGBValue<int> >
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createIRGBViffImage(upperleft, lowerright, a);
    }
};

template <>
struct CreateViffImage<RGBValue<float> >
{
    template <class ImageIterator, class Accessor>
    static ViffImage *
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a)
    {
        return createFRGBViffImage(upperleft, lowerright, a);
    }
};

//@}

} // namespace vigra

#endif // VIGRA_VIFF_HXX
