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
 
 
#ifndef VIGRA_STDIMAGE_HXX
#define VIGRA_STDIMAGE_HXX

#include "vigra/tuple.hxx"
#include "vigra/basicimage.hxx"
#include "vigra/iteratortraits.hxx"
#include "vigra/accessor.hxx"
#include "vigra/rgbvalue.hxx"

namespace vigra { 

/** \addtogroup StandardImageTypes Standard Image Types

    \brief The most common instantiations of the \ref vigra::BasicImage template
*/
//@{

    /** Byte (8-bit unsigned) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<unsigned char> BImage;



    /** Short integer (16-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<short> SImage;


    /** Integer (32-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<int> IImage;


    /** Float (float) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<float> FImage;


    /** Double (double) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
   
       <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
 */
typedef BasicImage<double> DImage;


    /** Byte (3x 8-bit unsigned) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<unsigned char>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<unsigned char> > BRGBImage;


    /** Integer (3x 32-bit signed) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<int>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<int> > IRGBImage;


    /** Floating-point (3x float) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<float>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<float> > FRGBImage;


    /** Double-precision floating-point (3x double) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<double>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<double> > DRGBImage;

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "TinyVector<float, 2>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<float, 2> > FVector2Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "TinyVector<float, 3>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<float, 3> > FVector3Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "TinyVector<float, 4>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<float, 4> > FVector4Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "TinyVector<double, 2>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<double, 2> > DVector2Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "TinyVector<double, 3>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
//typedef BasicImage<TinyVector<double, 3> > DVector3Image; 
typedef BasicImage<TinyVector<double, 3> > DVector3Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "TinyVector<double, 4>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<double, 4> > DVector4Image; 

//@}

} // namespace vigra

#endif // VIGRA_STDIMAGE_HXX
