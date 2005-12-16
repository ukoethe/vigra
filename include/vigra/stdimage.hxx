/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */                
/*                                                                      */
/************************************************************************/
 
 
#ifndef VIGRA_STDIMAGE_HXX
#define VIGRA_STDIMAGE_HXX

#include "vigra/sized_int.hxx"
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
typedef BasicImage<UInt8> BImage;



    /** Short integer (16-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<Int16> SImage;


    /** Integer (32-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<Int32> IImage;


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
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<vigra::UInt8>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<UInt8> > BRGBImage;


    /** Short (3x 16-bit signed) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<vigra::Int16>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<Int16> > SRGBImage;


    /** Integer (3x 32-bit signed) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<vigra::Int32>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<Int32> > IRGBImage;


    /** Floating-point (3x float) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<float>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<float> > FRGBImage;


    /** Double-precision floating-point (3x double) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<double>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<double> > DRGBImage;

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "vigra::TinyVector<float, 2>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<float, 2> > FVector2Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "vigra::TinyVector<float, 3>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<float, 3> > FVector3Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "vigra::TinyVector<float, 4>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<float, 4> > FVector4Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "vigra::TinyVector<double, 2>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<double, 2> > DVector2Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "vigra::TinyVector<double, 3>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
//typedef BasicImage<TinyVector<double, 3> > DVector3Image; 
typedef BasicImage<TinyVector<double, 3> > DVector3Image; 

    /** Floating-point TinyVector image.
        The pixel type is \ref vigra::TinyVector "vigra::TinyVector<double, 4>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::VectorAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<TinyVector<double, 4> > DVector4Image; 

//@}

} // namespace vigra

#endif // VIGRA_STDIMAGE_HXX
