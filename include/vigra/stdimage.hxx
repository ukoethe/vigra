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

#define VIGRA_DEFINE_ITERATORTRAITS(VALUETYPE, ACCESSOR, CONSTACCESSOR) \
    template<> \
    struct IteratorTraits< \
        BasicImageIterator<VALUETYPE, VALUETYPE **> > \
    { \
        typedef BasicImageIterator<VALUETYPE, VALUETYPE **> \
                                                     Iterator; \
        typedef Iterator                             iterator; \
        typedef iterator::iterator_category          iterator_category; \
        typedef iterator::value_type                 value_type; \
        typedef iterator::reference                  reference; \
        typedef iterator::index_reference            index_reference; \
        typedef iterator::pointer                    pointer; \
        typedef iterator::difference_type            difference_type; \
        typedef iterator::row_iterator               row_iterator; \
        typedef iterator::column_iterator            column_iterator; \
        typedef ACCESSOR<VALUETYPE >                 default_accessor; \
        typedef ACCESSOR<VALUETYPE >                 DefaultAccessor; \
    }; \
    template<> \
    struct IteratorTraits< \
        ConstBasicImageIterator<VALUETYPE, VALUETYPE **> > \
    { \
        typedef \
          ConstBasicImageIterator<VALUETYPE, VALUETYPE **> \
                                                     Iterator; \
        typedef Iterator                             iterator; \
        typedef iterator::iterator_category          iterator_category; \
        typedef iterator::value_type                 value_type; \
        typedef iterator::reference                  reference; \
        typedef iterator::index_reference            index_reference; \
        typedef iterator::pointer                    pointer; \
        typedef iterator::difference_type            difference_type; \
        typedef iterator::row_iterator               row_iterator; \
        typedef iterator::column_iterator            column_iterator; \
        typedef CONSTACCESSOR<VALUETYPE >            default_accessor; \
        typedef CONSTACCESSOR<VALUETYPE >            DefaultAccessor; \
    };

/** \addtogroup StandardImageTypes Standard Image Types

    \brief The most common instantiations of the \ref vigra::BasicImage template
*/
//@{

VIGRA_DEFINE_ITERATORTRAITS(unsigned char, StandardValueAccessor, StandardConstValueAccessor)

    /** Byte (8-bit unsigned) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<unsigned char> BImage;

VIGRA_DEFINE_ITERATORTRAITS(short, StandardValueAccessor, StandardConstValueAccessor)


    /** Short integer (16-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<short> SImage;

VIGRA_DEFINE_ITERATORTRAITS(int, StandardValueAccessor, StandardConstValueAccessor)

    /** Integer (32-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<int> IImage;

VIGRA_DEFINE_ITERATORTRAITS(float, StandardValueAccessor, StandardConstValueAccessor)

    /** Float (float) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<float> FImage;

VIGRA_DEFINE_ITERATORTRAITS(double, StandardValueAccessor, StandardConstValueAccessor)

    /** Double (double) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
   
       <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
 */
typedef BasicImage<double> DImage;

VIGRA_DEFINE_ITERATORTRAITS(RGBValue<unsigned char>, RGBAccessor, RGBAccessor)

    /** Byte (3x 8-bit unsigned) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<unsigned char>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<unsigned char> > BRGBImage;

VIGRA_DEFINE_ITERATORTRAITS(RGBValue<int>, RGBAccessor, RGBAccessor)

    /** Integer (3x 32-bit signed) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<int>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<int> > IRGBImage;

VIGRA_DEFINE_ITERATORTRAITS(RGBValue<float>, RGBAccessor, RGBAccessor)

    /** Floating-point (3x float) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<float>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<float> > FRGBImage;

VIGRA_DEFINE_ITERATORTRAITS(RGBValue<double>, RGBAccessor, RGBAccessor)

    /** Double-precision floating-point (3x double) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<double>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<double> > DRGBImage;

#define VIGRA_PIXELTYPE TinyVector<float, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor) 
#undef VIGRA_PIXELTYPE 
#define VIGRA_PIXELTYPE TinyVector<float, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor) 
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor) 
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor) 
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor) 
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor) 
#undef VIGRA_PIXELTYPE

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


#undef VIGRA_DEFINE_ITERATORTRAITS

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

// define traits for BasicImageIterator instanciations that
// were not explicitly defined above
template <class T>
struct IteratorTraits<BasicImageIterator<T, T **> >
{
    typedef BasicImageIterator<T, T **>          Iterator;
    typedef Iterator                             iterator;
    typedef typename iterator::iterator_category iterator_category;
    typedef typename iterator::value_type        value_type;
    typedef typename iterator::reference         reference;
    typedef typename iterator::index_reference   index_reference;
    typedef typename iterator::pointer           pointer;
    typedef typename iterator::difference_type   difference_type;
    typedef typename iterator::row_iterator      row_iterator;
    typedef typename iterator::column_iterator   column_iterator;
    typedef StandardAccessor<T>                  DefaultAccessor; 
    typedef StandardAccessor<T>                  default_accessor; 
};  

template <class T>
struct IteratorTraits<ConstBasicImageIterator<T, T **> >
{
    typedef ConstBasicImageIterator<T, T **> Iterator;
    typedef Iterator                               iterator;
    typedef typename iterator::iterator_category   iterator_category;
    typedef typename iterator::value_type          value_type;
    typedef typename iterator::reference           reference;
    typedef typename iterator::index_reference     index_reference;
    typedef typename iterator::pointer             pointer;
    typedef typename iterator::difference_type     difference_type;
    typedef typename iterator::row_iterator        row_iterator;
    typedef typename iterator::column_iterator     column_iterator;
    typedef StandardConstAccessor<T>               DefaultAccessor; 
    typedef StandardConstAccessor<T>               default_accessor; 
};  

#endif
    
//@}

} // namespace vigra

#endif // VIGRA_STDIMAGE_HXX
