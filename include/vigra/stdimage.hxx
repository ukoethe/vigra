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


/***********************************************************/

/** \page ArgumentObjectFactories Argument Object Factories
    
    Factory functions to create argument objects which simplify long argument lists.

    <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"
    
    <DL>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref ImageBasedArgumentObjectFactories
        <DD>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref IteratorBasedArgumentObjectFactories
        <DD>
    </DL>

    Long argument lists provide for greater flexibility of functions,
    but they are also tedious and error prone, when we don't need
    the flexibility. Thus, we define argument objects which
    automatically provide reasonable defaults for those arguments that we 
    didn't specify explicitly. 
    
    The argument objects are created via a number of factory functions.
    Since these functions have descriptive names, they also serve
    to improve readability: the name of each factory tells te purpose of its
    argument object. 
    
    Consider the following example. Without argument objects we had to 
    write something like this (cf. \ref copyImageIf()):
    
    \code
    vigra::BImage img1, img2, img3;
    
    // fill img1 and img2 ...
    
    vigra::copyImageIf(img1.upperLeft(), img1.lowerRight(), img1.accessor(),
                img2.upperLeft(), img2.accessor(),
                img3.upperLeft(), img3.accessor());
    \endcode
    
    Using the argument object factories, this becomes much shorter and
    more readable:
    
    \code
    vigra::copyImageIf(srcImageRange(img1),
                maskImage(img2),
                destImage(img3));
    \endcode
    
    The names of the factories clearly tell which image is source, mask, 
    and destination. In addition, the suffix <TT>Range</TT> must be used 
    for those argument objects that need to specify the lower right
    corner of the region of interest. Typically, this is only the first
    source argument, but sometimes the first destiniation argument must 
    also contain a range.
    
    The factory functions come in two flavours: Iterator based and 
    image based factories. Above we have seen the image based variant.
    The iterator based variant would look like this:
    
    \code
    vigra::copyImageIf(srcIterRange(img1.upperLeft(), img1.lowerRight()),
                maskIter(img2.upperLeft()),
                destIter(img3.upperLeft()));
    \endcode
    
    These factory functions contain the word <TT>Iter</TT> instead of the word 
    <TT>Image</TT>,  They would normally be used if we couldn't access the 
    images (for example, within a function which got passed iterators)
    or if we didn't want to operate on the entire image. The default 
    accessor is obtained via \ref vigra::IteratorTraits.
    
    All factory functions also allow to specify accessors explicitly. This
    is useful if we can't use the default accessor. This variant looks 
    like this:
    
    \code
    vigra::copyImageIf(srcImageRange(img1),
                maskImage(img2, MaskPredicateAccessor()),
                destImage(img3));    
    \endcode
    
    or
    
    \code
    vigra::copyImageIf(srcIterRange(img1.upperLeft(), img1.lowerRight()),
                maskIter(img2.upperLeft(), MaskPredicateAccessor()),
                destIter(img3.upperLeft()));
    \endcode
    
    All versions can be mixed freely within one explession.
    Technically, the argument objects are simply defined as 
    pairs and triples of iterators and accessor so that all algorithms 
    should declare a call interface version based on pairs and triples 
    (see for example \ref copyImageIf()).

  \section ImageBasedArgumentObjectFactories Image Based Argument Object Factories
        
    These factories can be used to create argument objects when we 
    are given instances or subclasses of \ref vigra::BasicImage (see
    \ref StandardImageTypes for instances defined per default).
    These factory functions access <TT>img.upperLeft()</TT>, 
    <TT>img.lowerRight()</TT>, and <TT>img.accessor()</TT> to obtain the iterators
    and accessor for the given image (unless the accessor is 
    given explicitly). The following factory functions are provided:
    
    <table>
    <tr><td>
        \htmlonly
        <th bgcolor="#f0e0c0" colspan=2 align=left>
        \endhtmlonly
        <TT>\ref vigra::BasicImage "vigra::BasicImage<SomeType>" img;</TT>
        \htmlonly
        </th>
        \endhtmlonly
    </td></tr>
    <tr><td>
        
    <TT>srcImageRange(img)</TT>
    </td><td>
        create argument object containing upper left, lower right, and
        default accessor of source image
        
    </td></tr>
    <tr><td>
        
    <TT>srcImageRange(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left, lower right
        of source image, and given accessor
        
    </td></tr>
    <tr><td>
        
    <TT>srcImage(img)</TT>
    </td><td>
        create argument object containing upper left, and
        default accessor of source image
        
    </td></tr>
    <tr><td>
        
    <TT>srcImage(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left
        of source image, and given accessor
        
    </td></tr>
    <tr><td>
    
    <TT>maskImage(img)</TT>
    </td><td>
        create argument object containing upper left, and
        default accessor of mask image
        
    </td></tr>
    <tr><td>
        
    <TT>maskImage(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left
        of mask image, and given accessor
        
    </td></tr>
    <tr><td>
    
    <TT>destImageRange(img)</TT>
    </td><td>
        create argument object containing upper left, lower right, and
        default accessor of destination image
        
    </td></tr>
    <tr><td>
        
    <TT>destImageRange(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left, lower right
        of destination image, and given accessor
        
    </td></tr>
    <tr><td>
        
    <TT>destImage(img)</TT>
    </td><td>
        create argument object containing upper left, and
        default accessor of destination image
        
    </td></tr>
    <tr><td>
        
    <TT>destImage(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left
        of destination image, and given accessor
        
    </td></tr>
    </table>


  \section IteratorBasedArgumentObjectFactories Iterator Based Argument Object Factories
        
    These factories can be used to create argument objects when we 
    are given \ref ImageIterators.
    These factory functions use \ref vigra::IteratorTraits to
    get the default accessor for the given iterator unless the 
    accessor is given explicitly. The following factory functions 
    are provided:
    
    <table>
    <tr><td>
        \htmlonly
        <th bgcolor="#f0e0c0" colspan=2 align=left>
        \endhtmlonly
        <TT>\ref vigra::BasicImage::Iterator "vigra::BasicImage<SomeType>::Iterator" i1, i2;</TT>
        \htmlonly
        </th>
        \endhtmlonly
    </td></tr>
    <tr><td>
        
    <TT>srcIterRange(i1, i2)</TT>
    </td><td>
        create argument object containing the given iterators and
        corresponding default accessor (for source image)
        
    </td></tr>
    <tr><td>
        
    <TT>srcIterRange(i1, i2, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterators and
        accessor (for source image)
        
    </td></tr>
    <tr><td>
        
    <TT>srcIter(i1)</TT>
    </td><td>
        create argument object containing the given iterator and
        corresponding default accessor (for source image)
        
    </td></tr>
    <tr><td>
        
    <TT>srcIter(i1, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterator and
        accessor (for source image)
        
    </td></tr>
    <tr><td>
    
    <TT>maskIter(i1)</TT>
    </td><td>
        create argument object containing the given iterator and
        corresponding default accessor (for mask image)
        
    </td></tr>
    <tr><td>
        
    <TT>maskIter(i1, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterator and
        accessor (for mask image)
        
    </td></tr>
    <tr><td>
    
    <TT>destIterRange(i1, i2)</TT>
    </td><td>
        create argument object containing the given iterators and
        corresponding default accessor (for destination image)
        
    </td></tr>
    <tr><td>
        
    <TT>destIterRange(i1, i2, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterators and
        accessor (for destination image)
        
    </td></tr>
    <tr><td>
        
    <TT>destIter(i1)</TT>
    </td><td>
        create argument object containing the given iterator and
        corresponding default accessor (for destination image)
        
    </td></tr>
    <tr><td>
        
    <TT>destIter(i1, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterator and
        accessor (for destination image)
        
    </td></tr>
    </table>
*/

template <class Iterator, class Accessor>
inline triple<Iterator, Iterator, Accessor>
srcIterRange(Iterator upperleft, Iterator lowerright, Accessor a)
{
    return triple<Iterator, Iterator, Accessor>(upperleft, lowerright, a);
}

template <class Iterator, class Accessor>
inline pair<Iterator, Accessor>
srcIter(Iterator upperleft, Accessor a)
{
    return pair<Iterator, Accessor>(upperleft, a);
}

template <class Iterator, class Accessor>
inline pair<Iterator, Accessor>
maskIter(Iterator upperleft, Accessor a)
{
    return pair<Iterator, Accessor>(upperleft, a);
}

template <class Iterator, class Accessor>
inline pair<Iterator, Accessor>
destIter(Iterator upperleft, Accessor a)
{
    return pair<Iterator, Accessor>(upperleft, a);
}


template <class Iterator, class Accessor>
inline triple<Iterator, Iterator, Accessor>
destIterRange(Iterator upperleft, Iterator lowerright, Accessor a)
{
    return triple<Iterator, Iterator, Accessor>(upperleft, lowerright, a);
}

/****************************************************************/

template <class PixelType, class Accessor>
inline triple<typename BasicImage<PixelType>::const_traverser, 
              typename BasicImage<PixelType>::const_traverser, Accessor>
srcImageRange(BasicImage<PixelType> const & img, Accessor a)
{
    return triple<typename BasicImage<PixelType>::const_traverser, 
                  typename BasicImage<PixelType>::const_traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::const_traverser, Accessor>
srcImage(BasicImage<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline triple<typename BasicImage<PixelType>::traverser, 
              typename BasicImage<PixelType>::traverser, Accessor>
destImageRange(BasicImage<PixelType> & img, Accessor a)
{
    return triple<typename BasicImage<PixelType>::traverser, 
                  typename BasicImage<PixelType>::traverser, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::traverser, Accessor>
destImage(BasicImage<PixelType> & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::const_traverser, Accessor>
maskImage(BasicImage<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
}

/****************************************************************/

template <class PixelType>
inline triple<typename BasicImage<PixelType>::const_traverser, 
              typename BasicImage<PixelType>::const_traverser, 
          typename BasicImage<PixelType>::ConstAccessor>
srcImageRange(BasicImage<PixelType> const & img)
{
    return triple<typename BasicImage<PixelType>::const_traverser, 
                  typename BasicImage<PixelType>::const_traverser, 
          typename BasicImage<PixelType>::ConstAccessor>(img.upperLeft(),
                                        img.lowerRight(),
                        img.accessor());
}

template <class PixelType>
inline pair< typename BasicImage<PixelType>::const_traverser, 
             typename BasicImage<PixelType>::ConstAccessor>
srcImage(BasicImage<PixelType> const & img)
{
    return pair<typename BasicImage<PixelType>::const_traverser, 
                typename BasicImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                         img.accessor());
}

template <class PixelType>
inline triple< typename BasicImage<PixelType>::traverser, 
               typename BasicImage<PixelType>::traverser, 
           typename BasicImage<PixelType>::Accessor>
destImageRange(BasicImage<PixelType> & img)
{
    return triple<typename BasicImage<PixelType>::traverser, 
                  typename BasicImage<PixelType>::traverser, 
          typename BasicImage<PixelType>::Accessor>(img.upperLeft(),
                                        img.lowerRight(),
                        img.accessor());
}

template <class PixelType>
inline pair< typename BasicImage<PixelType>::traverser, 
             typename BasicImage<PixelType>::Accessor>
destImage(BasicImage<PixelType> & img)
{
    return pair<typename BasicImage<PixelType>::traverser, 
                typename BasicImage<PixelType>::Accessor>(img.upperLeft(), 
                                         img.accessor());
}

template <class PixelType>
inline pair< typename BasicImage<PixelType>::const_traverser, 
             typename BasicImage<PixelType>::ConstAccessor>
maskImage(BasicImage<PixelType> const & img)
{
    return pair<typename BasicImage<PixelType>::const_traverser, 
                typename BasicImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                         img.accessor());
}

/****************************************************************/

template <class Iterator>
inline pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
srcIter(Iterator upperleft)
{
    return pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft,
              IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline triple<Iterator, Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
srcIterRange(Iterator upperleft, Iterator lowerright)
{
    return triple<Iterator, Iterator, 
                  typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft, lowerright, 
              IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
maskIter(Iterator upperleft)
{
    return pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft,
              IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
destIter(Iterator upperleft)
{
    return pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft,
              IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline triple<Iterator, Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
destIterRange(Iterator upperleft, Iterator lowerright)
{
    return triple<Iterator, Iterator, 
                  typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft, lowerright, 
              IteratorTraits<Iterator>::DefaultAccessor());
}

} // namespace vigra

#endif // VIGRA_STDIMAGE_HXX
