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
 
 
#ifndef VIGRA_STDIMAGE_HXX
#define VIGRA_STDIMAGE_HXX

#include "vigra/basicimage.hxx"
#include "vigra/iteratortraits.hxx"
#include "vigra/accessor.hxx"
#include "vigra/rgbvalue.hxx"

namespace vigra { 

/** \defgroup StandardImageTypes Standard Image Types

    \brief The most common instantiations of the \ref vigra::BasicImage template
*/
//@{
template<>
struct IteratorTraits<BasicImageIterator<unsigned char, unsigned char **> >
{
    typedef StandardAccessor<unsigned char> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<unsigned char, unsigned char **> >
{
    typedef StandardConstAccessor<unsigned char> DefaultAccessor;
};  

    /** Byte (8-bit unsigned) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<unsigned char> BImage;

template<>
struct IteratorTraits<BasicImageIterator<short, short **> >
{
    typedef StandardAccessor<short> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<short, short **> >
{
    typedef StandardConstAccessor<short> DefaultAccessor;
};  

    /** Short integer (16-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<short> SImage;

template<>
struct IteratorTraits<BasicImageIterator<int, int **> >
{
    typedef StandardAccessor<int> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<int, int **> >
{
    typedef StandardConstAccessor<int> DefaultAccessor;
};  

    /** Integer (32-bit signed) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<int> IImage;

template<>
struct IteratorTraits<BasicImageIterator<float, float **> >
{
    typedef StandardAccessor<float> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<float, float **> >
{
    typedef StandardConstAccessor<float> DefaultAccessor;
};  

    /** Float (float) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<float> FImage;

template<>
struct IteratorTraits<BasicImageIterator<double, double **> >
{
    typedef StandardAccessor<double> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<double, double **> >
{
    typedef StandardConstAccessor<double> DefaultAccessor;
};  

    /** Double (double) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
   
       <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
 */
typedef BasicImage<double> DImage;

template<>
struct IteratorTraits<BasicImageIterator<RGBValue<unsigned char>, RGBValue<unsigned char> **> >
{
    typedef RGBAccessor<RGBValue<unsigned char> > DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<RGBValue<unsigned char>, RGBValue<unsigned char> **> >
{
    typedef RGBAccessor<RGBValue<unsigned char> > DefaultAccessor;
};  

    /** Byte (3x 8-bit unsigned) RGB image.
        The pixel type is \ref vigra::RGBValue "vigra::RGBValue<unsigned char>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<unsigned char> > BRGBImage;

template<>
struct IteratorTraits<BasicImageIterator<RGBValue<int>, RGBValue<int> ** > >
{
    typedef RGBAccessor<RGBValue<int> > DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<RGBValue<int>, RGBValue<int> ** > >
{
    typedef RGBAccessor<RGBValue<int> > DefaultAccessor;
};  

    /** Integer (3x 32-bit signed) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<int>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<int> > IRGBImage;

template<>
struct IteratorTraits<BasicImageIterator<RGBValue<float>, RGBValue<float> **> >
{
    typedef RGBAccessor<RGBValue<float> > DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<RGBValue<float>, RGBValue<float> ** > >
{
    typedef RGBAccessor<RGBValue<float> > DefaultAccessor;
};  

    /** Floating-point (3x float) RGB image.
        The pixel type is \ref vigra::RGBValue "RGBValue<float>".
        It uses \ref vigra::BasicImageIterator and \ref vigra::RGBAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="stdimage_8hxx-source.html">vigra/stdimage.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<RGBValue<float> > FRGBImage;

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T>
struct IteratorTraits<BasicImageIterator<T, T **> >
{
    typedef StandardAccessor<T> DefaultAccessor;
};  

template<class T>
struct IteratorTraits<ConstBasicImageIterator<T, T **> >
{
    typedef StandardConstAccessor<T> DefaultAccessor;
    
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
inline triple<typename BasicImage<PixelType>::ConstIterator, 
              typename BasicImage<PixelType>::ConstIterator, Accessor>
srcImageRange(BasicImage<PixelType> const & img, Accessor a)
{
    return triple<typename BasicImage<PixelType>::ConstIterator, 
                  typename BasicImage<PixelType>::ConstIterator, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::ConstIterator, Accessor>
srcImage(BasicImage<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::ConstIterator, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline triple<typename BasicImage<PixelType>::Iterator, 
              typename BasicImage<PixelType>::Iterator, Accessor>
destImageRange(BasicImage<PixelType> & img, Accessor a)
{
    return triple<typename BasicImage<PixelType>::Iterator, 
                  typename BasicImage<PixelType>::Iterator, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::Iterator, Accessor>
destImage(BasicImage<PixelType> & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::Iterator, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::ConstIterator, Accessor>
maskImage(BasicImage<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::ConstIterator, 
                Accessor>(img.upperLeft(), a);
}

/****************************************************************/

template <class PixelType>
inline triple<typename BasicImage<PixelType>::ConstIterator, 
              typename BasicImage<PixelType>::ConstIterator, 
          typename BasicImage<PixelType>::ConstAccessor>
srcImageRange(BasicImage<PixelType> const & img)
{
    return triple<typename BasicImage<PixelType>::ConstIterator, 
                  typename BasicImage<PixelType>::ConstIterator, 
          typename BasicImage<PixelType>::ConstAccessor>(img.upperLeft(),
                                        img.lowerRight(),
                        img.accessor());
}

template <class PixelType>
inline pair< typename BasicImage<PixelType>::ConstIterator, 
             typename BasicImage<PixelType>::ConstAccessor>
srcImage(BasicImage<PixelType> const & img)
{
    return pair<typename BasicImage<PixelType>::ConstIterator, 
                typename BasicImage<PixelType>::ConstAccessor>(img.upperLeft(), 
                                         img.accessor());
}

template <class PixelType>
inline triple< typename BasicImage<PixelType>::Iterator, 
               typename BasicImage<PixelType>::Iterator, 
           typename BasicImage<PixelType>::Accessor>
destImageRange(BasicImage<PixelType> & img)
{
    return triple<typename BasicImage<PixelType>::Iterator, 
                  typename BasicImage<PixelType>::Iterator, 
          typename BasicImage<PixelType>::Accessor>(img.upperLeft(),
                                        img.lowerRight(),
                        img.accessor());
}

template <class PixelType>
inline pair< typename BasicImage<PixelType>::Iterator, 
             typename BasicImage<PixelType>::Accessor>
destImage(BasicImage<PixelType> & img)
{
    return pair<typename BasicImage<PixelType>::Iterator, 
                typename BasicImage<PixelType>::Accessor>(img.upperLeft(), 
                                         img.accessor());
}

template <class PixelType>
inline pair< typename BasicImage<PixelType>::ConstIterator, 
             typename BasicImage<PixelType>::ConstAccessor>
maskImage(BasicImage<PixelType> const & img)
{
    return pair<typename BasicImage<PixelType>::ConstIterator, 
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
