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

/** @heading Image Data Structures and Iterators

    Include-File:
    \URL[vigra/stdimage.hxx]{../include/vigra/stdimage.hxx}
    
    @memo Images, image iterators, and supporting types and functions
*/
//@{

//@Include: basicimage.hxx 

/** @name Standard Image Classes

    @memo Important instantiations of \Ref{BasicImage} template
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
        It uses \Ref{BasicImageIterator} and \Ref{StandardAccessor} and 
        their const counterparts to access the data.
        
        Namespace: vigra
    */
typedef BasicImage<unsigned char> BImage;

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
        It uses \Ref{BasicImageIterator} and \Ref{StandardAccessor} and 
        their const counterparts to access the data.
        
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
        It uses \Ref{BasicImageIterator} and \Ref{StandardAccessor} and 
        their const counterparts to access the data.
        
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
        It uses \Ref{BasicImageIterator} and \Ref{StandardAccessor} and 
    their const counterparts to access the data.
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
        The pixel type is \Ref{RGBValue}<unsigned char>.
        It uses \Ref{BasicImageIterator} and \Ref{RGBAccessor} and 
        their const counterparts to access the data.
        
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
        The pixel type is \Ref{RGBValue}<int>.
        It uses \Ref{BasicImageIterator} and \Ref{RGBAccessor} and 
        their const counterparts to access the data.
        
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
        The pixel type is \Ref{RGBValue}<float>.
        It uses \Ref{BasicImageIterator} and \Ref{RGBAccessor} and 
        their const counterparts to access the data.
        
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

//@Include: imageiterator.hxx 
//@Include: imageiteratoradapter.hxx
//@Include: accessor.hxx
//@Include: rgbvalue.hxx 


/***********************************************************/

/** @name Argument Object Factories

    Long argument list provide for greater flexibility of functions,
    but they are also tedious and error prone, when we don't need
    the flexibility. Thus, we define argument objects which
    automatically provide reasonable defaults for those arguments that we 
    didn't specify explicitly. 
    
    The argument objects are created via a number of factory functions.
    Since these functions have descriptive names, they also serve
    to improve readability: the name of each factory tells te purpose of its
    argument object. 
    
    Consider the following example. Without argument objects we had to 
    write something like this (cf. \Ref{copyImageIf}):
    
    \begin{verbatim}
    vigra::CImage img1, img2, img3;
    
    // fill img1 and img2 ...
    
    vigra::copyImageIf(img1.upperLeft(), img1.lowerRight(), img1.accessor(),
                img2.upperLeft(), img2.accessor(),
                img3.upperLeft(), img3.accessor());
    \end{verbatim}
    
    Using the argument object factories, this becomes much shorter and
    more readable:
    
    \begin{verbatim}
    vigra::copyImageIf(srcImageRange(img1),
                maskImage(img2),
                destImage(img3));
    \end{verbatim}
    
    The names of the factories clearly tell which image is source, mask, 
    and destination. In addition, the suffix #Range# must be used 
    for those argument objects that need to specify the lower right
    corner of the region of interest. Typically, this is only the first
    source argument, but sometimes the first destiniation argument must 
    also contain a range.
    
    The factory functions come in two flavours: Iterator based and 
    image based factories. Above we have seen the image based variant.
    The iterator based variant would look like this:
    
    \begin{verbatim}
    vigra::copyImageIf(srcIterRange(img1.upperLeft(), img1.lowerRight()),
                maskIter(img2.upperLeft()),
                destIter(img3.upperLeft()));
    \end{verbatim}
    
    These factory functions contain the word #Iter# instead of the word 
    #Image#,  They would normally be used if we couldn't access the 
    images (for example, within a function which got passed iterators)
    or if we didn't want to operate on the entire image. The default 
    accessor is obtained via \Ref{IteratorTraits}.
    
    All factory functions also allow to specify accessors explicitly. This
    is useful if we can't use the default accessor. This variant looks 
    like this:
    
    \begin{verbatim}
    vigra::copyImageIf(srcImageRange(img1),
                maskImage(img2, MaskPredicateAccessor()),
                destImage(img3));    
    \end{verbatim}
    
    or
    
    \begin{verbatim}
    vigra::copyImageIf(srcIterRange(img1.upperLeft(), img1.lowerRight()),
                maskIter(img2.upperLeft(), MaskPredicateAccessor()),
                destIter(img3.upperLeft()));
    \end{verbatim}
    
    All versions can be mixed freely within one explession.
    Technically, the argument objects are simply defined as 
    pairs and triples of iterators and accessor so that all algorithms 
    should declare a call interface version based on pairs and triples 
    (see for example \Ref{copyImageIf}).
    
    @memo Factory functions to create argument objects which simplify long argument lists
*/
//@{
    /** @name Image Based Argument Object Factories
        
    These factories can be used to create argument objects when we 
    are given instances or subclasses of \Ref{BasicImage} (see
    \Ref{Standard Image Classes} for instances defined per default).
    These factory functions access #img.upperLeft()#, 
    #img.lowerRight()#, and #img.accessor()# to obtain the iterators
    and accessor for the given image (unless the accessor is 
    given explicitly). The following factory functions are provided:
    
    \begin{tabular}{ll}
    #vigra::BasicImage<SomeType> img;# & 
    
        \\
        
    #srcImageRange(img)# &
        create argument object containing upper left, lower right, and
        default accessor of source image
        
        \\
        
    #srcImageRange(img, SomeAccessor())# &
        create argument object containing upper left, lower right
        of source image, and given accessor
        
        \\
        
    #srcImage(img)# &
        create argument object containing upper left, and
        default accessor of source image
        
        \\
        
    #srcImage(img, SomeAccessor())# &
        create argument object containing upper left
        of source image, and given accessor
        
        \\
    
    #maskImage(img)# &
        create argument object containing upper left, and
        default accessor of mask image
        
        \\
        
    #maskImage(img, SomeAccessor())# &
        create argument object containing upper left
        of mask image, and given accessor
        
        \\
    
    #destImageRange(img)# &
        create argument object containing upper left, lower right, and
        default accessor of destination image
        
        \\
        
    #destImageRange(img, SomeAccessor())# &
        create argument object containing upper left, lower right
        of destination image, and given accessor
        
        \\
        
    #destImage(img)# &
        create argument object containing upper left, and
        default accessor of destination image
        
        \\
        
    #destImage(img, SomeAccessor())# &
        create argument object containing upper left
        of destination image, and given accessor
        
        \\
    
    \end{tabular}
    */
    /** @name Iterator Based Argument Object Factories
        
    These factories can be used to create argument objects when we 
    are given \Ref{Image Iterators}.
    These factory functions use \Ref{IteratorTraits} to
    get the default accessor for the given iterator unless the 
    accessor is given explicitly. The following factory functions 
    are provided:
    
    \begin{tabular}{ll}
    #vigra::BasicImage<SomeType>::Iterator i1, i2;# & 
    
        \\
        
    #srcIterRange(i1, i2)# &
        create argument object containing the given iterators and
        corresponding default accessor (for source image)
        
        \\
        
    #srcIterRange(i1, i2, SomeAccessor())# &
        create argument object containing given iterators and
        accessor (for source image)
        
        \\
        
    #srcIter(i1)# &
        create argument object containing the given iterator and
        corresponding default accessor (for source image)
        
        \\
        
    #srcIter(i1, SomeAccessor())# &
        create argument object containing given iterator and
        accessor (for source image)
        
        \\
    
    #maskIter(i1)# &
        create argument object containing the given iterator and
        corresponding default accessor (for mask image)
        
        \\
        
    #maskIter(i1, SomeAccessor())# &
        create argument object containing given iterator and
        accessor (for mask image)
        
        \\
    
    #destIterRange(i1, i2)# &
        create argument object containing the given iterators and
        corresponding default accessor (for destination image)
        
        \\
        
    #destIterRange(i1, i2, SomeAccessor())# &
        create argument object containing given iterators and
        accessor (for destination image)
        
        \\
        
    #destIter(i1)# &
        create argument object containing the given iterator and
        corresponding default accessor (for destination image)
        
        \\
        
    #destIter(i1, SomeAccessor())# &
        create argument object containing given iterator and
        accessor (for destination image)
        
        \\
    
    \end{tabular}
    */
//@}

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
inline triple<typename BasicImage<PixelType>::Iterator, 
              typename BasicImage<PixelType>::Iterator, Accessor>
srcImageRange(BasicImage<PixelType> & img, Accessor a)
{
    return triple<typename BasicImage<PixelType>::Iterator, 
                  typename BasicImage<PixelType>::Iterator, 
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImage<PixelType>::Iterator, Accessor>
srcImage(BasicImage<PixelType> & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::Iterator, 
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
inline pair<typename BasicImage<PixelType>::Iterator, Accessor>
maskImage(BasicImage<PixelType> & img, Accessor a)
{
    return pair<typename BasicImage<PixelType>::Iterator, 
                Accessor>(img.upperLeft(), a);
}

/****************************************************************/

template <class PixelType>
inline triple<typename BasicImage<PixelType>::Iterator, 
              typename BasicImage<PixelType>::Iterator, 
          typename BasicImage<PixelType>::Accessor>
srcImageRange(BasicImage<PixelType> & img)
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
srcImage(BasicImage<PixelType> & img)
{
    return pair<typename BasicImage<PixelType>::Iterator, 
                typename BasicImage<PixelType>::Accessor>(img.upperLeft(), 
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
inline pair< typename BasicImage<PixelType>::Iterator, 
             typename BasicImage<PixelType>::Accessor>
maskImage(BasicImage<PixelType> & img)
{
    return pair<typename BasicImage<PixelType>::Iterator, 
                typename BasicImage<PixelType>::Accessor>(img.upperLeft(), 
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

//@}

} // namespace vigra

#endif // VIGRA_STDIMAGE_HXX
