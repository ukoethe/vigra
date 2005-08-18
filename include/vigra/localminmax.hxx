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
 
 
#ifndef VIGRA_LOCALMINMAX_HXX
#define VIGRA_LOCALMINMAX_HXX

#include <vector>
#include <functional>
#include "vigra/utilities.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/initimage.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/pixelneighborhood.hxx"

namespace vigra {

/** \addtogroup LocalMinMax Local Minima and Maxima

    Detect local minima and maxima of the gray level,
    including extremal plateaus larger than 1 pixel
*/
//@{

/********************************************************/
/*                                                      */
/*                       localMinima                    */
/*                                                      */
/********************************************************/

/** \brief Find local minima in an image.

    The minima are found only when the have a size of one pixel.
    Use \ref extendedLocalMinima() to find minimal plateaus. Minima are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be less-comparable.
    A pixel at the image border will never be marked as minimum. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type>
        void 
        localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da, 
                    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type>
        void 
        localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), minima(w,h);
    
    // init destiniation image
    minima = 0;
    
    vigra::localMinima(srcImageRange(src), destImage(minima));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    int w = slr.x - sul.x - 2;
    int h = slr.y - sul.y - 2;

    static const Diff2D dist[] = {
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1), 
        Diff2D( -1, 0), Diff2D( -1, 1), Diff2D( 0, 1), Diff2D( 1, 1)};
    
    int i,x,y;
    
    sul += Diff2D(1,1);
    dul += Diff2D(1,1);
    
    for(y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator  sx = sul;
        DestIterator dx = dul;
        
        for(x=0; x<w; ++x, ++sx.x, ++dx.x)
        {
            for(i=0; i<8; ++i)
            {
                if(!(sa(sx) < sa(sx, dist[i]))) break;
            }
        
            if(i == 8) da.set(marker, dx);
        }
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    localMinima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second,
                NumericTraits<typename DestAccessor::value_type>::one());
}

/********************************************************/
/*                                                      */
/*                       localMaxima                    */
/*                                                      */
/********************************************************/

/** \brief Find local maxima in an image.

    The maxima are found only when the have a size of one pixel.
    Use \ref extendedLocalMaxima() to find maximal plateaus. Maxima are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be less-comparable.
    A pixel at the image border will never be marked as maximum. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type>
        void 
        localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da, 
                    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type>
        void 
        localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), maxima(w,h);
    
    // init destiniation image
    maxima = 0;
    
    vigra::localMaxima(srcImageRange(src), destImage(maxima));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    int w = slr.x - sul.x - 2;
    int h = slr.y - sul.y - 2;

    static const Diff2D dist[] = {
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1), 
        Diff2D( -1, 0), Diff2D( -1, 1), Diff2D( 0, 1), Diff2D( 1, 1)};
    
    int i,x,y;
    
    sul += Diff2D(1,1);
    dul += Diff2D(1,1);
    
    for(y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator  sx = sul;
        DestIterator dx = dul;
        
        for(x=0; x<w; ++x, ++sx.x, ++dx.x)
        {
            for(i=0; i<8; ++i)
            {
                if(!(sa(sx, dist[i]) < sa(sx))) break;
            }
        
            if(i == 8) da.set(marker, dx);
        }
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    localMaxima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

namespace detail {

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue,
          class Compare, class Equal>
void 
extendedLocalMinmax(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Compare compare, Equal equal)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    int i,x,y;
    
    BasicImage<int> labels(w,h);
    
    int number_of_regions = 
        labelImage(sul, slr, sa, labels.upperLeft(), labels.accessor(), true, equal);

    // assume that a region is a extremum until the opposite is proved
    std::vector<unsigned char> isExtremum(number_of_regions+1, 1);
    
    BasicImage<int>::traverser ly = labels.upperLeft();
    
    for(y=0; y<h; ++y, ++sul.y, ++ly.y)
    {
        SrcIterator  sx = sul;
        BasicImage<int>::traverser lx(ly);
        
        for(x=0; x<w; ++x, ++sx.x, ++lx.x)
        {
            int lab = *lx;
            if(x == 0 || y == 0 || x == w-1 || y == h-1)
            {
                // mark all regions that touch the image border as non-extremum
                isExtremum[lab] = 0;
                continue;
            }
            
            SrcType v = sa(sx);
            NeighborhoodCirculator<SrcIterator, EightNeighborCode> sc(sx);
            NeighborhoodCirculator<BasicImage<int>::traverser, EightNeighborCode> lc(lx);
            for(i=0; i<8; ++i, ++sc, ++lc)
            {
                if(lab != *lc && compare(sa(sc),v))
                    isExtremum[lab] = 0;
            }
            
        }
    }

    ly = labels.upperLeft();
    for(y=0; y<h; ++y, ++dul.y, ++ly.y)
    {
        DestIterator  xd = dul;
        BasicImage<int>::Iterator lx(ly);
        
        for(x=0; x<w; ++x, ++xd.x, ++lx.x)
        {
            if(isExtremum[*lx]) 
                da.set(marker, xd);
        }
    }
}

} // namespace detail

/********************************************************/
/*                                                      */
/*                 extendedLocalMinima                  */
/*                                                      */
/********************************************************/

/** \brief Find local minimal regions in an image.

    This function finds regions of uniform pixel value
    whose neighboring regions are all have smaller values
    (minimal plateaus of arbitrary size). By default, the pixels
    in a plateau have exactly identical values. By passing an <tt>EqualityFunctor</tt>
    with tolerance, one can allow for plateaus that are not quite constant
    (this is often necessary with float pixel values). 
    
    
    Minimal regions are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be equality-comparable and
    less-comparable.
    A pixel or region touching the image border will never be marked as minimum or 
    minimal plateau. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void 
        extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                            DestIterator dul, DestAccessor da, 
                            DestValue marker = NumericTraits<DestValue>::one(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void 
        extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            DestValue marker = NumericTraits<DestValue>::one(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    
    // optional: define an equality functor
    template <class T>
    struct EqualWithToleranceFunctor
    {
        EqualWithToleranceFunctor(T tolerance)
        : t(tolerance)
        {}

        bool operator()(T l, T r) const
        {
            return vigra::abs(l-r) <= t;
        }

        T t;
    };

    vigra::BImage src(w,h), minima(w,h);
    
    // init destiniation image
    minima.init(0);
    
    vigra::extendedLocalMinima(srcImageRange(src), destImage(minima));
    
    // allow plateaus with tolerance
    minima.init(0);
    vigra::extendedLocalMinima(srcImageRange(src), destImage(minima), 1.0,
                               EqualWithToleranceFunctor<unsigned char>(1));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    EqualityFunctor equal;
    u == u  
    equal(u, u);
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue,
          class EqualityFunctor>
void 
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    detail::extendedLocalMinmax(sul, slr, sa, dul, da, 
                                marker, std::less<SrcType>(), equal);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    detail::extendedLocalMinmax(sul, slr, sa, dul, da, 
                                marker, std::less<SrcType>(), std::equal_to<SrcType>());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    extendedLocalMinima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue,
          class EqualityFunctor>
inline void 
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, EqualityFunctor equal)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, equal);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

/********************************************************/
/*                                                      */
/*                 extendedLocalMaxima                  */
/*                                                      */
/********************************************************/

/** \brief Find local maximal regions in an image.

    This function finds regions of uniform pixel value
    whose neighboring regions are all have smaller values
    (maximal plateaus of arbitrary size). By default, the pixels
    in a plateau have exactly identical values. By passing an <tt>EqualityFunctor</tt>
    with tolerance, one can allow for plateaus that are not quite constant
    (this is often necessary with float pixel values). 
    
    
    Maximal regions are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be equality-comparable and
    less-comparable.
    A pixel or region touching the image border will never be marked as maximum or 
    maximal plateau. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void 
        extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                            DestIterator dul, DestAccessor da, 
                            DestValue marker = NumericTraits<DestValue>::one(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor, 
                  class DestValue = DestAccessor::value_type,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void 
        extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            DestValue marker = NumericTraits<DestValue>::one(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    
    // optional: define an equality functor
    template <class T>
    struct EqualWithToleranceFunctor
    {
        EqualWithToleranceFunctor(T tolerance)
        : t(tolerance)
        {}

        bool operator()(T l, T r) const
        {
            return vigra::abs(l-r) <= t;
        }

        T t;
    };

    vigra::BImage src(w,h), maxima(w,h);
    
    // init destiniation image
    maxima.init(0);
    
    vigra::extendedLocalMaxima(srcImageRange(src), destImage(maxima));
    
    // allow plateaus with tolerance
    maxima.init(0);
    vigra::extendedLocalMaxima(srcImageRange(src), destImage(maxima), 1.0,
                               EqualWithToleranceFunctor<unsigned char>(1));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    EqualityFunctor equal;
    u == u  
    equal(u, u);
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue,
          class EqualityFunctor>
void 
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    detail::extendedLocalMinmax(sul, slr, sa, dul, da, 
                                marker, std::greater<SrcType>(), equal);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    detail::extendedLocalMinmax(sul, slr, sa, dul, da, 
                                marker, std::greater<SrcType>(), std::equal_to<SrcType>());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    extendedLocalMaxima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue,
          class EqualityFunctor>
inline void 
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, EqualityFunctor equal)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, equal);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

//@}

} // namespace vigra

#endif // VIGRA_LOCALMINMAX_HXX
