/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
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


#ifndef VIGRA_LOCALMINMAX_HXX
#define VIGRA_LOCALMINMAX_HXX

#include <vector>
#include <functional>
#include "utilities.hxx"
#include "stdimage.hxx"
#include "initimage.hxx"
#include "labelimage.hxx"
#include "pixelneighborhood.hxx"

namespace vigra {

/** \addtogroup LocalMinMax Local Minima and Maxima

    Detect local minima and maxima of the gray level,
    including extremal plateaus larger than 1 pixel
*/
//@{

namespace detail {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood,
          class Compare>
void
localMinMax(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                DestIterator dul, DestAccessor da,
                DestValue marker, Neighborhood,
                typename SrcAccessor::value_type threshold,
                Compare compare)
{
    int w = slr.x - sul.x - 2;
    int h = slr.y - sul.y - 2;

    int i,x,y;

    sul += Diff2D(1,1);
    dul += Diff2D(1,1);

    for(y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator  sx = sul;
        DestIterator dx = dul;

        for(x=0; x<w; ++x, ++sx.x, ++dx.x)
        {
			if(!compare(sa(sx), threshold)) continue;
            typename SrcAccessor::value_type v = sa(sx);
            NeighborhoodCirculator<SrcIterator, Neighborhood> sc(sx);
            for(i = 0; i < Neighborhood::DirectionCount; ++i, ++sc)
            {
                if(!compare(v, sa(sc)))
                    break;
            }

            if(i == Neighborhood::DirectionCount)
                da.set(marker, dx);
        }
    }
}

} // namespace detail


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
    Pass \ref vigra::EightNeighborCode or \ref vigra::FourNeighborCode
    to determine the neighborhood where pixel values are compared.
    If a threshold is given, only pixels with values less than 
    that threshold are marked.
    The function uses accessors.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class DestValue = DestAccessor::value_type,
                  class Neighborhood = EightNeighborCode>
        void
        localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da,
                    DestValue marker = NumericTraits<DestValue>::one(),
                    Neighborhood neighborhood = EightNeighborCode(),
                    SrcValue threshold = NumericTraits<SrcValue>::max())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class DestValue = DestAccessor::value_type,
                  class Neighborhood = EightNeighborCode>
        void
        localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    DestValue marker = NumericTraits<DestValue>::one(),
                    Neighborhood neighborhood = EightNeighborCode(),
                    SrcValue threshold = NumericTraits<SrcValue>::max())
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void localMinima)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, Neighborhood neighborhood,
            typename SrcAccessor::value_type threshold)
{
    detail::localMinMax(sul, slr, sa, dul, da, marker, neighborhood,
					threshold,
                    std::less<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, Neighborhood neighborhood)
{
    localMinima(sul, slr, sa, dul, da, marker, neighborhood,
				NumericTraits<typename SrcAccessor::value_type>::max());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker)
{
    localMinima(sul, slr, sa, dul, da, marker, EightNeighborCode(),
				NumericTraits<typename SrcAccessor::value_type>::max());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    localMinima(sul, slr, sa, dul, da,
                NumericTraits<typename DestAccessor::value_type>::one(),
                EightNeighborCode(),
                NumericTraits<typename SrcAccessor::value_type>::max());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood, 
            typename SrcAccessor::value_type threshold)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, threshold);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood,
                NumericTraits<typename SrcAccessor::value_type>::max());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, EightNeighborCode(), 
                NumericTraits<typename SrcAccessor::value_type>::max());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second,
                NumericTraits<typename DestAccessor::value_type>::one(),
                EightNeighborCode(), 
                NumericTraits<typename SrcAccessor::value_type>::max());
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
    If a threshold is given, only pixels with values greater than 
    that threshold are marked.
    The function uses accessors.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class DestValue = DestAccessor::value_type,
                  class Neighborhood = EightNeighborCode>
        void
        localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da,
                    DestValue marker = NumericTraits<DestValue>::one(),
                    Neighborhood neighborhood = EightNeighborCode(),
                    SrcValue threshold = NumericTraits<SrcValue>::min())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class DestValue = DestAccessor::value_type,
                  class Neighborhood = EightNeighborCode>
        void
        localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    DestValue marker = NumericTraits<DestValue>::one(),
                    Neighborhood neighborhood = EightNeighborCode(),
                    SrcValue threshold = NumericTraits<SrcValue>::min())
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void localMaxima)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, Neighborhood neighborhood,
            typename SrcAccessor::value_type threshold)
{
    detail::localMinMax(sul, slr, sa, dul, da, marker, neighborhood, threshold,
                    std::greater<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class DestValue, class Neighborhood>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, Neighborhood neighborhood)
{
    localMaxima(sul, slr, sa, dul, da, marker, neighborhood, 
				NumericTraits<typename SrcAccessor::value_type>::min());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker)
{
    localMaxima(sul, slr, sa, dul, da, marker, EightNeighborCode(),
				NumericTraits<typename SrcAccessor::value_type>::min());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    localMaxima(sul, slr, sa, dul, da,
                NumericTraits<typename DestAccessor::value_type>::one(),
                EightNeighborCode(), 
                NumericTraits<typename SrcAccessor::value_type>::min());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood,
            typename SrcAccessor::value_type threshold)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, threshold);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, 
                NumericTraits<typename SrcAccessor::value_type>::min());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, EightNeighborCode(),
                NumericTraits<typename SrcAccessor::value_type>::min());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second,
                NumericTraits<typename DestAccessor::value_type>::one(),
                EightNeighborCode(), 
                NumericTraits<typename SrcAccessor::value_type>::min());
}

namespace detail {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood, class Compare, class Equal>
void
extendedLocalMinMax(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Neighborhood /*neighborhood*/, Compare compare, Equal equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    int i,x,y;

    BasicImage<int> labels(w,h);

    int number_of_regions =
        labelImage(sul, slr, sa, labels.upperLeft(), labels.accessor(),
                   (Neighborhood::DirectionCount == 8), equal);

    // assume that a region is a extremum until the opposite is proved
    std::vector<unsigned char> isExtremum(number_of_regions+1, (unsigned char)1);

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
            NeighborhoodCirculator<SrcIterator, Neighborhood> sc(sx);
            NeighborhoodCirculator<BasicImage<int>::traverser, Neighborhood> lc(lx);
            for(i=0; i<Neighborhood::DirectionCount; ++i, ++sc, ++lc)
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
    (this is often necessary with float pixel values). Pass
    \ref vigra::EightNeighborCode or \ref vigra::FourNeighborCode
    to determine the neighborhood where pixel values are compared.


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
                  class Neighborhood = EightNeighborCode,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void
        extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                            DestIterator dul, DestAccessor da,
                            DestValue marker = NumericTraits<DestValue>::one(),
                            Neighborhood neighborhood = EightNeighborCode(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class DestValue = DestAccessor::value_type,
                  class Neighborhood = EightNeighborCode,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void
        extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            DestValue marker = NumericTraits<DestValue>::one(),
                            Neighborhood neighborhood = EightNeighborCode(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void extendedLocalMinima)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Neighborhood neighborhood, EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    detail::extendedLocalMinMax(sul, slr, sa, dul, da,
                                marker, neighborhood,
                                std::less<SrcType>(), equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood>
inline void
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Neighborhood neighborhood)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMinima(sul, slr, sa, dul, da,
                        marker, neighborhood, std::equal_to<SrcType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMinima(sul, slr, sa, dul, da,
                        marker, EightNeighborCode());
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
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood,
            EqualityFunctor equal)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood>
inline void
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood)
{
    extendedLocalMinima(src.first, src.second, src.third,
                        dest.first, dest.second, marker, neighborhood);
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
                        dest.first, dest.second);
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
    (this is often necessary with float pixel values). Pass
    \ref vigra::EightNeighborCode or \ref vigra::FourNeighborCode
    to determine the neighborhood where pixel values are compared.


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
                  class Neighborhood = EightNeighborCode,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void
        extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                            DestIterator dul, DestAccessor da,
                            DestValue marker = NumericTraits<DestValue>::one(),
                            Neighborhood neighborhood = EightNeighborCode(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class DestValue = DestAccessor::value_type,
                  class Neighborhood = EightNeighborCode,
                  class EqualityFunctor = std::equal_to<typename SrcAssessor::value_type> >
        void
        extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            DestValue marker = NumericTraits<DestValue>::one(),
                            Neighborhood neighborhood = EightNeighborCode(),
                            EqualityFunctor equal = EqualityFunctor())
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void extendedLocalMaxima)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Neighborhood neighborhood, EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    detail::extendedLocalMinMax(sul, slr, sa, dul, da,
                                marker, neighborhood,
                                std::greater<SrcType>(), equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood>
inline void
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Neighborhood neighborhood)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMaxima(sul, slr, sa, dul, da,
                        marker, neighborhood, std::equal_to<SrcType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMaxima(sul, slr, sa, dul, da,
                        marker, EightNeighborCode());
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
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood,
            EqualityFunctor equal)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood>
inline void
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, Neighborhood neighborhood)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                        dest.first, dest.second, marker, neighborhood);
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
                        dest.first, dest.second);
}

//@}

} // namespace vigra

#endif // VIGRA_LOCALMINMAX_HXX
