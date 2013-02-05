/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2010 by Ullrich Koethe                  */
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
#include "labelvolume.hxx"
#include "pixelneighborhood.hxx"
#include "voxelneighborhood.hxx"

namespace vigra
{

/** \addtogroup LocalMinMax Local Minima and Maxima

    Detect local minima and maxima in a gray level image,
    including extremal plateaus larger than 1 pixel
*/
//@{

namespace detail {

template <class SrcIterator, class SrcAccessor,
          class Neighborhood,
          class Compare>
inline bool
isLocalExtremum(SrcIterator is, SrcAccessor sa, Neighborhood,
                typename SrcAccessor::value_type threshold,
                Compare compare, AtImageBorder atBorder)
{
    typename SrcAccessor::value_type v = sa(is);
    
    if(!compare(v, threshold))
        return false;

    int directionCount = Neighborhood::nearBorderDirectionCount(atBorder);
    RestrictedNeighborhoodCirculator<SrcIterator, Neighborhood> sc(is, atBorder);
    for(int i = 0; i < directionCount; ++i, ++sc)
    {
        if(!compare(v, sa(sc)))
            return false;
    }
    return true;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue, class Neighborhood,
          class Compare>
void
localMinMax(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, Neighborhood neighborhood,
            typename SrcAccessor::value_type threshold,
            Compare compare,
            bool allowExtremaAtBorder = false)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    int x, y;

    if(allowExtremaAtBorder)
    {
        SrcIterator is = sul;
        DestIterator id = dul;
        
        for(x=0; x<w; ++x, ++is.x, ++id.x)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare, 
                                isAtImageBorder(x, 0, w, h)))
                da.set(marker, id);
        }
        
        is = sul + Diff2D(0,1);
        id = dul + Diff2D(0,1);
        
        for(y=1; y<h-1; ++y, ++is.y, ++id.y)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare, 
                                isAtImageBorder(0, y, w, h)))
                da.set(marker, id);
        }
        
        is = sul + Diff2D(w-1,1);
        id = dul + Diff2D(w-1,1);
        
        for(y=1; y<h-1; ++y, ++is.y, ++id.y)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare, 
                                isAtImageBorder(w-1, y, w, h)))
                da.set(marker, id);
        }
        
        is = sul + Diff2D(0,h-1);
        id = dul + Diff2D(0,h-1);
        
        for(x=0; x<w; ++x, ++is.x, ++id.x)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare, 
                                isAtImageBorder(x, h-1, w, h)))
                da.set(marker, id);
        }
    }

    w -= 2;
    h -= 2;
    sul += Diff2D(1,1);
    dul += Diff2D(1,1);

    for(y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator  sx = sul;
        DestIterator dx = dul;

        for(x=0; x<w; ++x, ++sx.x, ++dx.x)
        {
            typename SrcAccessor::value_type v = sa(sx);
            
            if(!compare(v, threshold))
                continue;

            int i;
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

template<class SrcIterator, class SrcShape, class SrcAccessor,
        class DestIterator, class DestAccessor, class DestValue,
        class Neighborhood, class Compare>
void 
localMinMax3D(SrcIterator sul, SrcShape shp, SrcAccessor sa, 
              DestIterator dul, DestAccessor da,
              DestValue marker,
              Neighborhood neighborhood,
              typename SrcAccessor::value_type threshold,
              Compare compare,
              bool allowExtremaAtBorder = false)
{
    int w = shp[0];
    int h = shp[1];
    int d = shp[2];

    int x, y, z;

    if (allowExtremaAtBorder)
    {
        throw std::runtime_error("not implemented!");
        /*
        SrcIterator is = sul;
        DestIterator id = dul;

        for(x=0; x<w; ++x, ++is.x, ++id.x)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare,
                        isAtImageBorder(x, 0, w, h)))
            da.set(marker, id);
        }

        is = sul + Diff2D(0,1);
        id = dul + Diff2D(0,1);

        for(y=1; y<h-1; ++y, ++is.y, ++id.y)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare,
                        isAtImageBorder(0, y, w, h)))
            da.set(marker, id);
        }

        is = sul + Diff2D(w-1,1);
        id = dul + Diff2D(w-1,1);

        for(y=1; y<h-1; ++y, ++is.y, ++id.y)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare,
                        isAtImageBorder(w-1, y, w, h)))
            da.set(marker, id);
        }

        is = sul + Diff2D(0,h-1);
        id = dul + Diff2D(0,h-1);

        for(x=0; x<w; ++x, ++is.x, ++id.x)
        {
            if(isLocalExtremum(is, sa, neighborhood, threshold, compare,
                        isAtImageBorder(x, h-1, w, h)))
            da.set(marker, id);
        }
    */
    }

    w -= 2;
    h -= 2;
    d -= 2;
    sul.dim0() += 1;
    sul.dim1() += 1;
    sul.dim2() += 1;
    dul += Diff3D(1, 1, 1);

    SrcIterator zs = sul;
    DestIterator zd = dul;

    for (z = 0; z != d; ++z, ++zs.dim2(), ++zd.dim2())
    {
        SrcIterator ys(zs);
        DestIterator yd(zd);

        for (y = 0; y != h; ++y, ++ys.dim1(), ++yd.dim1())
        {
            SrcIterator xs(ys);
            DestIterator xd(yd);

            for (x = 0; x != w; ++x, ++xs.dim0(), ++xd.dim0())
            {

                typename SrcAccessor::value_type v = sa(xs);
                if (!compare(v, threshold))
                    continue;

                int i;
                NeighborhoodCirculator<SrcIterator, Neighborhood> sc(xs);
                for (i = 0; i < Neighborhood::DirectionCount; ++i, ++sc)
                {
                    if(!compare(v, sa(sc)))
                        break;
                }

                if(i == Neighborhood::DirectionCount)
                    da.set(marker, xd);
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood, class Compare, class Equal>
void
extendedLocalMinMax(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da, DestValue marker,
                    Neighborhood /*neighborhood*/,
                    Compare compare, Equal equal, 
                    typename SrcAccessor::value_type threshold,
                    bool allowExtremaAtBorder = false)
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
            SrcType v = sa(sx);
            
            if(isExtremum[lab] == 0)
                continue;
                
            if(!compare(v, threshold))
            {
                // mark all regions that don't exceed the threshold as non-extremum
                isExtremum[lab] = 0;
                continue;
            }

            AtImageBorder atBorder = isAtImageBorder(x, y, w, h);
            if(atBorder == NotAtBorder)
            {
                NeighborhoodCirculator<SrcIterator, Neighborhood> sc(sx);
                NeighborhoodCirculator<BasicImage<int>::traverser, Neighborhood> lc(lx);
                for(i=0; i<Neighborhood::DirectionCount; ++i, ++sc, ++lc)
                {
                    if(lab != *lc && compare(sa(sc),v))
                    {
                        isExtremum[lab] = 0;
                        break;
                    }
                }
            }
            else
            {
                if(allowExtremaAtBorder)
                {
                    RestrictedNeighborhoodCirculator<SrcIterator, Neighborhood> 
                                                               sc(sx, atBorder), scend(sc);
                    do
                    {
                        if(lab != *(lx+sc.diff()) && compare(sa(sc),v))
                        {
                            isExtremum[lab] = 0;
                            break;
                        }
                    }
                    while(++sc != scend);
                }
                else
                {
                    isExtremum[lab] = 0;
                }
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

template<class SrcIterator, class SrcShape, class SrcAccessor,
        class DestIterator, class DestAccessor, class DestValue,
        class Neighborhood, class Compare, class Equal>
void 
extendedLocalMinMax3D(SrcIterator sul, SrcShape shp, SrcAccessor sa,
                      DestIterator dul, DestAccessor da,
                      DestValue marker,
                      Neighborhood neighbourhood,
                      Compare compare,
                      Equal equal,
                      typename SrcAccessor::value_type threshold,
                      bool allowExtremaAtBorder = false)
{
    typedef typename SrcAccessor::value_type SrcType;

    int w = shp[0];
    int h = shp[1];
    int d = shp[2];

    int i, x, y, z;

    MultiArray<3, int> labels(shp);

    int number_of_regions =
        labelVolume(sul, shp, sa, labels.traverser_begin(), 
                    typename AccessorTraits<int>::default_accessor(),
                    neighbourhood);
    
    MultiArray<3, int>::traverser zl(labels.traverser_begin());

    SrcIterator zs = sul;
    DestIterator zd = dul;

    // assume that a region is a extremum until the opposite is proved
    std::vector<unsigned char> isExtremum(number_of_regions + 1, (unsigned char)1);

    for (z = 0; z != d; ++z, ++zs.dim2(), ++zd.dim2(), ++zl.dim2())
    {
        SrcIterator ys(zs);
        DestIterator yd(zd);
        MultiArray<3, int>::traverser yl(zl);

        for (y = 0; y != h; ++y, ++ys.dim1(), ++yd.dim1(), ++yl.dim1())
        {
            SrcIterator xs(ys);
            DestIterator xd(yd);
            MultiArray<3, int>::traverser xl(yl);

            for (x = 0; x != w; ++x, ++xs.dim0(), ++xd.dim0(), ++xl.dim0())
            {

                int lab = *xl;
                SrcType v = sa(xs);

                if (isExtremum[lab] == 0)
                    continue;

                if (!compare(v, threshold))
                {
                    // mark all regions that don't exceed the threshold as non-extremum
                    isExtremum[lab] = 0;
                    continue;
                }

                AtVolumeBorder atBorder = isAtVolumeBorder(x, y, z, w, h, d);
                if (atBorder == NotAtBorder)
                {
                    NeighborhoodCirculator<SrcIterator, Neighborhood> sc(xs);
                    NeighborhoodCirculator<MultiArray<3, int>::traverser, Neighborhood> lc(xl);
                    for (i = 0; i < Neighborhood::DirectionCount; ++i, ++sc, ++lc)
                    {
                        if (lab != *lc && compare(sa(sc), v))
                        {

                            isExtremum[lab] = 0;
                            break;
                        }
                    }
                }
                else
                {
                    if (allowExtremaAtBorder)
                    {
                        RestrictedNeighborhoodCirculator<SrcIterator, Neighborhood>
                        sc(xs, atBorder), scend(sc);
                        do
                        {
                            if (lab != *(xl + sc.diff()) && compare(sa(sc), v))
                            {
                                isExtremum[lab] = 0;
                                break;
                            }
                        } 
                        while (++sc != scend);
                    }
                    else
                    {
                        isExtremum[lab] = 0;
                    }
                }
            }
        }
    }

    zl = labels.traverser_begin();
    zs = sul;
    zd = dul;

    for (z = 0; z != d; ++z, ++zs.dim2(), ++zd.dim2(), ++zl.dim2())
    {
        SrcIterator ys(zs);
        DestIterator yd(zd);
        MultiArray<3, int>::traverser yl(zl);

        for (y = 0; y != h; ++y, ++ys.dim1(), ++yd.dim1(), ++yl.dim1())
        {
            SrcIterator xs(ys);
            DestIterator xd(yd);
            MultiArray<3, int>::traverser xl(yl);

            for (x = 0; x != w; ++x, ++xs.dim0(), ++xd.dim0(), ++xl.dim0())
            {
                if(isExtremum[*xl])
                    da.set(marker, xd);
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue,
          class Neighborhood, class Compare, class Equal>
void
extendedLocalMinMaxOld(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker,
            Neighborhood /*neighborhood*/,
            Compare compare, Equal equal, 
            typename SrcAccessor::value_type threshold,
            bool allowExtremaAtBorder = false)
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
            SrcType v = sa(sx);
            if(x == 0 || y == 0 || x == w-1 || y == h-1 || !compare(v, threshold))
            {
                // mark all regions that touch the image border as non-extremum
                // likewise for all pixels that don't exceed the threshold
                isExtremum[lab] = 0;
                continue;
            }

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


/** \brief Options object for localMinima() and localMaxima().

    <b> Usage:</b>

    <b>\#include</b> \<vigra/localminmax.hxx\><br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), minima(w,h);
    ... // fill src

    // use 4-neighborhood, allow minima at the image border, 
    // and discard those where the gray value is not below 5
    vigra::localMinima(srcImageRange(src), destImage(minima),
                       vigra::LocalMinmaxOptions().neighborhood(4).allowAtBorder().threshold(5));

    \endcode
*/
class LocalMinmaxOptions
{
  public:
    double marker, thresh;
    int neigh;
    bool use_threshold, allow_at_border, allow_plateaus;
    
        /**\brief Construct default options object.
         *
            Defaults are: marker value '1', no threshold, indirect neighborhood, 
                          don't allow extrema at border and extremal plateaus.
         */
    LocalMinmaxOptions()
    : marker(1.0), 
      thresh(0.0),
      neigh(1),
      use_threshold(false),
      allow_at_border(false),
      allow_plateaus(false)
    {}
    
        /**\brief Use the given neighborhood. 
        
            The value '0' indicates direct neighborhood (i.e. 4-neighborhood 
            in 2D, 6-neighborhood in 3D, 2*N neighborhood in N-D), the value '1'
            indicates indirect neighborhood (i.e. 8-neighborhood in 2D, 
            26-neighborhood in 3D, 3<sup>N</sup>-1 neighborhood in N-D). The specific number 
            of neighbors for the desired dimension can also be used.
        
            Default: 1 (indirect neighborhood)
         */
    LocalMinmaxOptions & neighborhood(unsigned int n)
    {
        neigh = n;
        return *this;
    }
    
        /**\brief Mark extrema in the destination image with the given value.
        
            Default: 1
         */
    LocalMinmaxOptions & markWith(double m)
    {
        marker = m;
        return *this;
    }
    
        /**\brief Threshold the extrema.
        
           Discard minima whose gray value is not below the threshold.
           and maxima whose gray level is not above the threshold.
        
            Default: don't threshold (i.e. return all extrema)
         */
    LocalMinmaxOptions & threshold(double t)
    {
        use_threshold = true;
        thresh = t;
        return *this;
    }
    
        /**\brief Detect extrema at the image border.
        
            Default: false
         */
    LocalMinmaxOptions & allowAtBorder(bool f = true)
    {
        allow_at_border = f;
        return *this;
    }
    
        /**\brief Allow extremal plateaus.
        
            That is regions of constant gray value whose neighbors are all
            higher (minima) or lower than the value of the region.
        
            Default: false
         */
    LocalMinmaxOptions & allowPlateaus(bool f = true)
    {
        allow_plateaus = f;
        return *this;
    }
};


/********************************************************/
/*                                                      */
/*                       localMinima                    */
/*                                                      */
/********************************************************/

/** \brief Find local minima in an image or multi-dimensional array.

    Note: the function is not yet implemented for arbitrary dimensional
    arrays, but see \ref localMinima3D() for 3D.
    
    By default, minima are defined as points which are not 
    at the array border and whose value is lower than the value 
    of all indirect neighbors (i.e. 8-neighbors in 2D, 
    26-neighbors in 3D, 3<sup>N</sup>-1 neighbors in N-D). 
    The detected points will be marked 
    with the default value 1 in the destination array.
    
    The defaults can be overridden in various ways by providing 
    \ref LocalMinmaxOptions : you can switch to the direct neighborhood
    (i.e. 4-neighborhood in 2D, 6-neighborhood in 3D, 2*N neighborhood 
    in N-D), allow minima at the border, discard minima where the function 
    value is not below a given threshold, allow extended minima
    (i.e. minima that form minimal plateaus rather than isolated pixels --
    note that this option is only supported for 2D images), 
    and change the marker in the destination image. See usage examples below 
    for details. 
    
    There are also variants of the localMinima() function where parameters
    are passed explicitly rather than via an option object. These versions
    of the function are deprecated, but will be kept for compatibility.

    <b> Declarations:</b>

    use arbitrary-dimensional arrays:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class C1, class T2, class C2>
        void
        localMinima(MultiArrayView<N, T1, C1> src,
                    MultiArrayView<N, T2, C2> dest,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions());
    }
    \endcode

    pass image iterators explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions());
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions());
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/localminmax.hxx\><br>
    <b>\#include</b> \<vigra/multi_localminmax.hxx\><br>
    Namespace: vigra

    \code
    // 3D examples using MultiArray
    vigra::MultiArrayShape<3>::type shape(w,h,d);
    vigra::MultiArray<3, unsigned char> src(shape), minima(shape);
    ... // fill src

    // use default parameterisation
    vigra::localMinima(src, minima);

    // reset destination image
    minima = 0;

    // use 6-neighborhood and allow minima at the image border
    vigra::localMinima(src, minima,
                       vigra::LocalMinmaxOptions().neighborhood(6).allowAtBorder());
    \endcode

    \code
    // 2D examples using BasicImage
    vigra::BImage src(w,h), minima(w,h);
    ... // fill src

    // use default parameterisation
    vigra::localMinima(srcImageRange(src), destImage(minima));

    // reset destination image
    minima = 0;

    // use 4-neighborhood and allow minima at the image border
    vigra::localMinima(srcImageRange(src), destImage(minima),
                       vigra::LocalMinmaxOptions().neighborhood(4).allowAtBorder());

    // reset destination image
    minima = 0;

    // allow extended minima (minimal plateaus) and use value '255' as a marker
    vigra::localMinima(srcImageRange(src), destImage(minima),
                       vigra::LocalMinmaxOptions().allowPlateaus().markWith(255));
    \endcode

    <b> Required Interface:</b>

    \code
    SrcIterator src_upperleft, src_lowerright;
    DestIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u < u
    \endcode
*/
doxygen_overloaded_function(template <...> void localMinima)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    
    SrcType threshold = options.use_threshold
                           ? std::min(NumericTraits<SrcType>::max(), (SrcType)options.thresh)
                           : NumericTraits<SrcType>::max();
    DestType marker = (DestType)options.marker;
    
    if(options.allow_plateaus)
    {
        if(options.neigh == 0 || options.neigh == 4)
        {
            detail::extendedLocalMinMax(sul, slr, sa, dul, da, marker, FourNeighborCode(),
                                        std::less<SrcType>(), std::equal_to<SrcType>(), 
                                        threshold, options.allow_at_border);
        }
        else if(options.neigh == 1 || options.neigh == 8)
        {
            detail::extendedLocalMinMax(sul, slr, sa, dul, da, marker, EightNeighborCode(),
                                        std::less<SrcType>(), std::equal_to<SrcType>(), 
                                        threshold, options.allow_at_border);
        }
        else
            vigra_precondition(false, "localMinima(): neighborhood must be 4 or 8.");

    }
    else
    {
        if(options.neigh == 0 || options.neigh == 4)
        {
            detail::localMinMax(sul, slr, sa, dul, da, marker, FourNeighborCode(),
                                threshold, std::less<SrcType>(), options.allow_at_border);
        }
        else if(options.neigh == 1 || options.neigh == 8)
        {
            detail::localMinMax(sul, slr, sa, dul, da, marker, EightNeighborCode(),
                                threshold, std::less<SrcType>(), options.allow_at_border);
        }
        else
            vigra_precondition(false, "localMinima(): neighborhood must be 4 or 8.");
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, FourNeighborCode neighborhood)
{
    detail::localMinMax(sul, slr, sa, dul, da, marker, neighborhood,
                        NumericTraits<typename SrcAccessor::value_type>::max(),
                        std::less<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, EightNeighborCode neighborhood)
{
    detail::localMinMax(sul, slr, sa, dul, da, marker, neighborhood,
                        NumericTraits<typename SrcAccessor::value_type>::max(),
                        std::less<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker)
{
    localMinima(sul, slr, sa, dul, da, marker, EightNeighborCode());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, FourNeighborCode neighborhood)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, EightNeighborCode neighborhood)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, EightNeighborCode());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, options);
}

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                       localMinima3D                  */
/*                                                      */
/********************************************************/

/** \brief Find local minima in a 3D multi array.

 By default, minima are defined as points which are not
 at the array border and whose value is lower than the value
 of all indirect neighbors.
 The detected points will be marked. See localMinima() for more details.

 */
doxygen_overloaded_function(template <...> void localMinima3D)

template<class SrcIterator, class SrcAccessor, class SrcShape,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMinima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
              DestIterator dul, DestAccessor da,
              DestValue marker,
              NeighborCode3DTwentySix neighborhood)
{
    detail::localMinMax3D(sul, slr, sa, dul, da, marker, neighborhood,
                NumericTraits<typename SrcAccessor::value_type>::max(),
                std::less<typename SrcAccessor::value_type>());
}

template<class SrcIterator, class SrcAccessor, class SrcShape,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMinima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
              DestIterator dul, DestAccessor da,
              DestValue marker,
              NeighborCode3DSix neighborhood)
{
    detail::localMinMax3D(sul, slr, sa, dul, da, marker, neighborhood,
                NumericTraits<typename SrcAccessor::value_type>::max(),
                std::less<typename SrcAccessor::value_type>());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMinima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
              DestIterator dul, DestAccessor da,
              DestValue marker)
{
    localMinima3D(sul, slr, sa, dul, da, marker, NeighborCode3DSix());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMinima3D(triple<SrcIterator, SrcShape, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              DestValue marker,
              NeighborCode3DSix neighborhood)
{
    localMinima3D(src.first, src.second, src.third, dest.first, dest.second,
                marker, neighborhood);
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void
localMinima3D(triple<SrcIterator, SrcShape, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              DestValue marker,
              NeighborCode3DTwentySix neighborhood)
{
    localMinima3D(src.first, src.second, src.third, dest.first, dest.second,
                marker, neighborhood);
}

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                       localMaxima                    */
/*                                                      */
/********************************************************/

/** \brief Find local maxima in an image or multi-dimensional array.

    Note: the function is not yet implemented for arbitrary dimensional
    arrays, but see \ref localMaxima3D() for 3D.
    
    By default, maxima are defined as points which are not 
    at the array border and whose value is higher than the value 
    of all indirect neighbors (i.e. 8-neighbors in 2D, 
    26-neighbors in 3D, 3<sup>N</sup>-1 neighbors in N-D). 
    The detected points will be marked 
    with the default value 1 in the destination array.
    
    The defaults can be overridden in various ways by providing 
    \ref LocalMinmaxOptions : you can switch to the direct neighborhood
    (i.e. 4-neighborhood in 2D, 6-neighborhood in 3D, 2*N neighborhood 
    in N-D), allow maxima at the border, discard maxima where the function 
    value is not above a given threshold, allow extended maxima
    (i.e. maxima that form maximal plateaus rather than isolated pixels --
    note that this option is only supported for 2D images), 
    and change the marker in the destination image. See usage examples below 
    for details. 
    
    There are also variants of the localMaxima() function where parameters
    are passed explicitly rather than via an option object. These versions
    of the function are deprecated, but will be kept for compatibility.

    <b> Declarations:</b>

    use arbitrary-dimensional arrays:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class C1, class T2, class C2>
        void
        localMaxima(MultiArrayView<N, T1, C1> src,
                    MultiArrayView<N, T2, C2> dest,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions());
    }
    \endcode

    pass image iterators explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                    DestIterator dul, DestAccessor da,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions());
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions());
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<vigra/localminmax.hxx\><br>
        <b>\#include</b> \<vigra/multi_localminmax.hxx\><br>
    Namespace: vigra

    \code
    // 3D examples using MultiArray
    vigra::MultiArrayShape<3>::type shape(w,h,d);
    vigra::MultiArray<3, unsigned char> src(shape), maxima(shape);
    ... // fill src

    // use default parameterisation
    vigra::localMaxima(src, maxima);

    // reset destination image
    maxima = 0;

    // use 6-neighborhood and allow maxima at the image border
    vigra::localMaxima(src, maxima,
                       vigra::LocalMinmaxOptions().neighborhood(6).allowAtBorder());
    \endcode

    \code
    // 2D examples using BasicImage
    vigra::BImage src(w,h), maxima(w,h);
    ... // fill src

    // use default parameterisation
    vigra::localMaxima(srcImageRange(src), destImage(maxima));

    // reset destination image
    maxima = 0;

    // use 4-neighborhood and allow maxima at the image border
    vigra::localMaxima(srcImageRange(src), destImage(maxima),
                       vigra::LocalMinmaxOptions().neighborhood(4).allowAtBorder());

    // reset destination image
    maxima = 0;

    // allow extended maxima (maximal plateaus) and use value '255' as a marker
    vigra::localMaxima(srcImageRange(src), destImage(maxima),
                       vigra::LocalMinmaxOptions().allowPlateaus().markWith(255));
    \endcode

    <b> Required Interface:</b>

    \code
    SrcIterator src_upperleft, src_lowerright;
    DestIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u < u
    \endcode
*/
doxygen_overloaded_function(template <...> void localMaxima)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    
    SrcType threshold = options.use_threshold
                           ? std::max(NumericTraits<SrcType>::min(), (SrcType)options.thresh)
                           : NumericTraits<SrcType>::min();
    DestType marker = (DestType)options.marker;
    
    if(options.allow_plateaus)
    {
        if(options.neigh == 0 || options.neigh == 4)
        {
            detail::extendedLocalMinMax(sul, slr, sa, dul, da, marker, FourNeighborCode(),
                                        std::greater<SrcType>(), std::equal_to<SrcType>(), 
                                        threshold, options.allow_at_border);
        }
        else if(options.neigh == 1 || options.neigh == 8)
        {
            detail::extendedLocalMinMax(sul, slr, sa, dul, da, marker, EightNeighborCode(),
                                        std::greater<SrcType>(), std::equal_to<SrcType>(), 
                                        threshold, options.allow_at_border);
        }
        else
            vigra_precondition(false, "localMaxima(): neighborhood must be 4 or 8.");
    }
    else
    {
        if(options.neigh == 0 || options.neigh == 4)
        {
            detail::localMinMax(sul, slr, sa, dul, da, marker, FourNeighborCode(),
                                threshold, std::greater<SrcType>(), options.allow_at_border);
        }
        else if(options.neigh == 1 || options.neigh == 8)
        {
            detail::localMinMax(sul, slr, sa, dul, da, marker, EightNeighborCode(),
                                threshold, std::greater<SrcType>(), options.allow_at_border);
        }
        else
            vigra_precondition(false, "localMaxima(): neighborhood must be 4 or 8.");
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, FourNeighborCode neighborhood)
{
    detail::localMinMax(sul, slr, sa, dul, da, marker, neighborhood,
                        NumericTraits<typename SrcAccessor::value_type>::min(),
                        std::greater<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker, EightNeighborCode neighborhood)
{
    detail::localMinMax(sul, slr, sa, dul, da, marker, neighborhood,
                        NumericTraits<typename SrcAccessor::value_type>::min(),
                        std::greater<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da,
            DestValue marker)
{
    localMaxima(sul, slr, sa, dul, da, marker, EightNeighborCode());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, FourNeighborCode neighborhood)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class DestValue>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker, EightNeighborCode neighborhood)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DestValue marker)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, EightNeighborCode());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, options);
}

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                       localMaxima3D                  */
/*                                                      */
/********************************************************/

/** \brief Find local maxima in a 3D multi array.

 By default, maxima are defined as points which are not
 at the array border and whose value is higher than the value
 of all indirect neighbors.
 The detected points will be marked as specified. See localMaxima() for mor details.
 */
doxygen_overloaded_function(template <...> void localMaxima3D)

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
              DestIterator dul, DestAccessor da,
              DestValue marker,
              NeighborCode3DSix neighborhood)
{
    detail::localMinMax3D(sul, slr, sa, dul, da, marker, neighborhood,
                NumericTraits<typename SrcAccessor::value_type>::min(),
                std::greater<typename SrcAccessor::value_type>());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
              DestIterator dul, DestAccessor da,
              DestValue marker,
              NeighborCode3DTwentySix neighborhood)
{
    detail::localMinMax3D(sul, slr, sa, dul, da, marker, neighborhood,
                NumericTraits<typename SrcAccessor::value_type>::min(),
                std::greater<typename SrcAccessor::value_type>());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima3D(triple<SrcIterator, SrcShape, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              DestValue marker,
              NeighborCode3DTwentySix neighborhood)
{
    localMaxima3D(src.first, src.second, src.third, dest.first, dest.second,
                marker, neighborhood);
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima3D(vigra::triple<SrcIterator, SrcShape, SrcAccessor> src,
              std::pair<DestIterator, DestAccessor> dest,
              DestValue marker)
{
    localMaxima3D(src.first, src.second, src.third, dest.first, dest.second,
                marker, NeighborCode3DSix());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima3D(triple<SrcIterator, SrcShape, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              DestValue marker,
              NeighborCode3DSix neighborhood)
{
    localMaxima3D(src.first, src.second, src.third, dest.first, dest.second,
                marker, neighborhood);
}

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                 extendedLocalMinima                  */
/*                                                      */
/********************************************************/

/** \brief Find local minimal regions in an image or volume.

    Note: the function is not yet implemented for arbitrary dimensional
    arrays, but see \ref extendedLocalMinima3D() for 3D.
    
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
    less-comparable. A pixel or region touching the image border will 
    never be marked as minimum or minimal plateau. Use localMinima() with the 
    appropriate options if you need that functionality. Likewise if you want to
    apply a threshold onl the fly. In fact, all functionality
    except for 'equality with tolerance' can be accessed via that function in
    a more readable way, so localMinima() should be preferred.
    The function uses accessors.

    <b> Declarations:</b>

    use 3-dimensional arrays:
    \code
    namespace vigra {
        template <class T1, class C1, class T2, class C2,
                  class Neighborhood>
        void
        extendedLocalMinima(MultiArrayView<3, T1, C1> src,
                            MultiArrayView<3, T2, C2> dest,
                            LocalMinmaxOptions const & options = LocalMinmaxOptions());
    \endcode

    pass image iterators explicitly:
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
                            EqualityFunctor equal = EqualityFunctor());
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
                            EqualityFunctor equal = EqualityFunctor());
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<vigra/localminmax.hxx\><br>
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
                               FourNeighborCode(),
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
          class DestIterator, class DestAccessor,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, 
            typename DestAccessor::value_type marker,
            Neighborhood neighborhood, EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    detail::extendedLocalMinMax(sul, slr, sa, dul, da,
                                marker, neighborhood,
                                std::less<SrcType>(), equal, 
                                NumericTraits<typename SrcAccessor::value_type>::max());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline void
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, 
            typename DestAccessor::value_type marker,
            Neighborhood neighborhood)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMinima(sul, slr, sa, dul, da,
                        marker, neighborhood, std::equal_to<SrcType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, 
            typename DestAccessor::value_type marker)
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
          class DestIterator, class DestAccessor,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            typename DestAccessor::value_type marker, Neighborhood neighborhood,
            EqualityFunctor equal)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline void
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            typename DestAccessor::value_type marker, Neighborhood neighborhood)
{
    extendedLocalMinima(src.first, src.second, src.third,
                        dest.first, dest.second, marker, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            typename DestAccessor::value_type marker)
{
    extendedLocalMinima(src.first, src.second, src.third,
                        dest.first, dest.second, marker, EightNeighborCode());
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

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                 extendedLocalMinima3D                */
/*                                                      */
/********************************************************/

/** \brief Find local minimal regions in a volume.

 This function finds regions of uniform pixel value
 whose neighboring regions are all have smaller values
 (minimal plateaus of arbitrary size). By default, the pixels
 in a plateau have exactly identical values. By passing an <tt>EqualityFunctor</tt>
 with tolerance, one can allow for plateaus that are not quite constant
 (this is often necessary with float pixel values). Pass the neighborhood
 where pixel values are compared. See extendedLocalMinima() for more details.

*/
doxygen_overloaded_function(template <...> void extendedLocalMinima3D)

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class Neighborhood,
         class EqualityFunctor>
inline void 
extendedLocalMinima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
                      DestIterator dul, DestAccessor da,
                      typename DestAccessor::value_type marker,
                      Neighborhood neighborhood,
                      EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    detail::extendedLocalMinMax3D(sul, slr, sa, dul, da, marker, neighborhood,
                std::less<SrcType>(), equal,
                NumericTraits<typename SrcAccessor::value_type>::max());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class Neighborhood>
inline void 
extendedLocalMinima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
                      DestIterator dul, DestAccessor da,
                      typename DestAccessor::value_type marker,
                      Neighborhood neighborhood)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMinima3D(sul, slr, sa, dul, da, marker, neighborhood,
                          std::equal_to<SrcType>());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor>
inline void 
extendedLocalMinima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
                      DestIterator dul, DestAccessor da)
{
    extendedLocalMinima3D(sul, slr, sa, dul, da, 
                          NumericTraits<typename DestAccessor::value_type>::one(),
                          NeighborCode3DSix());
}

template<class SrcIterator, class SrcAccessor, class SrcShape,
         class DestIterator, class DestAccessor, class Neighborhood>
inline void 
extendedLocalMinima3D(triple<SrcIterator, SrcShape, SrcAccessor> src,
                      pair<DestIterator, DestAccessor> dest,
                      typename DestAccessor::value_type marker,
                      Neighborhood neighborhood)
{
    extendedLocalMinima3D(src.first, src.second, src.third, 
                          dest.first, dest.second, 
                          marker, neighborhood);
}

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                 extendedLocalMaxima                  */
/*                                                      */
/********************************************************/

/** \brief Find local maximal regions in an image or volume.

    Note: the function is not yet implemented for arbitrary dimensional
    arrays, but see \ref extendedLocalMaxima3D() for 3D.
    
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
    less-comparable. A pixel or region touching the image border will 
    never be marked as maximum or maximal plateau. Use localMaxima() with the 
    appropriate options if you need that functionality. Likewise if you want to
    apply a threshold onl the fly. In fact, all functionality
    except for 'equality with tolerance' can be accessed via that function in
    a more readable way, so localMaxima() should be preferred.
    The function uses accessors.

    <b> Declarations:</b>

    use 3-dimensional arrays:
    \code
    namespace vigra {
        template <class T1, class C1, class T2, class C2,
                  class Neighborhood>
        void
        extendedLocalMaxima(MultiArrayView<3, T1, C1> src,
                            MultiArrayView<3, T2, C2> dest,
                            LocalMinmaxOptions const & options = LocalMinmaxOptions());
    \endcode

    pass image iterators explicitly:
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

        <b>\#include</b> \<vigra/localminmax.hxx\><br>
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
                               FourNeighborCode(),
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
          class DestIterator, class DestAccessor,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, 
            typename DestAccessor::value_type marker,
            Neighborhood neighborhood, EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    detail::extendedLocalMinMax(sul, slr, sa, dul, da,
                                marker, neighborhood,
                                std::greater<SrcType>(), equal, 
                                NumericTraits<typename SrcAccessor::value_type>::min());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline void
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, 
            typename DestAccessor::value_type marker,
            Neighborhood neighborhood)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMaxima(sul, slr, sa, dul, da,
                        marker, neighborhood, std::equal_to<SrcType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, 
            typename DestAccessor::value_type marker)
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
          class DestIterator, class DestAccessor,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            typename DestAccessor::value_type marker, Neighborhood neighborhood,
            EqualityFunctor equal)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker, neighborhood, equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline void
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            typename DestAccessor::value_type marker, Neighborhood neighborhood)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                        dest.first, dest.second, marker, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            typename DestAccessor::value_type marker)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                        dest.first, dest.second, marker, EightNeighborCode());
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

/********************************************************/
/*                                                      */
/*                 extendedLocalMaxima3D                */
/*                                                      */
/********************************************************/

/** \brief Find local maximal regions in 3D multi array.

 This function finds regions of uniform pixel value
 whose neighboring regions are all have smaller values
 (maximal plateaus of arbitrary size). By default, the pixels
 in a plateau have exactly identical values. By passing an <tt>EqualityFunctor</tt>
 with tolerance, one can allow for plateaus that are not quite constant
 (this is often necessary with float pixel values). Pass
 the neighborhood where pixel values are compared. See extendedLocalMaxima() for more details.
 */

doxygen_overloaded_function(template <...> void extendedLocalMaxima3D)

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class Neighborhood,
         class EqualityFunctor>
inline void 
extendedLocalMaxima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
                      DestIterator dul, DestAccessor da,
                      typename DestAccessor::value_type marker,
                      Neighborhood neighborhood,
                      EqualityFunctor equal)
{
    typedef typename SrcAccessor::value_type SrcType;

    detail::extendedLocalMinMax3D(sul, slr, sa, dul, da, marker, neighborhood,
                                  std::greater<SrcType>(), equal, 
                                  NumericTraits<typename SrcAccessor::value_type>::min());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class Neighborhood>
inline void 
extendedLocalMaxima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
                      DestIterator dul, DestAccessor da,
                      typename DestAccessor::value_type marker,
                      Neighborhood neighborhood)
{
    typedef typename SrcAccessor::value_type SrcType;

    extendedLocalMaxima3D(sul, slr, sa, dul, da, 
                          marker, neighborhood, 
                          std::equal_to<SrcType>());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor>
inline void 
extendedLocalMaxima3D(SrcIterator sul, SrcShape slr, SrcAccessor sa,
                      DestIterator dul, DestAccessor da)
{
    extendedLocalMaxima3D(sul, slr, sa, dul, da, 
                          NumericTraits<typename DestAccessor::value_type>::one(),
                          NeighborCode3DSix());
}

template<class SrcIterator, class SrcShape, class SrcAccessor,
         class DestIterator, class DestAccessor, class Neighborhood>
inline void 
extendedLocalMaxima3D(triple<SrcIterator, SrcShape, SrcAccessor> src,
                      pair<DestIterator, DestAccessor> dest,
                      typename DestAccessor::value_type marker,
                      Neighborhood neighborhood)
{
    extendedLocalMaxima3D(src.first, src.second, src.third, 
                          dest.first, dest.second, 
                          marker, neighborhood);
}

//@}

} // namespace vigra

#endif // VIGRA_LOCALMINMAX_HXX
