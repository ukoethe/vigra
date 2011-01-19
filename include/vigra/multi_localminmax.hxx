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


#ifndef VIGRA_MULTI_LOCALMINMAX_HXX
#define VIGRA_MULTI_LOCALMINMAX_HXX

#include <vector>
#include <functional>
#include "multi_array.hxx"
#include "localminmax.hxx"

namespace vigra {

namespace detail {

// direct neighborhood
template <unsigned int M>
struct IsLocalExtremum2
{
    template <class T, class Shape, class Compare>
    static bool exec(T * v, Shape const & stride, Compare const & compare)
    {
        return compare(*v, *(v-stride[M])) && 
                compare(*v, *(v+stride[M])) && 
                IsLocalExtremum2<M-1>::exec(v, stride, compare);
    }

    template <class Shape, class T, class Compare>
    static bool execAtBorder(Shape const & point, Shape const & shape,
                  T * v, Shape const & stride, Compare const & compare)
    {
        return (point[M] > 0 && compare(*v, *(v-stride[M]))) && 
                (point[M] < shape[M]-1 && compare(*v, *(v+stride[M]))) && 
                IsLocalExtremum2<M-1>::exec(point, shape, v, stride, compare);
    }
};

template <>
struct IsLocalExtremum2<0>
{
    template <class T, class Shape, class Compare>
    static bool exec(T * v, Shape const & stride, Compare const & compare)
    {
        return compare(*v, *(v-stride[0])) &&  compare(*v, *(v+stride[0]));
    }

    template <class Shape, class T, class Compare>
    static bool execAtBorder(Shape const & point, Shape const & shape,
                  T * v, Shape const & stride, Compare const & compare)
    {
        return (point[0] > 0 && compare(*v, *(v-stride[0]))) && 
                (point[0] < shape[0]-1 && compare(*v, *(v+stride[0])));
    }
};

// indirect neighborhood
template <unsigned int M>
struct IsLocalExtremum3
{
    template <class T, class Shape, class Compare>
    static bool exec(T * v, Shape const & stride, Compare const & compare)
    {
        return exec(v, v, stride, compare, true);
    }
    
    template <class T, class Shape, class Compare>
    static bool exec(T * v, T * u, Shape const & stride, 
                     Compare const & compare, bool isCenter)
    {
        return IsLocalExtremum3<M-1>::exec(v, u-stride[M], stride, compare, false) &&
               IsLocalExtremum3<M-1>::exec(v, u, stride, compare, isCenter) &&
               IsLocalExtremum3<M-1>::exec(v, u+stride[M], stride, compare, false);
    }

    template <class Shape, class T, class Compare>
    static bool execAtBorder(Shape const & point, Shape const & shape,
                  T * v, Shape const & stride, Compare const & compare)
    {
        return execAtBorder(point, shape, v, v, stride, compare, true);
    }

    template <class Shape, class T, class Compare>
    static bool execAtBorder(Shape const & point, Shape const & shape,
                             T * v, T * u, Shape const & stride, 
                             Compare const & compare, bool isCenter)
    {
        return (point[M] > 0 && IsLocalExtremum3<M-1>::exec(v, u-stride[M], stride, compare, false)) && 
                IsLocalExtremum3<M-1>::exec(point, shape, v, u, stride, compare, isCenter) &&
                (point[M] < shape[M]-1 && IsLocalExtremum3<M-1>::exec(v, u+stride[M], stride, compare, false));
    }
};

template <>
struct IsLocalExtremum3<0>
{
    template <class T, class Shape, class Compare>
    static bool exec(T * v, Shape const & stride, Compare const & compare)
    {
        return compare(*v, *(v-stride[0])) &&  compare(*v, *(v+stride[0]));
    }

    template <class T, class Shape, class Compare>
    static bool exec(T * v, T * u, Shape const & stride, 
                     Compare const & compare, bool isCenter)
    {
        return compare(*v, *(u-stride[0])) &&  
                (!isCenter && compare(*v, *u)) &&
                compare(*v, *(u+stride[0]));
    }

    template <class Shape, class T, class Compare>
    static bool execAtBorder(Shape const & point, Shape const & shape,
                  T * v, Shape const & stride, Compare const & compare)
    {
        return (point[0] > 0 && compare(*v, *(v-stride[0]))) && 
                (point[0] < shape[0]-1 && compare(*v, *(v+stride[0])));
    }

    template <class Shape, class T, class Compare>
    static bool execAtBorder(Shape const & point, Shape const & shape,
                             T * v, T * u, Shape const & stride, 
                             Compare const & compare, bool isCenter)
    {
        return (point[M] > 0 && compare(*v, *(u-stride[0]))) && 
                (!isCenter && compare(*v, *u)) &&
                (point[M] < shape[M]-1 && compare(*v, *(u+stride[0])));
    }
};

template <unsigned int N, class T1, class C1, class T2, class C2, class Compare>
void
localMinMax(MultiArrayView<N, T1, C1> src,
            MultiArrayView<N, T2, C2> dest,
            T2 marker, unsigned int neighborhood,
            T1 threshold,
            Compare compare,
            bool allowExtremaAtBorder = false)
{
    typedef typename MultiArrayShape<N>::type Shape;
    typedef MultiCoordinateNavigator<N> Navigator;
    
    Shape shape = src.shape(),
          unit  = Shape(MultiArrayIndex(1));
    
    vigra_precondition(shape == dest.shape(),
        "localMinMax(): Shape mismatch between input and output.");
    vigra_precondition(neighborhood == 2*N || neighborhood == pow(3, N) - 1,
        "localMinMax(): Invalid neighborhood.");
    
    if(allowExtremaAtBorder)
    {
        for(unsigned int d=0; d<N; ++d)
        {
            Navigator nav(shape, d);
            for(; nav.hasMore(); ++nav)
            {
                Shape i = nav.begin();
                
                for(; i[d] < shape[d]; i[d] += shape[d]-1)
                {                    
                    if(!compare(src[i], threshold))
                        continue;
                    
                    if(neighborhood == 2*N)
                    {
                        if(IsLocalExtremum2<N>::execAtBorder(i, shape, &src[i], 
                                                              src.stride(), compare))
                            dest[i] = marker;
                    }
                    else
                    {
                        if(IsLocalExtremum3<N>::execAtBorder(i, shape, &src[i], 
                                                              src.stride(), compare))
                            dest[i] = marker;
                    }
                }
            }
        }
    }

    src = src.subarray(unit, shape - unit);
    dest = dest.subarray(unit, shape - unit);
    shape = src.shape();

    Navigator nav(shape, 0);
    for(; nav.hasMore(); ++nav)
    {
        Shape i = nav.begin();
        
        for(; i[0] < shape[0]; ++i[0])
        {                    
            if(!compare(src[i], threshold))
                continue;
            
            if(neighborhood == 2*N)
            {
                if(IsLocalExtremum2<N>::exec(&src[i], src.stride(), compare))
                    dest[i] = marker;
            }
            else
            {
                if(IsLocalExtremum3<N>::exec(&src[i], src.stride(), compare))
                    dest[i] = marker;
            }
        }
    }
}

template <class T1, class C1, class T2, class C2, 
          class Neighborhood, class Compare, class Equal>
void
extendLocalMinMax(MultiArrayView<3, T1, C1> src,
                  MultiArrayView<3, T2, C2> dest,
                  T2 marker, 
                  Neighborhood neighborhood,
                  Compare compare, Equal equal, 
                  T1 threshold,
                  bool allowExtremaAtBorder = false)
{
    typedef typename MultiArrayView<3, T1, C1>::traverser SrcIterator;
    typedef typename MultiArrayView<3, T2, C2>::traverser DestIterator;
    typedef typename MultiArray<3, int>::traverser LabelIterator;
    typedef MultiArrayShape<3>::type Shape;
    
    Shape shape = src.shape();
    MultiArrayIndex w = shape[0], h = shape[1], d = shape[2];
    
    vigra_precondition(shape == dest.shape(),
        "extendLocalMinMax(): Shape mismatch between input and output.");

    MultiArray<3, int> labels(shape);

    int number_of_regions = labelVolume(srcMultiArrayRange(src), destMultiArray(labels),
                                        neighborhood, equal);

    // assume that a region is a extremum until the opposite is proved
    ArrayVector<unsigned char> isExtremum(number_of_regions+1, (unsigned char)1);
    
    SrcIterator zs = src.traverser_begin();
    LabelIterator zl = labels.traverser_begin();
    for(MultiArrayIndex z = 0; z != d; ++z, ++zs.dim2(), ++zl.dim2())
    {
        SrcIterator ys(zs);
        LabelIterator yl(zl);

        for(MultiArrayIndex y = 0; y != h; ++y, ++ys.dim1(), ++yl.dim1())
        {
            SrcIterator xs(ys);
            LabelIterator xl(yl);

            for(MultiArrayIndex x = 0; x != w; ++x, ++xs.dim0(), ++xl.dim0())
            {
                int lab = *xl;
                T1 v = *xs;
                
                if(isExtremum[lab] == 0)
                    continue;
                    
                if(!compare(v, threshold))
                {
                    // mark all regions that don't exceed the threshold as non-extremum
                    isExtremum[lab] = 0;
                    continue;
                }

                AtVolumeBorder atBorder = isAtVolumeBorder(x, y, z, w, h, d);
                if(atBorder == NotAtBorder)
                {
                    NeighborhoodCirculator<SrcIterator, Neighborhood> cs(xs);
                    NeighborhoodCirculator<LabelIterator, Neighborhood> cl(xl);
                    for(i=0; i<Neighborhood::DirectionCount; ++i, ++cs, ++cl)
                    {
                        if(lab != *cl && compare(*cs,v))
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
                                                                   cs(xs, atBorder), scend(cs);
                        do
                        {
                            if(lab != *(xl+cs.diff()) && compare(cs,v))
                            {
                                isExtremum[lab] = 0;
                                break;
                            }
                        }
                        while(++cs != scend);
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
    DestIterator zd = dest.traverser_begin();
    for(MultiArrayIndex z = 0; z != d; ++z, ++zl.dim2(), ++zd.dim2())
    {
        LabelIterator yl(zl);
        DestIterator yd(zd);

        for(MultiArrayIndex y = 0; y != h; ++y, ++yl.dim1(), ++yd.dim1())
        {
            LabelIterator xl(yl);
            DestIterator xd(yd);

            for(MultiArrayIndex x = 0; x != w; ++x, ++xl.dim0(), ++xd.dim0())
            {
                if(isExtremum[*xl])
                    *xd = marker;
            }
        }
    }
}

} // namespace detail

/********************************************************/
/*                                                      */
/*                       localMinima                    */
/*                                                      */
/********************************************************/

// documentation is in localminmax.hxx
template <unsigned int N, class T1, class C1, class T2, class C2>
void
localMinima(MultiArrayView<N, T1, C1> src,
            MultiArrayView<N, T2, C2> dest,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    T1 threshold = options.use_threshold
                           ? std::min(NumericTraits<T1>::max(), (T1)options.thresh)
                           : NumericTraits<T1>::max();
    T2 marker = (T2)options.marker;
    
    vigra_precondition(!options.allow_plateaus,
        "localMinima(): Option 'allowPlateaus' is not implemented for arbitrary dimensions,\n"
        "               use extendedLocalMinima() for 2D and 3D problems.");

    if(options.neigh == 0)
        options.neigh = 2*N;
    if(options.neigh == 1)
        options.neigh = pow(3, N) - 1;
    
    detail::localMinMax(src, dest, marker, options.neigh, 
                        threshold, std::less<T1>(), options.allow_at_border);
}

/********************************************************/
/*                                                      */
/*                       localMaxima                    */
/*                                                      */
/********************************************************/

// documentation is in localminmax.hxx
template <unsigned int N, class T1, class C1, class T2, class C2>
void
localMaxima(MultiArrayView<N, T1, C1> src,
            MultiArrayView<N, T2, C2> dest,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    T1 threshold = options.use_threshold
                           ? std::max(NumericTraits<T1>::min(), (T1)options.thresh)
                           : NumericTraits<T1>::min();
    T2 marker = (T2)options.marker;
    
    vigra_precondition(!options.allow_plateaus,
        "localMaxima(): Option 'allowPlateaus' is not implemented for arbitrary dimensions,\n"
        "               use extendedLocalMinima() for 2D and 3D problems.");

    if(options.neigh == 0)
        options.neigh = 2*N;
    if(options.neigh == 1)
        options.neigh = pow(3, N) - 1;
    
    detail::localMinMax(src, dest, marker, options.neigh, 
                        threshold, std::greater<T1>(), options.allow_at_border);
}

/**************************************************************************/

/********************************************************/
/*                                                      */
/*                 extendedLocalMinima                  */
/*                                                      */
/********************************************************/

// documentation is in localminmax.hxx
template <class T1, class C1, class T2, class C2,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMinima(MultiArrayView<3, T1, C1> src,
                    MultiArrayView<3, T2, C2> dest,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    T1 threshold = options.use_threshold
                           ? std::min(NumericTraits<T1>::max(), (T1)options.thresh)
                           : NumericTraits<T1>::max();
    T2 marker = (T2)options.marker;
    
    if(options.neigh == 0 || options.neigh == 6)
    {
        detail::extendedLocalMinMax(src, dest, marker, NeighborCode3DSix(), 
                                    threshold, std::less<T1>(), std::equal_to<T1>(), 
                                    options.allow_at_border);
    }
    else if(options.neigh == 1 || options.neigh == 26)
    {
        detail::extendedLocalMinMax(src, dest, marker, NeighborCode3DTwentySix(), 
                                    threshold, std::less<T1>(), std::equal_to<T1>(), 
                                    options.allow_at_border);
    }
    else
        vigra_precondition(false,
          "extendedLocalMinima(): Invalid neighborhood.");
}

/********************************************************/
/*                                                      */
/*                 extendedLocalMaxima                  */
/*                                                      */
/********************************************************/

// documentation is in localminmax.hxx
template <class T1, class C1, class T2, class C2,
          class Neighborhood, class EqualityFunctor>
inline void
extendedLocalMaxima(MultiArrayView<3, T1, C1> src,
                    MultiArrayView<3, T2, C2> dest,
                    LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    T1 threshold = options.use_threshold
                           ? std::max(NumericTraits<T1>::min(), (T1)options.thresh)
                           : NumericTraits<T1>::min();
    T2 marker = (T2)options.marker;
    
    if(options.neigh == 0 || options.neigh == 6)
    {
        detail::extendedLocalMinMax(src, dest, marker, NeighborCode3DSix(), 
                                    threshold, std::greater<T1>(), std::equal_to<T1>(), 
                                    options.allow_at_border);
    }
    else if(options.neigh == 1 || options.neigh == 26)
    {
        detail::extendedLocalMinMax(src, dest, marker, NeighborCode3DTwentySix(), 
                                    threshold, std::greater<T1>(), std::equal_to<T1>(), 
                                    options.allow_at_border);
    }
    else
        vigra_precondition(false,
          "extendedLocalMaxima(): Invalid neighborhood.");
}

} // namespace vigra

#endif // VIGRA_MULTI_LOCALMINMAX_HXX
