/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_RESIZEIMAGE_HXX
#define VIGRA_RESIZEIMAGE_HXX

#include <vector>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/imageiteratoradapter.hxx"
#include "vigra/recursiveconvolution.hxx"

namespace vigra {

/** \addtogroup GeometricTransformations Geometric Transformations
    Zoom up and down by repeating pixels, or using linear or spline interpolation
    
    <b>\#include</b> "<a href="stdimagefunctions_8hxx-source.html">vigra/stdimagefunctions.hxx</a>"<br>   
    <b>or</b><br>
    <b>\#include</b> "<a href="resizeimage_8hxx-source.html">vigra/resizeimage.hxx</a>"<br>   
*/
//@{

/********************************************************/
/*                                                      */
/*                   resizeBandLimitLine                */
/*         (internally used by resize functions)        */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor, 
          class TempVector,
          class DestIterator, class DestAccessor>
void resizeBandLimitLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                         TempVector & line,
                         DestIterator id, DestAccessor ad, double width)
{
    typedef typename TempVector::value_type TempType;
    
    int w = isend - is;
    
    int x;
    
    double b = exp(-1.0/width);
    double norm = (1.0 - b) / (1.0 + b);
    TempType old = (1.0 / (1.0 - b)) * as(is);

    // left side of filter
    for(x=0; x<w; ++x, ++is)
    {
        old = as(is) + b * old;
        line[x] = old;
    }
    
    // right side of the filter
    --is;
    old = (1.0 / (1.0 - b)) * as(is);
    id += w;
    ++is;
    
    for(x=w-1; x>=0; --x)
    {    
        --is;
        --id;
        
        old = b * old;
        ad.set(norm * (line[x] + old), id);

        old = as(is) + old;
    }
}

/********************************************************/
/*                                                      */
/*           resizeCalculateSplineCoefficients          */
/*         (internally used by resize functions)        */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor, class VALUETYPE>
void 
resizeCalculateSplineCoefficients(SrcIterator i1, SrcIterator iend, 
                SrcAccessor a, VALUETYPE * i2)
{
    int n = iend - i1;
    
    if(n <= 0) return;
    
    VALUETYPE zero = NumericTraits<VALUETYPE>::zero();
    VALUETYPE two = 2.0 * NumericTraits<VALUETYPE>::one();
    VALUETYPE half = 0.5 * NumericTraits<VALUETYPE>::one();
    
    *i2 = zero;
    if(n == 1) return;

    std::vector<VALUETYPE> vec(n);
    typename std::vector<VALUETYPE>::iterator u = vec.begin();
    
    *u = zero;
    
    for(++i1, ++i2, ++u, --iend; i1 != iend; ++i1, ++i2, ++u)
    {
        VALUETYPE p = 0.5 * i2[-1] + two;
        *i2 = half / p;
        *u = 3.0 *(a(i1,1) - 2.0 * a(i1) + a(i1, -1)) - 0.5 * u[-1] / p;
    }
    
    *i2 = zero;
    
    for(--i2, --u; u != vec; --u, --i2)
    {
        *i2 = *i2 * i2[1] + *u;
    }
}

/********************************************************/
/*                                                      */
/*         resizeImageInternalSplineGradient            */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor, 
          class DoubleIterator, class TempIterator, class DestIterator>
void 
resizeImageInternalSplineGradient(SrcIterator in, SrcIterator inend, SrcAccessor sa,
                         DoubleIterator tmp, TempIterator r, DestIterator id)
{
    int w = inend - in;
    
    int x;
    
    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;

    // calculate border derivatives
    SrcIterator xs = in;
    TMPTYPE p0 = -11.0/6.0 * sa(xs);  ++xs;
            p0 += 3.0 * sa(xs);  ++xs;
            p0 += -1.5 * sa(xs);  ++xs;
            p0 += 1.0/3.0 * sa(xs);
        
    xs = in + w-1;
    TMPTYPE pw = 11.0/6.0 * sa(xs);  --xs;
            pw += -3.0 * sa(xs);  --xs;
            pw +=  1.5 * sa(xs);  --xs;
            pw += -1.0/3.0 * sa(xs);
        
    xs = in + 2;
    SrcIterator xs1 = in;
        
    for(x=1; x<w-1; ++x, ++xs, ++xs1)
    {
        r[x] = 3.0 * (sa(xs) - sa(xs1));
    }
    
    r[1] -= p0;
    r[w-2] -= pw;

    double q = 0.25;
    
    id[0] = p0;
    id[w-1] = pw;
    id[1] = 0.25 * r[1];

    for(x=2; x<w-1; ++x)
    {
        tmp[x] = q;
        q = 1.0 / (4.0 - q);
        id[x] = q * (r[x] - id[x-1]);
    }
    
    for(x=w-3; x>=1; --x)
    {
        id[x] -= tmp[x+1]*id[x+1];
    }
}

/********************************************************/
/*                                                      */
/*         resizeImageInternalSplineInterpolation       */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
resizeImageInternalSplineInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    int w = iend.x - is.x;
    int h = iend.y - is.y;
    
    int wnew = idend.x - id.x;
    int hnew = idend.y - id.y;
        
    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
    typedef typename BasicImage<TMPTYPE>::Iterator TMPITER;
    typedef 
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    BasicImage<TMPTYPE> dx(w,h);
    BasicImage<TMPTYPE> dy(w,h);
    BasicImage<TMPTYPE> dxy(w,h);
    BasicImage<TMPTYPE> W(4,4), W1(4,4);
    std::vector<TMPTYPE> R(w > h ? w : h);
    std::vector<double> tmp(w > h ? w : h);
    
    typename BasicImage<TMPTYPE>::Accessor ta;
    
    SrcIterator in = is;
    
    TMPITER idx = dx.upperLeft();
    TMPITER idy = dy.upperLeft();
    TMPITER idxy = dxy.upperLeft();
    typename std::vector<TMPTYPE>::iterator r = R.begin();
    typename std::vector<double>::iterator it = tmp.begin();
    
    double ig[] = { 1.0, 0.0, -3.0,  2.0,
                    0.0, 1.0, -2.0,  1.0,
                    0.0, 0.0,  3.0, -2.0,
                    0.0, 0.0, -1.0,  1.0 };
    
    int x, y, i, j, k;

    
    // calculate x derivatives
    for(y=0; y<h; ++y, ++in.y, ++idx.y)
    {
        RowIterator<SrcIterator> sr = in;
        RowIterator<TMPITER> dr = idx;
        resizeImageInternalSplineGradient(sr, sr+w, accessorAdapter(sr, sa),
                                          it, r, dr);
    }
    
    in = is;
        
    // calculate y derivatives
    for(x=0; x<w; ++x, ++in.x, ++idy.x)
    {
        ColumnIterator<SrcIterator> sc = in;
        ColumnIterator<TMPITER> dc = idy;
        resizeImageInternalSplineGradient(sc, sc+h, accessorAdapter(sc, sa),  
                                          it, r, dc);
    }
        
    in = is;
    idy = dy.upperLeft();

    // calculate mixed derivatives
    for(y=0; y<h; ++y, ++idy.y, ++idxy.y)
    {
        RowIterator<TMPITER> sr = idy;
        RowIterator<TMPITER> dr = idxy;
        resizeImageInternalSplineGradient(sr, sr+w, accessorAdapter(sr, ta), 
                                          it, r, dr);
    }
    
    double du = (double)(w-1) / (wnew-1);
    double dv = (double)(h-1) / (hnew-1);
    double ov = 0.0;
    int oy = 0;
    int yy = oy;
    
    DestIterator xxd = id, yyd = id;

    static Diff2D down(0,1), right(1,0), downright(1,1);
   
    for(y=0; y<h-1; ++y, ++in.y, ov -= 1.0)
    {
        if(y < h-2 && ov >= 1.0) continue;
        int y1 = y+1;
        double v = ov;
        double ou = 0.0;
        int ox = 0;
        int xx = ox;
        
        SrcIterator xs = in;
        for(x=0; x<w-1; ++x, ++xs.x, ou -= 1.0)
        {
            if(x < w-2 && ou >= 1.0) continue;
            int x1 = x+1;
            double u = ou;
            
            DestIterator xd = id + Diff2D(ox,oy);
            W[0][0] = sa(xs);
            W[0][1] = dy(x, y);
            W[0][2] = sa(xs, down);
            W[0][3] = dy(x, y1);
            W[1][0] = dx(x, y);
            W[1][1] = dxy(x, y);
            W[1][2] = dx(x, y1);
            W[1][3] = dxy(x, y1);
            W[2][0] = sa(xs, right);
            W[2][1] = dy(x1,y);
            W[2][2] = sa(xs, downright);
            W[2][3] = dy(x1, y1);
            W[3][0] = dx(x1, y);
            W[3][1] = dxy(x1, y);
            W[3][2] = dx(x1, y1);
            W[3][3] = dxy(x1, y1);
            
            for(i=0; i<4; ++i)
            {
                for(j=0; j<4; ++j)
                {
                    W1[j][i] = ig[j] * W[0][i];
                    for(k=1; k<4; ++k)
                    {
                        W1[j][i] += ig[j+4*k] * W[k][i];
                    }
                }
            }
            for(i=0; i<4; ++i)
            {
                for(j=0; j<4; ++j)
                {
                    W[j][i] = ig[i] * W1[j][0];
                    for(k=1; k<4; ++k)
                    {
                       W[j][i] += ig[4*k+i] * W1[j][k];
                    }
                }
            }
            
            TMPTYPE a1,a2,a3,a4;
            
            yyd = xd;
            for(v=ov, yy=oy; v<1.0; v+=dv, ++yyd.y, ++yy)
            {
                a1 = W[0][0] + v * (W[0][1] + 
                               v * (W[0][2] + v * W[0][3]));
                a2 = W[1][0] + v * (W[1][1] + 
                               v * (W[1][2] + v * W[1][3]));
                a3 = W[2][0] + v * (W[2][1] + 
                               v * (W[2][2] + v * W[2][3]));
                a4 = W[3][0] + v * (W[3][1] + 
                               v * (W[3][2] + v * W[3][3]));
                
                xxd = yyd;
                for(u=ou, xx=ox; u<1.0; u+=du, ++xxd.x, ++xx)
                {
                    da.set(DestTraits::fromRealPromote(a1 + u * (a2 + u * (a3 + u * a4))), xxd);
                }

                if(xx == wnew-1)
                {
                    da.set(DestTraits::fromRealPromote(a1 + a2 + a3 + a4), xxd);
                }
            }
            
            if(yy == hnew-1)
            {
                a1 = W[0][0] + W[0][1] + W[0][2] + W[0][3];
                a2 = W[1][0] + W[1][1] + W[1][2] + W[1][3];
                a3 = W[2][0] + W[2][1] + W[2][2] + W[2][3];
                a4 = W[3][0] + W[3][1] + W[3][2] + W[3][3];

                DestIterator xxd = yyd;
                for(u=ou, xx=ox; u<1.0; u+=du, ++xxd.x, ++xx)
                {
                    da.set(DestTraits::fromRealPromote(a1 + u * (a2 + u * (a3 + u * a4))), xxd);
                }
                
                if(xx == wnew-1)
                {
                    da.set(DestTraits::fromRealPromote(a1 + a2 + a3 + a4), xxd);
                }
            }
            
            ou = u;
            ox = xx;
        }
        ov = v;
        oy = yy;
    }
}
    
/********************************************************/
/*                                                      */
/*           resizeImageSplineInterpolation             */
/*                                                      */
/********************************************************/

/** \brief Resize image using bi-cubic spline interpolation.

    The function uses the bi-cubic, non-separable spline algorithm described in 
    [Hoschek/Lasser: 
    {\it "Grundlagen der geometrischen Datenverarbeitung"}, Teubner, 1992] to obtain 
    optimal interpolation quality.
    
    The range of both the input and output images (resp. regions)
    must be given. The input image must have a size of at
    least 4x4, the destination of at least 2x2. The scaling factors are then calculated
    accordingly. If the source image is larger than the destination, it
    is smoothed (band limited) using a recursive
    exponential filter. The source value_type (SrcAccessor::value_type) must 
    be a linear algebra, i.e. it must support addition, subtraction,
    and multiplication (+, -, *), multiplication with a scalar 
    real number and \ref NumericTraits "NumericTraits".
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void 
        resizeImageSplineInterpolation(
              SrcImageIterator is, SrcImageIterator iend, SrcAccessor sa,
          DestImageIterator id, DestImageIterator idend, DestAccessor da)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void 
        resizeImageSplineInterpolation(
              triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
          triple<DestImageIterator, DestImageIterator, DestAccessor> dest)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="resizeimage_8hxx-source.html">vigra/resizeimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::resizeImageSplineInterpolation(
               src.upperLeft(), src.lowerRight(), src.accessor(), 
               dest.upperLeft(), dest.lowerRight(), dest.accessor());
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft, src_lowerright;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote 
                             u = src_accessor(src_upperleft),
                 v = src_accessor(src_upperleft, 1);
    double d;
    
    u = d * v;
    u = u + v;
    u = u - v;
    u = u * v;
    u += v;
    u -= v;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(u), 
    dest_upperleft);

    \endcode

    <b> Preconditions:</b>
    
    \code
    src_lowerright.x - src_upperleft.x > 3
    src_lowerright.y - src_upperleft.y > 3
    dest_lowerright.x - dest_upperleft.x > 1
    dest_lowerright.y - dest_upperleft.y > 1
    \endcode
    
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
resizeImageSplineInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    int w = iend.x - is.x;
    int h = iend.y - is.y;
    
    int wnew = idend.x - id.x;
    int hnew = idend.y - id.y;
    
    vigra_precondition((w > 3) && (h > 3),
                 "resizeImageSplineInterpolation(): "
         "Source image to small.\n");
    vigra_precondition((wnew > 1) && (hnew > 1),
                 "resizeImageSplineInterpolation(): "
         "Destination image to small.\n");
    
    double scale = 2.0;
    
    if(wnew < w || hnew < h)
    {
        typedef typename SrcAccessor::value_type SRCVT;
        typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
        typedef typename BasicImage<TMPTYPE>::Iterator TMPITER;
    
        BasicImage<TMPTYPE> t(w,h);
        TMPITER it = t.upperLeft();
    
        if(wnew < w)
        {
            recursiveSmoothX(is, iend, sa, 
                    it, t.accessor(), (double)w/wnew/scale);
            
            if(hnew < h)
            {
               recursiveSmoothY(it, t.lowerRight(), t.accessor(), 
                    it, t.accessor(), (double)h/hnew/scale);
            }
        }
        else
        {
           recursiveSmoothY(is, iend, sa, 
                    it, t.accessor(), (double)h/hnew/scale);
        }
        
        resizeImageInternalSplineInterpolation(it, t.lowerRight(), t.accessor(),
                                               id, idend, da);
    }
    else
    {
        resizeImageInternalSplineInterpolation(is, iend, sa, id, idend, da);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void 
resizeImageSplineInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageSplineInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

/********************************************************/
/*                                                      */
/*             resizeLineLinearInterpolation            */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void 
resizeLineLinearInterpolation(SrcIterator i1, SrcIterator iend, SrcAccessor as,
                           DestIterator id, DestIterator idend, DestAccessor ad)
{
    int wold = iend - i1;
    int wnew = idend - id;
    
    if((wold <= 1) || (wnew <= 1)) return; // oder error ? 
    
    typedef 
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    ad.set(DestTraits::fromRealPromote(as(i1)), id);
    ++id;
    
    --iend, --idend;
    ad.set(DestTraits::fromRealPromote(as(iend)), idend);
    
    double dx = (double)(wold - 1) / (wnew - 1);
    double x = dx;
    
    for(; id != idend; ++id, x += dx)
    {
        if(x >= 1.0) 
        {
            int xx = (int)x;
            i1 += xx;
            x -= (double)xx;
        }
        double x1 = 1.0 - x;

        ad.set(DestTraits::fromRealPromote(x1 * as(i1) + x * as(i1, 1)), id);
    }
}

/********************************************************/
/*                                                      */
/*           resizeImageLinearInterpolation             */
/*                                                      */
/********************************************************/

/** \brief Resize image using linear interpolation.

    The function uses the standard separable bilinear interpolation algorithm to 
    obtain a good compromize between quality and speed.
    
    The range must of both the input and output images (resp. regions)
    must be given. Both images must have a size of at
    least 2x2. The scaling factors are then calculated
    accordingly. If the source image is larger than the destination, it
    is smoothed (band limited) using a recursive
    exponential filter. The source value_type (SrcAccessor::value_type) must 
    be a linear space, i.e. it must support addition, multiplication 
    with a scalar real number and \ref NumericTraits "NumericTraits".
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void 
        resizeImageLinearInterpolation(
              SrcImageIterator is, SrcImageIterator iend, SrcAccessor sa,
          DestImageIterator id, DestImageIterator idend, DestAccessor da)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void 
        resizeImageLinearInterpolation(
              triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
          triple<DestImageIterator, DestImageIterator, DestAccessor> dest)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="resizeimage_8hxx-source.html">vigra/resizeimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::resizeImageLinearInterpolation(
               src.upperLeft(), src.lowerRight(), src.accessor(), 
               dest.upperLeft(), dest.lowerRight(), dest.accessor());
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft, src_lowerright;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote 
                             u = src_accessor(src_upperleft),
                 v = src_accessor(src_upperleft, 1);
    double d;
    
    u = d * v;
    u = u + v;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(u), 
    dest_upperleft);

    \endcode
    
    <b> Preconditions:</b>
    
    \code
    src_lowerright.x - src_upperleft.x > 1
    src_lowerright.y - src_upperleft.y > 1
    dest_lowerright.x - dest_upperleft.x > 1
    dest_lowerright.y - dest_upperleft.y > 1
    \endcode
    
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
resizeImageLinearInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    int w = iend.x - is.x;
    int h = iend.y - is.y;
    
    int wnew = idend.x - id.x;
    int hnew = idend.y - id.y;
    
    vigra_precondition((w > 1) && (h > 1),
                 "resizeImageLinearInterpolation(): "
         "Source image to small.\n");
    vigra_precondition((wnew > 1) && (hnew > 1),
                 "resizeImageLinearInterpolation(): "
         "Destination image to small.\n");
        
    double const scale = 2.0;
    
    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
    
    BasicImage<TMPTYPE> tmp(w, hnew);
    std::vector<TMPTYPE> coeffs((h > w) ? h : w);
    BasicImage<TMPTYPE> line((h > w) ? h : w, 1);
    
    int x,y;
    
    typename BasicImage<TMPTYPE>::Iterator yt = tmp.upperLeft();
    RowIterator<typename BasicImage<TMPTYPE>::Iterator> 
                                           lt(line.upperLeft());
    
    for(x=0; x<w; ++x, ++is.x, ++yt.x)
    {
        ColumnIterator<SrcIterator> c1(is);
        ColumnIterator<typename BasicImage<TMPTYPE>::Iterator> ct(yt);
    
        if(hnew < h)
        {
            resizeBandLimitLine(c1, c1 + h, accessorAdapter(c1, sa), coeffs, 
                 lt, accessorAdapter(lt, line.accessor()), (double)h/hnew/scale);

            resizeLineLinearInterpolation(lt, lt + h, accessorAdapter(lt, line.accessor()), 
                                          ct, ct + hnew, accessorAdapter(ct, tmp.accessor()));
        }
        else
        {
            resizeLineLinearInterpolation(c1, c1 + h, accessorAdapter(c1, sa), 
                                          ct, ct + hnew, accessorAdapter(ct, tmp.accessor()));
        }
    }
    
    yt = tmp.upperLeft();
    
    for(y=0; y < hnew; ++y, ++yt.y, ++id.y)
    {
        RowIterator<DestIterator> rd(id);
        RowIterator<typename BasicImage<TMPTYPE>::Iterator> rt(yt);
    
        if(wnew < w)
        {
            resizeBandLimitLine(rt, rt + w, accessorAdapter(rt, tmp.accessor()), coeffs,
                              lt, accessorAdapter(lt, line.accessor()), (double)w/wnew/scale);

            resizeLineLinearInterpolation(lt, lt + w, accessorAdapter(lt, line.accessor()), 
                                          rd, rd + wnew, accessorAdapter(rd, da));
        }
        else
        {
            resizeLineLinearInterpolation(rt, rt + w, accessorAdapter(rt, tmp.accessor()),  
                                          rd, rd + wnew, accessorAdapter(rd, da));
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void 
resizeImageLinearInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageLinearInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

/********************************************************/
/*                                                      */
/*               resizeLineNoInterpolation              */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void 
resizeLineNoInterpolation(SrcIterator i1, SrcIterator iend, SrcAccessor as,
                           DestIterator id, DestIterator idend, DestAccessor ad)
{
    int wold = iend - i1;
    int wnew = idend - id;
    
    if((wold <= 1) || (wnew <= 1)) return; // oder error ? 
    
    ad.set(as(i1), id);
    ++id;
    
    --iend, --idend;
    ad.set(as(iend), idend);
    
    double dx = (double)(wold - 1) / (wnew - 1);
    double x = dx;
    
    for(; id != idend; ++id, x += dx)
    {
    if(x >= 1.0) 
    {
        int xx = (int)x;
        i1 += xx;
        x -= (double)xx;
    }

    ad.set(as(i1), id);
    }
}

/********************************************************/
/*                                                      */
/*              resizeImageNoInterpolation              */
/*                                                      */
/********************************************************/

/** \brief Resize image by repeating the nearest pixel values.

    This algorithm is very fast and does not require any arithmetic on the pixel types.
    
    The range must of both the input and output images (resp. regions)
    must be given. Both images must have a size of at
    least 2x2. The scaling factors are then calculated
    accordingly. Destiniation pixels are directly copied from the appropriate
    source pixels. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void 
        resizeImageNoInterpolation(
              SrcImageIterator is, SrcImageIterator iend, SrcAccessor sa,
          DestImageIterator id, DestImageIterator idend, DestAccessor da)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void 
        resizeImageNoInterpolation(
              triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
          triple<DestImageIterator, DestImageIterator, DestAccessor> dest)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="resizeimage_8hxx-source.html">vigra/resizeimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::resizeImageNoInterpolation(
               src.upperLeft(), src.lowerRight(), src.accessor(), 
               dest.upperLeft(), dest.lowerRight(), dest.accessor());
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft, src_lowerright;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    dest_accessor.set(src_accessor(src_upperleft), dest_upperleft);

    \endcode
    
    <b> Preconditions:</b>
    
    \code
    src_lowerright.x - src_upperleft.x > 1
    src_lowerright.y - src_upperleft.y > 1
    dest_lowerright.x - dest_upperleft.x > 1
    dest_lowerright.y - dest_upperleft.y > 1
    \endcode
    
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
resizeImageNoInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    int w = iend.x - is.x;
    int h = iend.y - is.y;
    
    int wnew = idend.x - id.x;
    int hnew = idend.y - id.y;
    
    vigra_precondition((w > 1) && (h > 1),
                 "resizeImageNoInterpolation(): "
         "Source image to small.\n");
    vigra_precondition((wnew > 1) && (hnew > 1),
                 "resizeImageNoInterpolation(): "
         "Destination image to small.\n");
        
    typedef typename SrcAccessor::value_type SRCVT;

    BasicImage<SRCVT> tmp(w, hnew);
    
    int x,y;
    
    typename BasicImage<SRCVT>::Iterator yt = tmp.upperLeft();

    for(x=0; x<w; ++x, ++is.x, ++yt.x)
    {
        ColumnIterator<SrcIterator> c1(is);
        ColumnIterator<typename BasicImage<SRCVT>::Iterator> ct(yt);
    
    resizeLineNoInterpolation(c1, c1 + h, accessorAdapter(c1, sa), 
                              ct, ct + hnew, accessorAdapter(ct, tmp.accessor()));
    }
    
    yt = tmp.upperLeft();
    
    for(y=0; y < hnew; ++y, ++yt.y, ++id.y)
    {
        RowIterator<DestIterator> rd(id);
        RowIterator<typename BasicImage<SRCVT>::Iterator> rt(yt);
    
    resizeLineNoInterpolation(rt, rt + w, accessorAdapter(rt, tmp.accessor()),  
                  rd, rd + wnew, accessorAdapter(rd, da));
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void 
resizeImageNoInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageNoInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

//@}

} // namespace vigra

#endif // VIGRA_RESIZEIMAGE_HXX
