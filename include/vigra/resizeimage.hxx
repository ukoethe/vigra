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


#ifndef VIGRA_RESIZEIMAGE_HXX
#define VIGRA_RESIZEIMAGE_HXX

#include <vector>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdimage.hxx"
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
    typedef BasicImage<SRCVT> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<SRCVT> tmp(w, hnew);

    int x,y;

    typename BasicImage<SRCVT>::Iterator yt = tmp.upperLeft();

    for(x=0; x<w; ++x, ++is.x, ++yt.x)
    {
        typename SrcIterator::column_iterator c1 = is.columnIterator();
        typename TmpImageIterator::column_iterator ct = yt.columnIterator();

        resizeLineNoInterpolation(c1, c1 + h, sa, ct, ct + hnew, tmp.accessor());
    }

    yt = tmp.upperLeft();

    for(y=0; y < hnew; ++y, ++yt.y, ++id.y)
    {
        typename DestIterator::row_iterator rd = id.rowIterator();
        typename TmpImageIterator::row_iterator rt = yt.rowIterator();

        resizeLineNoInterpolation(rt, rt + w, tmp.accessor(), rd, rd + wnew, da);
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
    typedef BasicImage<TMPTYPE> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<TMPTYPE> tmp(w, hnew);
    BasicImage<TMPTYPE> line((h > w) ? h : w, 1);

    int x,y;

    typename BasicImage<TMPTYPE>::Iterator yt = tmp.upperLeft();
    typename TmpImageIterator::row_iterator lt = line.upperLeft().rowIterator();

    for(x=0; x<w; ++x, ++is.x, ++yt.x)
    {
        typename SrcIterator::column_iterator c1 = is.columnIterator();
        typename TmpImageIterator::column_iterator ct = yt.columnIterator();

        if(hnew < h)
        {
            recursiveSmoothLine(c1, c1 + h, sa,
                 lt, line.accessor(), (double)h/hnew/scale);

            resizeLineLinearInterpolation(lt, lt + h, line.accessor(),
                                          ct, ct + hnew, tmp.accessor());
        }
        else
        {
            resizeLineLinearInterpolation(c1, c1 + h, sa,
                                          ct, ct + hnew, tmp.accessor());
        }
    }

    yt = tmp.upperLeft();

    for(y=0; y < hnew; ++y, ++yt.y, ++id.y)
    {
        typename DestIterator::row_iterator rd = id.rowIterator();
        typename TmpImageIterator::row_iterator rt = yt.rowIterator();

        if(wnew < w)
        {
            recursiveSmoothLine(rt, rt + w, tmp.accessor(),
                              lt, line.accessor(), (double)w/wnew/scale);

            resizeLineLinearInterpolation(lt, lt + w, line.accessor(),
                                          rd, rd + wnew, da);
        }
        else
        {
            resizeLineLinearInterpolation(rt, rt + w, tmp.accessor(),
                                          rd, rd + wnew, da);
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
/*              CubicFIRInterpolationKernel             */
/*                                                      */
/********************************************************/

class CubicFIRInterpolationKernel
{
public:
    double operator[] (double x) const
    {
        x = fabs(x);
        if (x <= 1.0)
        {
            return 1.0 + x * x * (-2.5 + 1.5 * x);
        }
        else if (x >= 2.0)
        {
            return 0.0;
        }
        else
        {
            return 2.0 + x * (-4.0 + x * (2.5 -0.5 * x));
        }
    }

    int radius() const
        {return 2;}
};

/***************************************************************/
/*                                                             */
/*               resizeLineCubicFIRInterpolation               */
/*                                                             */
/***************************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void resizeLineCubicFIRInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                           DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TMPTYPE;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    int src_width = src_iter_end - src_iter;
    int dest_width = dest_iter_end - dest_iter;
    double dx =  (double)(src_width - 1) / (dest_width - 1);

    CubicFIRInterpolationKernel kernel;

    dest_acc.set(src_acc(src_iter), dest_iter);
    dest_iter++;
    for (int i = 1; i < dest_width-1; i++, dest_iter++)
    {
        double x = dx * i;
        int i_old = (int)x;
        double t = x - i_old;
        TMPTYPE value;

        if (i_old == 0)
        {
            value = kernel[t] * src_acc(src_iter, i_old) + kernel[1.0-t] * src_acc(src_iter, i_old + 1)
                    + kernel[2.0-t] * src_acc(src_iter, i_old + 2) + kernel[1.0 + t] * src_acc(src_iter, i_old + 1);
        }
        else if (i_old == src_width-2)
        {
            value = kernel[t] * src_acc(src_iter, i_old) + kernel[1.0-t] * src_acc(src_iter, i_old + 1)
                    + kernel[2.0-t] * src_acc(src_iter, i_old) + kernel[1.0 + t] * src_acc(src_iter, i_old - 1);
        }
        else
        {
            value = kernel[t] * src_acc(src_iter, i_old) + kernel[1.0-t] * src_acc(src_iter, i_old + 1)
                    + kernel[2.0-t] * src_acc(src_iter, i_old + 2) + kernel[1.0 + t] * src_acc(src_iter, i_old - 1);
        }

        dest_acc.set(DestTraits::fromRealPromote(value), dest_iter);
    }
    dest_acc.set(src_acc(--src_iter_end), --dest_iter_end);
}



/*****************************************************************/
/*                                                               */
/*              resizeImageCubicFIRInterpolation                 */
/*                                                               */
/*****************************************************************/



template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeImageCubicFIRInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                      DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc)
{
    int width_old = src_iter_end.x - src_iter.x;
    int height_old = src_iter_end.y - src_iter.y;

    int width_new = dest_iter_end.x - dest_iter.x;
    int height_new = dest_iter_end.y - dest_iter.y;
    //    double dx =  (double)(width_old - 1) / (width_new - 1);
    //    double dy =  (double)(height_old - 1) / (height_new - 1);
    double const scale = 2.0;

    vigra_precondition((width_old > 1) && (height_old > 1),
                 "resizeImageCubicFIRInterpolation(): "
                 "Source image to small.\n");

    vigra_precondition((width_new > 1) && (height_new > 1),
                 "resizeImageCubicFIRInterpolation(): "
                  "Destination image to small.\n");

    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
    typedef BasicImage<TMPTYPE> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<TMPTYPE> tmp(width_old, height_new);
    BasicImage<TMPTYPE> line((height_old > width_old) ? height_old : width_old, 1);
    typename BasicImage<TMPTYPE>::Accessor tmp_acc = tmp.accessor();

    int x,y;

    typename BasicImage<TMPTYPE>::Iterator y_tmp = tmp.upperLeft();
    typename TmpImageIterator::row_iterator line_tmp = line.upperLeft().rowIterator();

    for(x=0; x<width_old; ++x, ++src_iter.x, ++y_tmp.x)
    {
        typename SrcIterator::column_iterator c_src = src_iter.columnIterator();
        typename TmpImageIterator::column_iterator c_tmp = y_tmp.columnIterator();

        if(height_new < height_old)
        {
            recursiveSmoothLine(c_src, c_src + height_old, src_acc,
                 line_tmp, line.accessor(), (double)height_old/height_new/scale);

            resizeLineCubicFIRInterpolation(line_tmp, line_tmp + height_old, line.accessor(),
                                            c_tmp, c_tmp + height_new, tmp_acc);
        }
        else
        {

            resizeLineCubicFIRInterpolation(c_src, c_src + height_old, src_acc,
                                            c_tmp, c_tmp + height_new, tmp_acc);
        }
    }

    y_tmp = tmp.upperLeft();

    typename BasicImage<SRCVT>::Iterator dest = dest_iter ;

    for(y=0; y < height_new; ++y, ++y_tmp.y, ++dest_iter.y)
    {
        typename DestIterator::row_iterator r_dest = dest_iter.rowIterator();
        typename TmpImageIterator::row_iterator r_tmp = y_tmp.rowIterator();

        if(width_new < width_old)
        {
            recursiveSmoothLine(r_tmp, r_tmp + width_old, tmp.accessor(),
                              line_tmp, line.accessor(), (double)width_old/width_new/scale);

            resizeLineCubicFIRInterpolation(line_tmp, line_tmp + width_old, line.accessor(),
                                            r_dest, r_dest + width_new, dest_acc);
        }
        else
        {
            resizeLineCubicFIRInterpolation(r_tmp, r_tmp + width_old, tmp_acc,
                                            r_dest, r_dest + width_new, dest_acc);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void
resizeImageCubicFIRInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageCubicFIRInterpolation(src.first, src.second, src.third,
                                     dest.first, dest.second, dest.third);
}

/********************************************************/
/*                                                      */
/*                  CubicBSplineKernel                  */
/*                                                      */
/********************************************************/
class CubicBSplineKernel
{
  public:
    double operator[] (double x) const
    {
        x = fabs(x);
        if (x < 1.0)
        {
            return 2.0/3.0 - x*x*(1.0 - x/2.0);
        }
        else if (x >= 2.0)
        {
            return 0.0;
        }
        else
        {
            double t = 2.0 - x;
            return t*t*t/6.0;
        }
    }

    double dx(double x) const
    {
        double ax = fabs(x);
        if (ax < 1.0)
        {
            return  x*(1.5*ax - 2.0);
        }
        else if (ax >= 2.0)
        {
            return 0.0;
        }
        else
        {
            double t = 2.0 - ax;
            return x < 0.0
                ?  t*t/2.0
                : -t*t/2.0;
        }
    }

    double dxx(double x) const
    {
        x = fabs(x);
        if (x < 1.0)
        {
            return  3.0*x - 2.0;
        }
        else if (x >= 2.0)
        {
            return 0.0;
        }
        else
        {
            return 2.0 - x;
        }
    }

    int radius() const
        {return 2;}
};

/******************************************************************/
/*                                                                */
/*             resizeLineCubicIIRInterpolation                    */
/*                                                                */
/******************************************************************/


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void resizeLineCubicIIRInterpolation(
    SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
    DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TMPTYPE;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    int src_width = src_iter_end - src_iter;
    int dest_width = dest_iter_end - dest_iter;

    double dx = (double)(src_width-1)/(dest_width-1);

    std::vector<TMPTYPE> tmp(src_width);
    typename std::vector<TMPTYPE>::iterator tmp_iter = tmp.begin();
    typename std::vector<TMPTYPE>::iterator tmp_iter_end = tmp.end();
    StandardAccessor<TMPTYPE> tmp_acc;

    recursiveFilterLine(src_iter, src_iter_end, src_acc, tmp_iter, tmp_acc,
                        sqrt(3.0) - 2.0, BORDER_TREATMENT_REFLECT);

    CubicBSplineKernel kernel;
    dest_acc.set(DestTraits::fromRealPromote(kernel[0.0] * tmp[0] + 2.0 * kernel[1.0] * tmp[1]),
                 dest_iter);
    dest_iter++;
    for (int i = 1; i < dest_width-1; i++, dest_iter++ )
    {
        double x = dx * i;
        int i_old = (int)x;
        double t = x - i_old;
        TMPTYPE value;

        if (i_old == 0)
        {
            value = kernel[t] * tmp[i_old] + kernel[1.0-t] * tmp[i_old + 1]
                    + kernel[2.0-t] * tmp[i_old + 2] + kernel[1.0 + t] * tmp[i_old + 1];

        }
        else if (i_old == tmp.size()-2)
        {
            value = kernel[t] * tmp[i_old] + kernel[1.0-t] * tmp[i_old + 1]
                    + kernel[2.0-t] * tmp[i_old] + kernel[1.0 + t] * tmp[i_old - 1];

        }
        else
        {
            value = kernel[t] * tmp[i_old] + kernel[1.0-t] * tmp[i_old + 1]
                    + kernel[2.0-t] * tmp[i_old + 2] + kernel[1.0 + t] * tmp[i_old - 1];

        }
        dest_acc.set(DestTraits::fromRealPromote(value), dest_iter);
    }
    dest_acc.set(DestTraits::fromRealPromote(kernel[0.0] * tmp[src_width-1]
                 + 2.0 * kernel[1.0] * tmp[src_width-2]), (--dest_iter_end));
}


/*****************************************************************/
/*                                                               */
/*            resizeImageCubicIIRInterpolation                   */
/*                                                               */
/*****************************************************************/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeImageCubicIIRInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                      DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc)
{

    int width_old = src_iter_end.x - src_iter.x;
    int height_old = src_iter_end.y - src_iter.y;

    int width_new = dest_iter_end.x - dest_iter.x;
    int height_new = dest_iter_end.y - dest_iter.y;
    //    double dx =  (double)(width_old - 1) / (width_new - 1);
    //    double dy =  (double)(height_old - 1) / (height_new - 1);
    double const scale = 2.0;

    vigra_precondition((width_old > 1) && (height_old > 1),
                 "resizeImageCubicFIRInterpolation(): "
                 "Source image to small.\n");

    vigra_precondition((width_new > 1) && (height_new > 1),
                 "resizeImageCubicFIRInterpolation(): "
                 "Destination image to small.\n");

    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
    typedef BasicImage<TMPTYPE> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<TMPTYPE> tmp(width_old, height_new);

    BasicImage<TMPTYPE> line((height_old > width_old) ? height_old : width_old, 1);
    typename BasicImage<TMPTYPE>::Accessor tmp_acc = tmp.accessor();

    int x,y;

    typename BasicImage<TMPTYPE>::Iterator y_tmp = tmp.upperLeft();
    typename TmpImageIterator::row_iterator line_tmp = line.upperLeft().rowIterator();

    for(x=0; x<width_old; ++x, ++src_iter.x, ++y_tmp.x)
    {

        typename SrcIterator::column_iterator c_src = src_iter.columnIterator();
        typename TmpImageIterator::column_iterator c_tmp = y_tmp.columnIterator();
        if(height_new < height_old)
        {
            recursiveSmoothLine(c_src, c_src + height_old, src_acc,
                 line_tmp, line.accessor(), (double)height_old/height_new/scale);

            resizeLineCubicIIRInterpolation(line_tmp, line_tmp + height_old, line.accessor(),
                                            c_tmp, c_tmp + height_new, tmp_acc);
        }
        else
        {

            resizeLineCubicIIRInterpolation(c_src, c_src + height_old, src_acc,
                                            c_tmp, c_tmp + height_new, tmp_acc);
        }
    }

    y_tmp = tmp.upperLeft();

    typename BasicImage<SRCVT>::Iterator dest = dest_iter ;

    for(y=0; y < height_new; ++y, ++y_tmp.y, ++dest_iter.y)
    {
        typename DestIterator::row_iterator r_dest = dest_iter.rowIterator();
        typename TmpImageIterator::row_iterator r_tmp = y_tmp.rowIterator();
        if(width_new < width_old)
        {
            recursiveSmoothLine(r_tmp, r_tmp + width_old, tmp.accessor(),
                              line_tmp, line.accessor(), (double)width_old/width_new/scale);

            resizeLineCubicIIRInterpolation(line_tmp, line_tmp + width_old, line.accessor(),
                                            r_dest, r_dest + width_new, dest_acc);
        }
        else
        {

            resizeLineCubicIIRInterpolation(r_tmp, r_tmp + width_old, tmp_acc,
                                            r_dest, r_dest + width_new, dest_acc);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void
resizeImageCubicIIRInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageCubicIIRInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
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
        typename SrcIterator::row_iterator sr = in.rowIterator();
        typename TMPITER::row_iterator dr = idx.rowIterator();
        resizeImageInternalSplineGradient(sr, sr+w, sa,
                                          it, r, dr);
    }

    in = is;

    // calculate y derivatives
    for(x=0; x<w; ++x, ++in.x, ++idy.x)
    {
        typename SrcIterator::column_iterator sc = in.columnIterator();
        typename TMPITER::column_iterator dc = idy.columnIterator();
        resizeImageInternalSplineGradient(sc, sc+h, sa,
                                          it, r, dc);
    }

    in = is;
    idy = dy.upperLeft();

    // calculate mixed derivatives
    for(y=0; y<h; ++y, ++idy.y, ++idxy.y)
    {
        typename TMPITER::row_iterator sr = idy.rowIterator();
        typename TMPITER::row_iterator dr = idxy.rowIterator();
        resizeImageInternalSplineGradient(sr, sr+w, ta,
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
    <i>"Grundlagen der geometrischen Datenverarbeitung"</i>, Teubner, 1992] to obtain
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


//@}

} // namespace vigra

#endif // VIGRA_RESIZEIMAGE_HXX
