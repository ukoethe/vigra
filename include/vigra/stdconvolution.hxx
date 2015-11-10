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


#ifndef VIGRA_STDCONVOLUTION_HXX
#define VIGRA_STDCONVOLUTION_HXX

#include <cmath>
#include "stdimage.hxx"
#include "bordertreatment.hxx"
#include "separableconvolution.hxx"
#include "utilities.hxx"
#include "sized_int.hxx"
#include "multi_iterator.hxx"
#include "multi_shape.hxx"

namespace vigra {

template <class ARITHTYPE>
class Kernel2D;

/** \addtogroup CommonConvolutionFilters
*/
//@{

    // documentation is in convolution.hxx
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void convolveImage(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                   DestIterator dest_ul, DestAccessor dest_acc,
                   KernelIterator ki, KernelAccessor ak,
                   Diff2D kul, Diff2D klr, BorderTreatmentMode border)
{
    vigra_precondition((border == BORDER_TREATMENT_CLIP    ||
                        border == BORDER_TREATMENT_AVOID   ||
                        border == BORDER_TREATMENT_REFLECT ||
                        border == BORDER_TREATMENT_REPEAT  ||
                        border == BORDER_TREATMENT_WRAP    ||
                        border == BORDER_TREATMENT_ZEROPAD),
                       "convolveImage():\n"
                       "  Border treatment must be one of follow treatments:\n"
                       "  - BORDER_TREATMENT_CLIP\n"
                       "  - BORDER_TREATMENT_AVOID\n"
                       "  - BORDER_TREATMENT_REFLECT\n"
                       "  - BORDER_TREATMENT_REPEAT\n"
                       "  - BORDER_TREATMENT_WRAP\n"
                       "  - BORDER_TREATMENT_ZEROPAD\n");

    vigra_precondition(kul.x <= 0 && kul.y <= 0,
                       "convolveImage(): coordinates of "
                       "kernel's upper left must be <= 0.");
    vigra_precondition(klr.x >= 0 && klr.y >= 0,
                       "convolveImage(): coordinates of "
                       "kernel's lower right must be >= 0.");

    // use traits to determine SumType as to prevent possible overflow
    typedef typename
        PromoteTraits<typename SrcAccessor::value_type,
                      typename KernelAccessor::value_type>::Promote SumType;
    typedef typename
        NumericTraits<typename KernelAccessor::value_type>::RealPromote KernelSumType;
    typedef typename DestAccessor::value_type DestType;

    // calculate width and height of the image
    int w = src_lr.x - src_ul.x;
    int h = src_lr.y - src_ul.y;

    // calculate width and height of the kernel
    int kernel_width  = klr.x - kul.x + 1;
    int kernel_height = klr.y - kul.y + 1;

    vigra_precondition(w >= std::max(klr.x, -kul.x) + 1 && h >= std::max(klr.y, -kul.y) + 1,
                       "convolveImage(): kernel larger than image.");

    KernelSumType norm = KernelSumType();
    if(border == BORDER_TREATMENT_CLIP)
    {
        // calculate the sum of the kernel elements for renormalization
        KernelIterator yk  = ki + klr;

        // determine sum within kernel (= norm)
        for(int y = 0; y < kernel_height; ++y, --yk.y)
        {
            KernelIterator xk  = yk;
            for(int x = 0; x < kernel_width; ++x, --xk.x)
            {
                norm += ak(xk);
            }
        }
        vigra_precondition(norm != NumericTraits<KernelSumType>::zero(),
            "convolveImage(): Cannot use BORDER_TREATMENT_CLIP with a DC-free kernel");
    }

    DestIterator yd = dest_ul;
    SrcIterator ys = src_ul;

    // iterate over the interior part
    for(int y=0; y<h; ++y, ++ys.y, ++yd.y)
    {
        // create x iterators
        DestIterator xd(yd);
        SrcIterator xs(ys);

        for(int x=0; x < w; ++x, ++xs.x, ++xd.x)
        {
            // init the sum
            SumType sum = NumericTraits<SumType>::zero();
            KernelIterator ykernel  = ki + klr;

            if(x >= klr.x && y >= klr.y && x < w + kul.x && y < h + kul.y)
            {
                // kernel is entirely inside the image
                SrcIterator yys = xs - klr;
                SrcIterator yyend = xs - kul;

                for(; yys.y <= yyend.y; ++yys.y, --ykernel.y)
                {
                    typename SrcIterator::row_iterator xxs = yys.rowIterator();
                    typename SrcIterator::row_iterator xxe = xxs + kernel_width;
                    typename KernelIterator::row_iterator xkernel= ykernel.rowIterator();

                    for(; xxs < xxe; ++xxs, --xkernel)
                    {
                        sum += ak(xkernel) * src_acc(xxs);
                    }
                }
            }
            else if(border == BORDER_TREATMENT_REPEAT)
            {
                Diff2D diff;
                for(int yk = klr.y; yk >= kul.y; --yk, --ykernel.y)
                {
                    diff.y = std::min(std::max(y - yk, 0), h-1);
                    typename KernelIterator::row_iterator xkernel  = ykernel.rowIterator();

                    for(int xk = klr.x; xk >= kul.x; --xk, --xkernel)
                    {
                        diff.x = std::min(std::max(x - xk, 0), w-1);
                        sum += ak(xkernel) * src_acc(src_ul, diff);
                    }
                }
            }
            else if(border == BORDER_TREATMENT_REFLECT)
            {
                Diff2D diff;
                for(int yk = klr.y; yk >= kul.y; --yk , --ykernel.y)
                {
                    diff.y = abs(y - yk);
                    if(diff.y >= h)
                        diff.y = 2*h - 2 - diff.y;
                    typename KernelIterator::row_iterator xkernel  = ykernel.rowIterator();

                    for(int xk = klr.x; xk >= kul.x; --xk, --xkernel)
                    {
                        diff.x = abs(x - xk);
                        if(diff.x >= w)
                            diff.x = 2*w - 2 - diff.x;
                        sum += ak(xkernel) * src_acc(src_ul, diff);
                    }
                }
            }
            else if(border == BORDER_TREATMENT_WRAP)
            {
                Diff2D diff;
                for(int yk = klr.y; yk >= kul.y; --yk, --ykernel.y)
                {
                    diff.y = (y - yk + h) % h;
                    typename KernelIterator::row_iterator xkernel  = ykernel.rowIterator();

                    for(int xk = klr.x; xk >= kul.x; --xk, --xkernel)
                    {
                        diff.x = (x - xk + w) % w;
                        sum += ak(xkernel) * src_acc(src_ul, diff);
                    }
                }
            }
            else if(border == BORDER_TREATMENT_CLIP)
            {
                KernelSumType ksum = NumericTraits<KernelSumType>::zero();
                Diff2D diff;
                for(int yk = klr.y; yk >= kul.y; --yk, --ykernel.y)
                {
                    diff.y = y - yk;
                    if(diff.y < 0 || diff.y >= h)
                        continue;
                    typename KernelIterator::row_iterator xkernel  = ykernel.rowIterator();

                    for(int xk = klr.x; xk >= kul.x; --xk, --xkernel)
                    {
                        diff.x = x - xk;
                        if(diff.x < 0 || diff.x >= w)
                            continue;
                        ksum += ak(xkernel);
                        sum += ak(xkernel) * src_acc(src_ul, diff);
                    }
                }

                sum *= norm / ksum;
            }
            else if(border == BORDER_TREATMENT_ZEROPAD)
            {
                Diff2D diff;
                for(int yk = klr.y; yk >= kul.y; --yk, --ykernel.y)
                {
                    diff.y = y - yk;
                    if(diff.y < 0 || diff.y >= h)
                        continue;
                    typename KernelIterator::row_iterator xkernel  = ykernel.rowIterator();

                    for(int xk = klr.x; xk >= kul.x; --xk, --xkernel)
                    {
                        diff.x = x - xk;
                        if(diff.x < 0 || diff.x >= w)
                            continue;
                        sum += ak(xkernel) * src_acc(src_ul, diff);
                    }
                }
            }
            else if(border == BORDER_TREATMENT_AVOID)
            {
                continue;
            }

            // store convolution result in destination pixel
            dest_acc.set(detail::RequiresExplicitCast<DestType>::cast(sum), xd);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
inline void
convolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
              BorderTreatmentMode> kernel)
{
    convolveImage(src.first, src.second, src.third,
                  dest.first, dest.second,
                  kernel.first, kernel.second, kernel.third,
                  kernel.fourth, kernel.fifth);
}

template <class T1, class S1,
          class T2, class S2,
          class T3>
inline void
convolveImage(MultiArrayView<2, T1, S1> const & src,
              MultiArrayView<2, T2, S2> dest,
              Kernel2D<T3> const & kernel)
{
    vigra_precondition(src.shape() == dest.shape(),
        "convolveImage(): shape mismatch between input and output.");
    convolveImage(srcImageRange(src),
                  destImage(dest),
                  kernel2d(kernel));
}

/** \brief Performs a 2-dimensional normalized convolution, i.e. convolution with a mask image.

    This functions computes
    <a href ="http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/PIRODDI1/NormConv/NormConv.html">normalized
    convolution</a> as defined in
    Knutsson, H. and Westin, C-F.: <i>Normalized and differential convolution:
    Methods for Interpolation and Filtering of incomplete and uncertain data</i>.
    Proc. of the IEEE Conf. on Computer Vision and Pattern Recognition, 1993, 515-523.

    The mask image must be binary and encodes which pixels of the original image
    are valid. It is used as follows:
    Only pixel under the mask are used in the calculations. Whenever a part of the
    kernel lies outside the mask, it is ignored, and the kernel is renormalized to its
    original norm (analogous to the CLIP \ref BorderTreatmentMode). Thus, a useful convolution
    result is computed whenever <i>at least one valid pixel is within the current window</i>
    Thus, destination pixels not under the mask still receive a value if they are <i>near</i>
    the mask. Therefore, this algorithm is useful as an interpolator of sparse input data.
    If you are only interested in the destination values under the mask, you can perform
    a subsequent \ref copyImageIf().

    The KernelIterator must point to the center of the kernel, and
    the kernel's size is given by its upper left (x and y of distance <= 0) and
    lower right (distance >= 0) corners. The image must always be larger than the
    kernel. At those positions where the kernel does not completely fit
    into the image, the specified \ref BorderTreatmentMode is
    applied. Only BORDER_TREATMENT_CLIP and BORDER_TREATMENT_AVOID are currently
    supported.

    The images's pixel type (SrcAccessor::value_type) must be a
    linear space over the kernel's value_type (KernelAccessor::value_type),
    i.e. addition of source values, multiplication with kernel values,
    and NumericTraits must be defined.
    The kernel's value_type must be an algebraic field,
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must
    be defined.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class TM, class SM,
                  class T3>
        void
        normalizedConvolveImage(MultiArrayView<2, T1, S1> const & src,
                                MultiArrayView<2, TM, SM> const & mask,
                                MultiArrayView<2, T2, S2> dest,
                                Kernel2D<T3> const & kernel);
    }
    \endcode

    \deprecatedAPI{normalizedConvolveImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void
        normalizedConvolveImage(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                                MaskIterator mul, MaskAccessor am,
                                DestIterator dest_ul, DestAccessor dest_acc,
                                KernelIterator ki, KernelAccessor ak,
                                Diff2D kul, Diff2D klr, BorderTreatmentMode border);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void normalizedConvolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     pair<MaskIterator, MaskAccessor> mask,
                                     pair<DestIterator, DestAccessor> dest,
                                     tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
                                     BorderTreatmentMode> kernel);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/stdconvolution.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, float>          src(w,h), dest(w,h);
    MultiArray<2, unsigned char>  mask(w,h);
    ...
    // define 3x3 binomial filter
    vigra::Kernel2D<float> binom;
    binom.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =   // upper left and lower right
                         0.0625, 0.125, 0.0625,
                         0.125,  0.25,  0.125,
                         0.0625, 0.125, 0.0625;

    normalizedConvolveImage(src, mask, dest, binom);
    \endcode

    \deprecatedUsage{normalizedConvolveImage}
    \code
    vigra::FImage src(w,h), dest(w,h);
    vigra::CImage mask(w,h);
    ...

    // define 3x3 binomial filter
    vigra::Kernel2D<float> binom;

    binom.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =   // upper left and lower right
                         0.0625, 0.125, 0.0625,
                         0.125,  0.25,  0.125,
                         0.0625, 0.125, 0.0625;

    vigra::normalizedConvolveImage(srcImageRange(src), maskImage(mask), destImage(dest), kernel2d(binom));
    \endcode
    <b> Required Interface:</b>
    \code
    ImageIterator src_ul, src_lr;
    ImageIterator mul;
    ImageIterator dest_ul;
    ImageIterator ik;

    SrcAccessor src_accessor;
    MaskAccessor mask_accessor;
    DestAccessor dest_accessor;
    KernelAccessor kernel_accessor;

    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(src_ul);

    s = s + s;
    s = kernel_accessor(ik) * s;
    s -= s;

    if(mask_accessor(mul)) ...;

    dest_accessor.set(
    NumericTraits<DestAccessor::value_type>::fromRealPromote(s), dest_ul);

    NumericTraits<KernelAccessor::value_type>::RealPromote k = kernel_accessor(ik);

    k += k;
    k -= k;
    k = k / k;

    \endcode
    \deprecatedEnd

    <b> Preconditions:</b>

    <ul>
    <li> The image must be longer than the kernel radius: <tt>w > std::max(kernel.lowerRight().x, -kernel.upperLeft().x)</tt> and
         <tt>h > std::max(kernel.lowerRight().y, -kernel.upperLeft().y)</tt>.
    <li> The sum of kernel elements must be != 0.
    <li> <tt>border == BORDER_TREATMENT_CLIP || border == BORDER_TREATMENT_AVOID</tt>
    </ul>
*/
doxygen_overloaded_function(template <...> void normalizedConvolveImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
void
normalizedConvolveImage(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                        MaskIterator mul, MaskAccessor am,
                        DestIterator dest_ul, DestAccessor dest_acc,
                        KernelIterator ki, KernelAccessor ak,
                        Diff2D kul, Diff2D klr, BorderTreatmentMode border)
{
    vigra_precondition((border == BORDER_TREATMENT_CLIP  ||
                        border == BORDER_TREATMENT_AVOID),
                       "normalizedConvolveImage(): "
                       "Border treatment must be BORDER_TREATMENT_CLIP or BORDER_TREATMENT_AVOID.");

    vigra_precondition(kul.x <= 0 && kul.y <= 0,
                       "normalizedConvolveImage(): left borders must be <= 0.");
    vigra_precondition(klr.x >= 0 && klr.y >= 0,
                       "normalizedConvolveImage(): right borders must be >= 0.");

    // use traits to determine SumType as to prevent possible overflow
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SumType;
    typedef typename
        NumericTraits<typename KernelAccessor::value_type>::RealPromote KSumType;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    // calculate width and height of the image
    int w = src_lr.x - src_ul.x;
    int h = src_lr.y - src_ul.y;
    int kernel_width = klr.x - kul.x + 1;
    int kernel_height = klr.y - kul.y + 1;

    int x,y;
    int ystart = (border == BORDER_TREATMENT_AVOID) ?  klr.y : 0;
    int yend   = (border == BORDER_TREATMENT_AVOID) ? h+kul.y : h;
    int xstart = (border == BORDER_TREATMENT_AVOID) ?  klr.x : 0;
    int xend   = (border == BORDER_TREATMENT_AVOID) ? w+kul.x : w;

    // create y iterators
    DestIterator yd = dest_ul + Diff2D(xstart, ystart);
    SrcIterator ys = src_ul + Diff2D(xstart, ystart);
    MaskIterator ym = mul + Diff2D(xstart, ystart);

    KSumType norm = ak(ki);
    int xx, yy;
    KernelIterator yk  = ki + klr;
    for(yy=0; yy<kernel_height; ++yy, --yk.y)
    {
        KernelIterator xk  = yk;

        for(xx=0; xx<kernel_width; ++xx, --xk.x)
        {
            norm += ak(xk);
        }
    }
    norm -= ak(ki);


    for(y=ystart; y < yend; ++y, ++ys.y, ++yd.y, ++ym.y)
    {
        // create x iterators
        DestIterator xd(yd);
        SrcIterator xs(ys);
        MaskIterator xm(ym);

        for(x=xstart; x < xend; ++x, ++xs.x, ++xd.x, ++xm.x)
        {
            // how much of the kernel fits into the image ?
            int x0, y0, x1, y1;

            y0 = (y<klr.y) ? -y : -klr.y;
            y1 = (h-y-1<-kul.y) ? h-y-1 : -kul.y;
            x0 = (x<klr.x) ? -x : -klr.x;
            x1 = (w-x-1<-kul.x) ? w-x-1 : -kul.x;

            bool first = true;
            // init the sum
            SumType sum = NumericTraits<SumType>::zero();
            KSumType ksum = NumericTraits<KSumType>::zero();

            SrcIterator yys = xs + Diff2D(x0, y0);
            MaskIterator yym = xm + Diff2D(x0, y0);
            KernelIterator yk  = ki - Diff2D(x0, y0);

            int kernel_width, kernel_height;
            kernel_width = x1 - x0 + 1;
            kernel_height = y1 - y0 + 1;
            for(yy=0; yy<kernel_height; ++yy, ++yys.y, --yk.y, ++yym.y)
            {
                typename SrcIterator::row_iterator xxs = yys.rowIterator();
                typename SrcIterator::row_iterator xxend = xxs + kernel_width;
                typename MaskIterator::row_iterator xxm = yym.rowIterator();
                typename KernelIterator::row_iterator xk  = yk.rowIterator();

                for(xx=0; xxs < xxend; ++xxs, --xk, ++xxm)
                {
                    if(!am(xxm)) continue;

                    if(first)
                    {
                        sum = detail::RequiresExplicitCast<SumType>::cast(ak(xk) * src_acc(xxs));
                        ksum = ak(xk);
                        first = false;
                    }
                    else
                    {
                        sum = detail::RequiresExplicitCast<SumType>::cast(sum + ak(xk) * src_acc(xxs));
                        ksum += ak(xk);
                    }
                }
            }
            // store average in destination pixel
            if(ksum != NumericTraits<KSumType>::zero())
            {
                dest_acc.set(DestTraits::fromRealPromote(
                             detail::RequiresExplicitCast<SumType>::cast((norm / ksum) * sum)), xd);
            }
        }
    }
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
inline void
normalizedConvolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        pair<MaskIterator, MaskAccessor> mask,
                        pair<DestIterator, DestAccessor> dest,
                        tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
                        BorderTreatmentMode> kernel)
{
    normalizedConvolveImage(src.first, src.second, src.third,
                            mask.first, mask.second,
                            dest.first, dest.second,
                            kernel.first, kernel.second, kernel.third,
                            kernel.fourth, kernel.fifth);
}

template <class T1, class S1,
          class T2, class S2,
          class TM, class SM,
          class T3>
inline void
normalizedConvolveImage(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, TM, SM> const & mask,
                        MultiArrayView<2, T2, S2> dest,
                        Kernel2D<T3> const & kernel)
{
    vigra_precondition(src.shape() == mask.shape() && src.shape() == dest.shape(),
        "normalizedConvolveImage(): shape mismatch between input and output.");
    normalizedConvolveImage(srcImageRange(src),
                            maskImage(mask),
                            destImage(dest),
                            kernel2d(kernel));
}

/** \brief Deprecated name of 2-dimensional normalized convolution, i.e. convolution with a mask image.

    See \ref normalizedConvolveImage() for documentation.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void
        convolveImageWithMask(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                              MaskIterator mul, MaskAccessor am,
                              DestIterator dest_ul, DestAccessor dest_acc,
                              KernelIterator ki, KernelAccessor ak,
                              Diff2D kul, Diff2D klr, BorderTreatmentMode border);
    }
    \endcode

    \deprecatedAPI{convolveImageWithMask}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void
        convolveImageWithMask(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                              MaskIterator mul, MaskAccessor am,
                              DestIterator dest_ul, DestAccessor dest_acc,
                              KernelIterator ki, KernelAccessor ak,
                              Diff2D kul, Diff2D klr, BorderTreatmentMode border);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void convolveImageWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                   pair<MaskIterator, MaskAccessor> mask,
                                   pair<DestIterator, DestAccessor> dest,
                                   tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
                                   BorderTreatmentMode> kernel);
    }
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void convolveImageWithMask)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
inline void
convolveImageWithMask(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                      MaskIterator mul, MaskAccessor am,
                      DestIterator dest_ul, DestAccessor dest_acc,
                      KernelIterator ki, KernelAccessor ak,
                      Diff2D kul, Diff2D klr, BorderTreatmentMode border)
{
    normalizedConvolveImage(src_ul, src_lr, src_acc,
                            mul, am,
                            dest_ul, dest_acc,
                            ki, ak, kul, klr, border);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
inline
void convolveImageWithMask(
                           triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           pair<MaskIterator, MaskAccessor> mask,
                           pair<DestIterator, DestAccessor> dest,
                           tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
                           BorderTreatmentMode> kernel)
{
    normalizedConvolveImage(src.first, src.second, src.third,
                            mask.first, mask.second,
                            dest.first, dest.second,
                            kernel.first, kernel.second, kernel.third,
                            kernel.fourth, kernel.fifth);
}

//@}

/********************************************************/
/*                                                      */
/*                      Kernel2D                        */
/*                                                      */
/********************************************************/

/** \brief Generic 2 dimensional convolution kernel.

    This kernel may be used for convolution of 2 dimensional signals.

    Convolution functions access the kernel via an ImageIterator
    which they get by calling \ref center(). This iterator
    points to the center of the kernel. The kernel's size is given by its upperLeft()
    (upperLeft().x <= 0, upperLeft().y <= 0)
    and lowerRight() (lowerRight().x >= 0, lowerRight().y >= 0) methods.
    The desired border treatment mode is returned by borderTreatment().

    The different init functions create a kernel with the specified
    properties. The requirements for the kernel's value_type depend
    on the init function used. At least NumericTraits must be defined.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/stdconvolution.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), dest(w,h);
    ...

    // define horizontal Sobel filter
    vigra::Kernel2D<float> sobel;
    sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
                         0.125, 0.0, -0.125,
                         0.25,  0.0, -0.25,
                         0.125, 0.0, -0.125;

    convolveImage(src, dest, sobel);
    \endcode

    <b> Required Interface:</b>

    \code
    value_type v = NumericTraits<value_type>::one();
    \endcode

    See also the init functions.
*/
template <class ARITHTYPE = double>
class Kernel2D
{
public:
        /** the kernel's value type
         */
    typedef ARITHTYPE value_type;

        /** 2D random access iterator over the kernel's values
         */
    typedef typename BasicImage<value_type>::traverser Iterator;

        /** const 2D random access iterator over the kernel's values
         */
    typedef typename BasicImage<value_type>::const_traverser ConstIterator;

        /** the kernel's accessor
         */
    typedef typename BasicImage<value_type>::Accessor Accessor;

        /** the kernel's const accessor
         */
    typedef typename BasicImage<value_type>::ConstAccessor ConstAccessor;

    struct InitProxy
    {
        typedef typename
        BasicImage<value_type>::ScanOrderIterator Iterator;

        InitProxy(Iterator i, int count, value_type & norm)
            : iter_(i), base_(i),
              count_(count), sum_(count),
              norm_(norm)
        {}

        ~InitProxy()
#ifndef _MSC_VER
            throw(PreconditionViolation)
#elif _MSC_VER >= 1900
            noexcept(false)
#endif
        {
            vigra_precondition(count_ == 1 || count_ == sum_,
                               "Kernel2D::initExplicitly(): "
                               "Too few init values.");
        }

        InitProxy & operator,(value_type const & v)
        {
            if(count_ == sum_)  norm_ = *iter_;

            --count_;
            vigra_precondition(count_ > 0,
                               "Kernel2D::initExplicitly(): "
                               "Too many init values.");

            norm_ += v;

            ++iter_;
            *iter_ = v;

            return *this;
        }

        Iterator iter_, base_;
        int count_, sum_;
        value_type & norm_;
    };

    static value_type one() { return NumericTraits<value_type>::one(); }

        /** Default constructor.
            Creates a kernel of size 1x1 which would copy the signal
            unchanged.
        */
    Kernel2D()
        : kernel_(1, 1, one()),
          left_(0, 0),
          right_(0, 0),
          norm_(one()),
          border_treatment_(BORDER_TREATMENT_REFLECT)
    {}

        /** Copy constructor.
         */
    Kernel2D(Kernel2D const & k)
        : kernel_(k.kernel_),
          left_(k.left_),
          right_(k.right_),
          norm_(k.norm_),
          border_treatment_(k.border_treatment_)
    {}

        /** Copy assignment.
         */
    Kernel2D & operator=(Kernel2D const & k)
    {
        if(this != &k)
        {
        kernel_ = k.kernel_;
            left_ = k.left_;
            right_ = k.right_;
            norm_ = k.norm_;
        border_treatment_ = k.border_treatment_;
        }
        return *this;
    }

        /** Initialization.
            This initializes the kernel with the given constant. The norm becomes
            v*width()*height().

            Instead of a single value an initializer list of length width()*height()
            can be used like this:

            \code
            vigra::Kernel2D<float> binom;

            binom.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =
            0.0625, 0.125, 0.0625,
            0.125,  0.25,  0.125,
            0.0625, 0.125, 0.0625;
            \endcode

            In this case, the norm will be set to the sum of the init values.
            An initializer list of wrong length will result in a run-time error.
        */
    InitProxy operator=(value_type const & v)
    {
        int size = (right_.x - left_.x + 1) *
                   (right_.y - left_.y + 1);
        kernel_ = v;
        norm_ = (double)size*v;

        return InitProxy(kernel_.begin(), size, norm_);
    }

        /** Destructor.
         */
    ~Kernel2D()
    {}

        /** Init the 2D kernel as the cartesian product of two 1D kernels
            of type \ref Kernel1D. The norm becomes the product of the two original
            norms.

            <b> Required Interface:</b>

            The kernel's value_type must be a linear algebra.

            \code
            vigra::Kernel2D<...>::value_type v;
            v = v * v;
            \endcode
        */
    void initSeparable(Kernel1D<value_type> const & kx,
                             Kernel1D<value_type> const & ky)
    {
        left_ = Diff2D(kx.left(), ky.left());
        right_ = Diff2D(kx.right(), ky.right());
        int w = right_.x - left_.x + 1;
        int h = right_.y - left_.y + 1;
        kernel_.resize(w, h);

        norm_ = kx.norm() * ky.norm();

        typedef typename Kernel1D<value_type>::const_iterator KIter;
        typename Kernel1D<value_type>::Accessor ka;

        KIter kiy = ky.center() + left_.y;
        Iterator iy = center() + left_;

        for(int y=left_.y; y<=right_.y; ++y, ++kiy, ++iy.y)
        {
            KIter kix = kx.center() + left_.x;
            Iterator ix = iy;
            for(int x=left_.x; x<=right_.x; ++x, ++kix, ++ix.x)
            {
                *ix = ka(kix) * ka(kiy);
            }
        }
    }

        /** Init the 2D kernel as the cartesian product of two 1D kernels
            given explicitly by iterators and sizes. The norm becomes the
            sum of the resulting kernel values.

            <b> Required Interface:</b>

            The kernel's value_type must be a linear algebra.

            \code
            vigra::Kernel2D<...>::value_type v;
            v = v * v;
            v += v;
            \endcode

            <b> Preconditions:</b>

            \code
            xleft <= 0;
            xright >= 0;
            yleft <= 0;
            yright >= 0;
            \endcode
        */
    template <class KernelIterator>
    void initSeparable(KernelIterator kxcenter, int xleft, int xright,
                             KernelIterator kycenter, int yleft, int yright)
    {
        vigra_precondition(xleft <= 0 && yleft <= 0,
                           "Kernel2D::initSeparable(): left borders must be <= 0.");
        vigra_precondition(xright >= 0 && yright >= 0,
                           "Kernel2D::initSeparable(): right borders must be >= 0.");

        left_ = Point2D(xleft, yleft);
        right_ = Point2D(xright, yright);

        int w = right_.x - left_.x + 1;
        int h = right_.y - left_.y + 1;
        kernel_.resize(w, h);

        KernelIterator kiy = kycenter + left_.y;
        Iterator iy = center() + left_;

        for(int y=left_.y; y<=right_.y; ++y, ++kiy, ++iy.y)
        {
            KernelIterator kix = kxcenter + left_.x;
            Iterator ix = iy;
            for(int x=left_.x; x<=right_.x; ++x, ++kix, ++ix.x)
            {
                *ix = *kix * *kiy;
            }
        }

        typename BasicImage<value_type>::iterator i = kernel_.begin();
        typename BasicImage<value_type>::iterator iend = kernel_.end();
        norm_ = *i;
        ++i;

        for(; i!= iend; ++i)
        {
            norm_ += *i;
        }
    }

        /** Init as a 2D box filter with given radius.
         */
    void initAveraging(int radius)
    {
        Kernel1D<value_type> avg;
        avg.initAveraging(radius);
        return initSeparable(avg, avg);
    }

        /** Init as a 2D Gaussian function with given standard deviation and norm.
         */
    void initGaussian(double std_dev, value_type norm)
    {
        Kernel1D<value_type> gauss;
        gauss.initGaussian(std_dev, norm);
        return initSeparable(gauss, gauss);
    }

        /** Init as a 2D Gaussian function with given standard deviation and unit norm.
         */
    void initGaussian(double std_dev)
    {
        return initGaussian(std_dev, NumericTraits<value_type>::one());
    }

        /** Init the 2D kernel as a circular averaging filter. The norm will be
            calculated as
            <TT>NumericTraits<value_type>::one() / (number of non-zero kernel values)</TT>.
            The kernel's value_type must be a linear space.

            <b> Required Interface:</b>

            \code
            value_type v = vigra::NumericTraits<value_type>::one();

            double d;
            v = d * v;
            \endcode

            <b> Precondition:</b>

            \code
            radius > 0;
            \endcode
        */
    void initDisk(int radius)
    {
        vigra_precondition(radius > 0,
                           "Kernel2D::initDisk(): radius must be > 0.");

        left_ = Point2D(-radius, -radius);
        right_ = Point2D(radius, radius);
        int w = right_.x - left_.x + 1;
        int h = right_.y - left_.y + 1;
        kernel_.resize(w, h);
        norm_ = NumericTraits<value_type>::one();

        kernel_ = NumericTraits<value_type>::zero();
        double count = 0.0;

        Iterator k = center();
        double r2 = (double)radius*radius;

        int i;
        for(i=0; i<= radius; ++i)
        {
            double r = (double) i - 0.5;
            int w = (int)(VIGRA_CSTD::sqrt(r2 - r*r) + 0.5);
            for(int j=-w; j<=w; ++j)
            {
                k(j, i) = NumericTraits<value_type>::one();
                k(j, -i) = NumericTraits<value_type>::one();
                count += (i != 0) ? 2.0 : 1.0;
            }
        }

        count = 1.0 / count;

        for(int y=-radius; y<=radius; ++y)
        {
            for(int x=-radius; x<=radius; ++x)
            {
                k(x,y) = count * k(x,y);
            }
        }
    }

        /** Init the kernel by an explicit initializer list.
            The upper left and lower right corners (inclusive) of the kernel must be passed
            either as <tt>Shape2</tt> or <tt>Diff2D</tt> objects. A comma-separated initializer
            list for the kernel's weights is given after the assignment operator like this:

            \code
            // define horizontal Sobel filter
            vigra::Kernel2D<float> sobel;

            sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =
            0.125, 0.0, -0.125,
            0.25,  0.0, -0.25,
            0.125, 0.0, -0.125;
            \endcode

            The norm is set to the sum of the initializer values. If the wrong number of
            values is given, a run-time error results. It is, however, possible to give
            just one initializer. This creates an averaging filter with the given constant:

            \code
            vigra::Kernel2D<float> average3x3;

            average3x3.initExplicitly(Shape2(-1,-1), Shape2(1,1)) = 1.0/9.0;
            \endcode

            Here, the norm is set to value*width()*height().

            <b> Preconditions:</b>

            \code
            1. upperleft.x <= 0;
            2. upperleft.y <= 0;
            3. lowerright.x >= 0;
            4. lowerright.y >= 0;
            5. the number of values in the initializer list
            is 1 or equals the size of the kernel.
            \endcode
        */
    Kernel2D & initExplicitly(Shape2 const & upperleft, Shape2 const & lowerright)
    {
        vigra_precondition(upperleft[0] <= 0 && upperleft[1] <= 0,
                           "Kernel2D::initExplicitly(): left borders must be <= 0.");
        vigra_precondition(lowerright[0] >= 0 && lowerright[1] >= 0,
                           "Kernel2D::initExplicitly(): right borders must be >= 0.");

        left_ = Point2D(upperleft[0], upperleft[1]);
        right_ = Point2D(lowerright[0], lowerright[1]);

        int w = right_.x - left_.x + 1;
        int h = right_.y - left_.y + 1;
        kernel_.resize(w, h);

        return *this;
    }

    Kernel2D & initExplicitly(Diff2D const & upperleft, Diff2D const & lowerright)
    {
        return initExplicitly(Shape2(upperleft), Shape2(lowerright));
    }

        /** Coordinates of the upper left corner of the kernel.
         */
    Point2D upperLeft() const { return left_; }

        /** Coordinates of the lower right corner of the kernel.
         */
    Point2D lowerRight() const { return right_; }

        /** Width of the kernel.
         */
    int width() const { return right_.x - left_.x + 1; }

        /** Height of the kernel.
         */
    int height() const { return right_.y - left_.y + 1; }

        /** ImageIterator that points to the center of the kernel (coordinate (0,0)).
         */
    Iterator center() { return kernel_.upperLeft() - left_; }

        /** ImageIterator that points to the center of the kernel (coordinate (0,0)).
         */
    ConstIterator center() const { return kernel_.upperLeft() - left_; }

        /** Access kernel entry at given position.
         */
    value_type & operator()(int x, int y)
    { return kernel_[Diff2D(x,y) - left_]; }

        /** Read kernel entry at given position.
         */
    value_type operator()(int x, int y) const
    { return kernel_[Diff2D(x,y) - left_]; }

        /** Access kernel entry at given position.
         */
    value_type & operator[](Diff2D const & d)
    { return kernel_[d - left_]; }

        /** Read kernel entry at given position.
         */
    value_type operator[](Diff2D const & d) const
    { return kernel_[d - left_]; }

        /** Norm of the kernel (i.e. sum of its elements).
         */
    value_type norm() const { return norm_; }

        /** The kernels default accessor.
         */
    Accessor accessor() { return Accessor(); }

        /** The kernels default const accessor.
         */
    ConstAccessor accessor() const { return ConstAccessor(); }

        /** Normalize the kernel to the given value. (The norm is the sum of all kernel
            elements.) The kernel's value_type must be a division algebra or
            algebraic field.

            <b> Required Interface:</b>

            \code
            value_type v = vigra::NumericTraits<value_type>::one(); // if norm is not
                                                                    // given explicitly

            v += v;
            v = v * v;
            v = v / v;
            \endcode
        */
    void normalize(value_type norm)
    {
        typename BasicImage<value_type>::iterator i = kernel_.begin();
        typename BasicImage<value_type>::iterator iend = kernel_.end();
        typename NumericTraits<value_type>::RealPromote sum = *i;
        ++i;

        for(; i!= iend; ++i)
        {
            sum += *i;
        }

        sum = norm / sum;
        i = kernel_.begin();
        for(; i != iend; ++i)
        {
            *i = *i * sum;
        }

        norm_ = norm;
    }

        /** Normalize the kernel to norm 1.
         */
    void normalize()
    {
        normalize(one());
    }

        /** current border treatment mode
         */
    BorderTreatmentMode borderTreatment() const
    { return border_treatment_; }

        /** Set border treatment mode.
            Only <TT>BORDER_TREATMENT_CLIP</TT> and <TT>BORDER_TREATMENT_AVOID</TT> are currently
            allowed.
        */
    void setBorderTreatment( BorderTreatmentMode new_mode)
    {
        vigra_precondition((new_mode == BORDER_TREATMENT_CLIP    ||
                            new_mode == BORDER_TREATMENT_AVOID   ||
                            new_mode == BORDER_TREATMENT_REFLECT ||
                            new_mode == BORDER_TREATMENT_REPEAT  ||
                            new_mode == BORDER_TREATMENT_WRAP),
                           "convolveImage():\n"
                           "  Border treatment must be one of follow treatments:\n"
                           "  - BORDER_TREATMENT_CLIP\n"
                           "  - BORDER_TREATMENT_AVOID\n"
                           "  - BORDER_TREATMENT_REFLECT\n"
                           "  - BORDER_TREATMENT_REPEAT\n"
                           "  - BORDER_TREATMENT_WRAP\n");

        border_treatment_ = new_mode;
    }


private:
    BasicImage<value_type> kernel_;
    Point2D left_, right_;
    value_type norm_;
    BorderTreatmentMode border_treatment_;
};

/**************************************************************/
/*                                                            */
/*         Argument object factories for Kernel2D             */
/*                                                            */
/*     (documentation: see vigra/convolution.hxx)             */
/*                                                            */
/**************************************************************/

template <class KernelIterator, class KernelAccessor>
inline
tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(KernelIterator ik, KernelAccessor ak, Diff2D kul, Diff2D klr,
         BorderTreatmentMode border)

{
    return
        tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, BorderTreatmentMode> (
                                                             ik, ak, kul, klr, border);
}

template <class T>
inline
tuple5<typename Kernel2D<T>::ConstIterator,
       typename Kernel2D<T>::ConstAccessor,
       Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(Kernel2D<T> const & k)

{
    return
        tuple5<typename Kernel2D<T>::ConstIterator,
               typename Kernel2D<T>::ConstAccessor,
               Diff2D, Diff2D, BorderTreatmentMode>(
            k.center(),
            k.accessor(),
            k.upperLeft(), k.lowerRight(),
            k.borderTreatment());
}

template <class T>
inline
tuple5<typename Kernel2D<T>::ConstIterator,
       typename Kernel2D<T>::ConstAccessor,
       Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(Kernel2D<T> const & k, BorderTreatmentMode border)

{
    return
        tuple5<typename Kernel2D<T>::ConstIterator,
               typename Kernel2D<T>::ConstAccessor,
               Diff2D, Diff2D, BorderTreatmentMode>(
            k.center(),
            k.accessor(),
            k.upperLeft(), k.lowerRight(),
            border);
}


} // namespace vigra

#endif // VIGRA_STDCONVOLUTION_HXX
