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


#ifndef VIGRA_STDCONVOLUTION_HXX
#define VIGRA_STDCONVOLUTION_HXX

#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/bordertreatment.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/utilities.hxx"

namespace vigra {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void internalPixelEvaluationByClip(int x, int y, int w, int h, SrcIterator xs,
                                   SrcAccessor src_acc, DestIterator xd, DestAccessor dest_acc,
                                   KernelIterator ki, Diff2D kul, Diff2D klr, KernelAccessor ak)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SumType;
    typedef typename
        NumericTraits<typename KernelAccessor::value_type>::RealPromote KSumType;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    // calculate width and height of the kernel
    int kernel_width = klr.x - kul.x + 1;
    int kernel_height = klr.y - kul.y + 1;

    KSumType norm = NumericTraits<KSumType>::zero();
    SumType sum = NumericTraits<SumType>::zero();
    int xx, yy;

    //klr (kernel_lowerright) ist Diff2D !!!
    KernelIterator yk  = ki + klr;

    //Die Summe der Punkte im Kernel wird ermittelt (= norm)
    for(yy=0; yy<kernel_height; ++yy, --yk.y){
        KernelIterator xk  = yk;
        for(xx=0; xx<kernel_width; ++xx, --xk.x){
            norm += ak(xk);
        }
    }

    int x0, y0, x1, y1;

    y0 = (y<klr.y) ?  -y : -klr.y;
    y1 = (h-y-1<-kul.y) ? h-y-1 : -kul.y;

    x0 = (x<klr.x) ? -x : -klr.x;
    x1 = (w-x-1<-kul.x) ? w-x-1 : -kul.x;

    SrcIterator yys = xs + Diff2D(x0, y0);
    yk  = ki - Diff2D(x0, y0);

    KSumType ksum = NumericTraits<KSumType>::zero();
    kernel_width = x1 - x0 + 1;
    kernel_height = y1 - y0 + 1;

    //es wird zuerst abgeschnitten und dann gespigelt!

    for(yy=0; yy<kernel_height; ++yy, ++yys.y, --yk.y){
        SrcIterator xxs = yys;
        KernelIterator xk  = yk;

        for(xx=0; xx<kernel_width; ++xx, ++xxs.x, --xk.x){
            sum += ak(xk) * src_acc(xxs);
            ksum += ak(xk);
        }
    }

    //                      store average in destination pixel
    dest_acc.set(DestTraits::fromRealPromote((norm / ksum) * sum), xd);

}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void internalPixelEvaluationByWrapReflectRepeat(int x, int y, int src_width, int src_height, SrcIterator xs,
                                                SrcAccessor src_acc, DestIterator xd, DestAccessor dest_acc,
                                                KernelIterator ki, Diff2D kul, Diff2D klr, KernelAccessor ak,
                                                BorderTreatmentMode border)
{

    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SumType;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    SumType sum = NumericTraits<SumType>::zero();

    SrcIterator src_ul = xs - Diff2D(x, y);
    SrcIterator src_lr = src_ul + Diff2D(src_width, src_height);

    SrcIterator yys = xs;
    KernelIterator yk  = ki;

    // calculate width and height of the kernel
    int kernel_width = klr.x - kul.x + 1;
    int kernel_height = klr.y - kul.y + 1;

    //Zeigt an wo der Kernel über die Grenzen hinausgeht
    bool top_to_much = (y<klr.y) ? true : false;
    bool down_to_much = (src_height-y-1<-kul.y)? true : false;
    bool left_to_much = (x<klr.x)? true : false;
    bool right_to_much = (src_width-x-1<-kul.x)? true : false;

    //Die Richtung x und y !!!
    //in der bei der Iteration über das aktuelle Bereich im Bild
    //iteriert wird. Also wenn von ur->ll dann (-1, +1) und wenn lr->ul
    //dann (-1, -1).
    Diff2D way_increment;

    /* iteriert wird immer aus dem gültigen in den ungültigen
       Bereich! dieser Tupel setzt sich wie folgt zusammen:
       1. Wird bei der Iteration in X-Richtung ungültiger Bereich
       erreicht so wird mit border_increment.first gesprungen und
       mit border_increment.third weiter iteriert.
       2. Wird bei der Iteration in Y-Richtung ungültiger Bereich
       erreicht so wird mit border_increment.second gesprungen und
       mit border_increment.fourth weiter iteriert.
    */
    tuple4<int, int, int, int> border_increment;
    if (border == BORDER_TREATMENT_REPEAT){
        border_increment = tuple4<int, int, int, int>(1, 1, 0, 0);
    }else if (border == BORDER_TREATMENT_REFLECT){
        border_increment = tuple4<int, int, int, int>(2, 2, -1, -1);
    }else{ // BORDER_TREATMENT_WRAP
        border_increment = tuple4<int, int, int, int>(src_width, src_height, 1, 1);
    }

    pair<int, int> valid_step_count;

    if(left_to_much && !top_to_much && !down_to_much){
        yys += klr;
        yk += kul;
        way_increment = Diff2D(-1, -1);
        border_increment.third = -border_increment.third;
        border_increment.fourth = -border_increment.fourth;
        valid_step_count = make_pair((yys - src_ul).x + 1, kernel_height);
    }else if(top_to_much && !left_to_much && !right_to_much){
        yys += klr;
        yk += kul;
        way_increment = Diff2D(-1, -1);
        border_increment.third = -border_increment.third;
        border_increment.fourth = -border_increment.fourth;
        valid_step_count = make_pair(kernel_width, (yys - src_ul).y + 1);
    }else if(right_to_much && !top_to_much && !down_to_much){
        yys += kul;
        yk += klr;
        way_increment = Diff2D(1, 1);
        border_increment.first = -border_increment.first;
        border_increment.second = -border_increment.second;
        valid_step_count = make_pair((src_lr - yys).x, kernel_height);
    }else if(down_to_much && !left_to_much && !right_to_much){
        yys += kul;
        yk += klr;
        way_increment = Diff2D(1, 1);
        border_increment.first = -border_increment.first;
        border_increment.second = -border_increment.second;
        valid_step_count = make_pair(kernel_width, (src_lr - yys).y);
    }else if(down_to_much && left_to_much){
        yys += kul + Diff2D(kernel_width - 1, 0);
        yk += kul + Diff2D(0, kernel_height - 1);
        way_increment = Diff2D(-1, 1);
        border_increment.second = -border_increment.second;
        border_increment.third = -border_increment.third;
        valid_step_count = make_pair((yys - src_ul).x + 1, (src_lr - yys).y);
    }else if(down_to_much && right_to_much){
        yys += kul;
        yk += klr;
        way_increment = Diff2D(1, 1);
        border_increment.first = -border_increment.first;
        border_increment.second = -border_increment.second;
        valid_step_count = make_pair((src_lr - yys).x, (src_lr - yys).y);
    }else if(top_to_much && left_to_much){
        yys += klr;
        yk += kul;
        way_increment = Diff2D(-1, -1);
        border_increment.third = -border_increment.third;
        border_increment.fourth = -border_increment.fourth;
        valid_step_count = make_pair((yys - src_ul).x + 1, (yys - src_ul).y + 1);
    }else{ //top_to_much && right_to_much
        yys += kul + Diff2D(0, kernel_height - 1);
        yk += kul + Diff2D(kernel_width - 1, 0);
        way_increment = Diff2D(1, -1);
        border_increment.first = -border_increment.first;
        border_increment.fourth = -border_increment.fourth;
        valid_step_count = make_pair((src_lr - yys).x, (yys - src_ul).y + 1);
    }

    int yy = 0, xx;

    //laeuft den zulässigen Bereich in y-Richtung durch
    for(; yy < valid_step_count.second; ++yy, yys.y += way_increment.y, yk.y -= way_increment.y ){
        SrcIterator xxs = yys;
        KernelIterator xk  = yk;

        //laeuft den zulässigen Bereich in x-Richtung durch
        for(xx = 0; xx < valid_step_count.first; ++xx, xxs.x += way_increment.x, xk.x -= way_increment.x){
            sum += ak(xk) * src_acc(xxs);
        }

        //Nächstes ++xxs.x wuerde in unzulässigen Bereich
        //bringen => Sprung in zulaessigen Bereich
        xxs.x += border_increment.first;

        for( ; xx < kernel_width; ++xx, xxs.x += border_increment.third, xk.x -= way_increment.x ){
            sum += ak(xk) * src_acc(xxs);
        }
    }

    //Nächstes ++yys.y wuerde in unzulässigen Bereich
    //bringen => Sprung in zulaessigen Bereich
    yys.y += border_increment.second;

    for( ; yy < kernel_height; ++yy, yys.y += border_increment.third, yk.y -= way_increment.y){
        SrcIterator xxs = yys;
        KernelIterator xk  = yk;

        for(xx=0; xx < valid_step_count.first; ++xx, xxs.x += way_increment.x, xk.x -= way_increment.x){
            sum += ak(xk) * src_acc(xxs);
        }

        //Sprung in den zulaessigen Bereich
        xxs.x += border_increment.first;

        for( ; xx < kernel_width; ++xx, xxs.x += border_increment.third, xk.x -= way_increment.x ){
            sum += ak(xk) * src_acc(xxs);
        }
    }

    // store average in destination pixel
    dest_acc.set(DestTraits::fromRealPromote(sum), xd);

}// end of internalPixelEvaluationByWrapReflectRepeat


/** \addtogroup StandardConvolution Two-dimensional convolution functions

Perform 2D non-separable convolution, with and without ROI mask.

These generic convolution functions implement
the standard 2D convolution operation for images that fit
into the required interface. Arbitrary ROI's are supported
by the mask version of the algorithm.
The functions need a suitable 2D kernel to operate.
*/
//@{

/** \brief Performs a 2 dimensional convolution of the source image using the given
    kernel.

    The KernelIterator must point to the center of the kernel, and
    the kernel's size is given by its upper left (x and y of distance <= 0) and
    lower right (distance >= 0) corners. The image must always be larger than the
    kernel. At those positions where the kernel does not completely fit
    into the image, the specified \ref BorderTreatmentMode is
    applied. You can choice between following BorderTreatmentModes:
    <ul>
    <li>BORDER_TREATMENT_CLIP</li>
    <li>BORDER_TREATMENT_AVOID</li>
    <li>BORDER_TREATMENT_WRAP</li>
    <li>BORDER_TREATMENT_REFLECT</li>
    <li>BORDER_TREATMENT_REPEAT</li>
    </ul><br>
    The images's pixel type (SrcAccessor::value_type) must be a
    linear space over the kernel's value_type (KernelAccessor::value_type),
    i.e. addition of source values, multiplication with kernel values,
    and NumericTraits must be defined.
    The kernel's value_type must be an algebraic field,
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must
    be defined.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
    template <class SrcIterator, class SrcAccessor,
    class DestIterator, class DestAccessor,
    class KernelIterator, class KernelAccessor>
    void convolveImage(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
    DestIterator dest_ul, DestAccessor dest_acc,
    KernelIterator ki, KernelAccessor ak,
    Diff2D kul, Diff2D klr, BorderTreatmentMode border)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
    template <class SrcIterator, class SrcAccessor,
    class DestIterator, class DestAccessor,
    class KernelIterator, class KernelAccessor>
    void convolveImage(
    triple<SrcIterator, SrcIterator, SrcAccessor> src,
    pair<DestIterator, DestAccessor> dest,
    tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
    BorderTreatmentMode> kernel)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="stdconvolution_8hxx-source.html">vigra/stdconvolution.hxx</a>"<br>
    Namespace: vigra


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // define horizontal Sobel filter
    vigra::Kernel2D<float> sobel;

    sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
    0.125, 0.0, -0.125,
    0.25,  0.0, -0.25,
    0.125, 0.0, -0.125;

    vigra::convolveImage(srcImageRange(src), destImage(dest), kernel2d(sobel));
    \endcode

    <b> Required Interface:</b>

    \code
    ImageIterator src_ul, src_lr;
    ImageIterator dest_ul;
    ImageIterator ik;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    KernelAccessor kernel_accessor;

    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(src_ul);

    s = s + s;
    s = kernel_accessor(ik) * s;
    s -= s;

    dest_accessor.set(
    NumericTraits<DestAccessor::value_type>::fromRealPromote(s), dest_ul);

    NumericTraits<KernelAccessor::value_type>::RealPromote k = kernel_accessor(ik);

    k += k;
    k -= k;
    k = k / k;

    \endcode

    <b> Preconditions:</b>

    \code
    kul.x <= 0
    kul.y <= 0
    klr.x >= 0
    klr.y >= 0
    src_lr.x - src_ul.x >= klr.x + kul.x + 1
    src_lr.y - src_ul.y >= klr.y + kul.y + 1
    border == BORDER_TREATMENT_CLIP || border == BORDER_TREATMENT_AVOID
    \endcode

    If border == BORDER_TREATMENT_CLIP: Sum of kernel elements must be
    != 0.

*/
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
                        border == BORDER_TREATMENT_WRAP),
                       "convolveImage():\n"
                       "  Border treatment must be one of follow treatments:\n"
                       "  - BORDER_TREATMENT_CLIP\n"
                       "  - BORDER_TREATMENT_AVOID\n"
                       "  - BORDER_TREATMENT_REFLECT\n"
                       "  - BORDER_TREATMENT_REPEAT\n"
                       "  - BORDER_TREATMENT_WRAP\n");

    vigra_precondition(kul.x <= 0 && kul.y <= 0,
                       "convolveImage(): coordinates of "
                       "kernel's upper left must be <= 0.");
    vigra_precondition(klr.x >= 0 && klr.y >= 0,
                       "convolveImage(): coordinates of "
                       "kernel's lower right must be >= 0.");


    // use traits to determine SumType as to prevent possible overflow
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SumType;
    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;

    // calculate width and height of the image
    int w = src_lr.x - src_ul.x;
    int h = src_lr.y - src_ul.y;

    // calculate width and height of the kernel
    int kernel_width = klr.x - kul.x + 1;
    int kernel_height = klr.y - kul.y + 1;

    vigra_precondition(w >= kernel_width && h >= kernel_height,
                       "convolveImage(): kernel larger than image.");

    int x,y;

    // The start and endpoints of the image, that will be modified.
    // It's ever (0,0) and (w, h)
    // except by AVOID treatment mode.
    int ystart = (border == BORDER_TREATMENT_AVOID) ?  klr.y   : 0;
    int yend   = (border == BORDER_TREATMENT_AVOID) ?  h+kul.y : h;
    int xstart = (border == BORDER_TREATMENT_AVOID) ?  klr.x   : 0;
    int xend   = (border == BORDER_TREATMENT_AVOID) ?  w+kul.x : w;

    // create y iterators (Ausser AVOID bleibt alles bei *_ul)
    DestIterator yd = dest_ul + Diff2D(xstart, ystart);
    SrcIterator ys = src_ul + Diff2D(xstart, ystart);

    // Durchlauf über das ganze Bild (bei AVOID nur die "Mitte")
    for(y=ystart; y < yend; ++y, ++ys.y, ++yd.y){
        // create x iterators
        DestIterator xd(yd);
        SrcIterator xs(ys);

        for(x=xstart; x < xend; ++x, ++xs.x, ++xd.x){
            // init the sum
            SumType sum = NumericTraits<SumType>::zero();

            // how much of the kernel fits into the image ?
            bool nearBorder = false;

            nearBorder = (y<klr.y) || (h-y-1<-kul.y) || (x<klr.x) || (w-x-1<-kul.x);

            if(!nearBorder){
                SrcIterator yys = xs - klr;
                KernelIterator yk  = ki + klr;

                int xx, yy;
                for(yy=0; yy<kernel_height; ++yy, ++yys.y, --yk.y){
                    SrcIterator xxs = yys;
                    KernelIterator xk  = yk;

                    for(xx=0; xx<kernel_width; ++xx, ++xxs.x, --xk.x){
                        sum += ak(xk) * src_acc(xxs);
                    }
                }

                // store average in destination pixel
                dest_acc.set(DestTraits::fromRealPromote(sum), xd);
            }
            else{
                if (border == BORDER_TREATMENT_CLIP){

                    internalPixelEvaluationByClip(x, y, w, h, xs, src_acc, xd, dest_acc, ki, kul, klr, ak);

                }else{

                    internalPixelEvaluationByWrapReflectRepeat(x, y, w, h, xs, src_acc, xd, dest_acc, ki, kul, klr, ak, border);

                }
            }
        }
    }
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
inline
void convolveImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
                   BorderTreatmentMode> kernel)
{
    convolveImage(src.first, src.second, src.third,
                  dest.first, dest.second,
                  kernel.first, kernel.second, kernel.third,
                  kernel.fourth, kernel.fifth);
}


/** \brief Performs a 2 dimensional convolution of the source image within the
    given ROI mask using the given kernel.

    The ROI is applied as follows:
    Only pixel under the ROI are used in the calculations. Whenever a part of the
    kernel lies outside the ROI, the kernel is renormalized to its original
    norm (analogous to the CLIP \ref BorderTreatmentMode). An convolution result is
    calculated whenever at the current kernel position {\it at least one pixel of the
    kernel is within the ROI}. I.e., pixels not under the ROI may nevertheless
    be assigned a value if they are {\it near} the ROI. Thus, this algorithm is also
    useful as an interpolator. To get rid of the results outside the ROI mask, a
    subsequent \ref copyImageIf() must be performed.

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

    pass arguments explicitly:
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
    Diff2D kul, Diff2D klr, BorderTreatmentMode border)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
    template <class SrcIterator, class SrcAccessor,
    class MaskIterator, class MaskAccessor,
    class DestIterator, class DestAccessor,
    class KernelIterator, class KernelAccessor>
    inline
    void convolveImageWithMask(
    triple<SrcIterator, SrcIterator, SrcAccessor> src,
    pair<MaskIterator, MaskAccessor> mask,
    pair<DestIterator, DestAccessor> dest,
    tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
    BorderTreatmentMode> kernel)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="stdconvolution_8hxx-source.html">vigra/stdconvolution.hxx</a>"<br>
    Namespace: vigra


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

    vigra::convolveImage(srcImageRange(src), maskImage(mask), destImage(dest), kernel2d(binom));
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

    <b> Preconditions:</b>

    \code
    kul.x <= 0
    kul.y <= 0
    klr.x >= 0
    klr.y >= 0
    src_lr.x - src_ul.x >= klr.x + kul.x + 1
    src_lr.y - src_ul.y >= klr.y + kul.y + 1
    border == BORDER_TREATMENT_CLIP || border == BORDER_TREATMENT_AVOID
    \endcode

    Sum of kernel elements must be != 0.

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
void
convolveImageWithMask(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                      MaskIterator mul, MaskAccessor am,
                      DestIterator dest_ul, DestAccessor dest_acc,
                      KernelIterator ki, KernelAccessor ak,
                      Diff2D kul, Diff2D klr, BorderTreatmentMode border)
{
    vigra_precondition((border == BORDER_TREATMENT_CLIP  ||
                        border == BORDER_TREATMENT_AVOID),
                       "convolveImageWithMask(): "
                       "Border treatment must be BORDER_TREATMENT_CLIP or BORDER_TREATMENT_AVOID.");

    vigra_precondition(kul.x <= 0 && kul.y <= 0,
                       "convolveImageWithMask(): left borders must be <= 0.");
    vigra_precondition(klr.x >= 0 && klr.y >= 0,
                       "convolveImageWithMask(): right borders must be >= 0.");

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
                    SumType sum;
                    KSumType ksum;

                    SrcIterator yys = xs + Diff2D(x0, y0);
                    MaskIterator yym = xm + Diff2D(x0, y0);
                    KernelIterator yk  = ki - Diff2D(x0, y0);

                    int xx, yy, kernel_width, kernel_height;
                    kernel_width = x1 - x0 + 1;
                    kernel_height = y1 - y0 + 1;
                    for(yy=0; yy<kernel_height; ++yy, ++yys.y, --yk.y, ++yym.y)
                        {
                            SrcIterator xxs = yys;
                            MaskIterator xxm = yym;
                            KernelIterator xk  = yk;

                            for(xx=0; xx<kernel_width; ++xx, ++xxs.x, --xk.x, ++xxm.x)
                                {
                                    if(!am(xxm)) continue;

                                    if(first)
                                        {
                                            sum = ak(xk) * src_acc(xxs);
                                            ksum = ak(xk);
                                            first = false;
                                        }
                                    else
                                        {
                                            sum += ak(xk) * src_acc(xxs);
                                            ksum += ak(xk);
                                        }
                                }
                        }
                    // store average in destination pixel
                    if(!first &&
                       ksum != NumericTraits<KSumType>::zero())
                        {
                            dest_acc.set(DestTraits::fromRealPromote((norm / ksum) * sum), xd);
                        }
                }
        }
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
    convolveImageWithMask(src.first, src.second, src.third,
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
(Note that the \ref StandardConvolution "2D convolution functions" don't currently
support all modes.)

The different init functions create a kernel with the specified
properties. The requirements for the kernel's value_type depend
on the init function used. At least NumericTraits must be defined.

The kernel defines a factory function kernel2d() to create an argument object
(see \ref KernelArgumentObjectFactories).

<b> Usage:</b>

<b>\#include</b> "<a href="stdconvolution_8hxx-source.html">vigra/stdconvolution.hxx</a>"<br>
Namespace: vigra

\code
vigra::FImage src(w,h), dest(w,h);
...

// define horizontal Sobel filter
vigra::Kernel2D<float> sobel;

sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
0.125, 0.0, -0.125,
0.25,  0.0, -0.25,
0.125, 0.0, -0.125;

vigra::convolveImage(srcImageRange(src), destImage(dest), kernel2d(sobel));
\endcode

<b> Required Interface:</b>

\code
value_type v = NumericTraits<value_type>::one();
\endcode

See also the init functions.

*/

template <class ARITHTYPE>
class Kernel2D
{
public:
        /** the kernel's value type
         */
    typedef ARITHTYPE value_type;

        /** 2D random access iterator over the kernel's values
         */
    typedef typename BasicImage<value_type>::Iterator Iterator;

        /** the kernel's accessor
         */
    typedef typename BasicImage<value_type>::Accessor Accessor;

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
        : kernel_(1, 1, Kernel2D<ARITHTYPE>::one()),
          left_(0, 0),
          right_(0, 0),
          border_treatment_(BORDER_TREATMENT_CLIP)
    {}

        /** Copy constructor.
         */
    Kernel2D(Kernel2D const & k)
        : left_(k.left_),
          right_(k.right_),
          norm_(k.norm_),
          kernel_(k.kernel_),
          border_treatment_(k.border_treatment_)
    {}

        /** Copy assignment.
         */
    Kernel2D & operator=(Kernel2D const & k)
    {
        if(this != &k)
            {
                left_ = k.left_;
                right_ = k.right_;
                norm_ = k.norm_;
                kernel_ = k.kernel_;
            }
        return *this;
    }

        /** Initialisation.
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
    void initSeparable(Kernel1D<value_type> & kx,
                       Kernel1D<value_type> & ky)
    {
        left_ = Diff2D(kx.left(), ky.left());
        right_ = Diff2D(kx.right(), ky.right());
        int w = right_.x - left_.x + 1;
        int h = right_.y - left_.y + 1;
        kernel_.resize(w, h);

        norm_ = kx.norm() * ky.norm();

        typedef typename Kernel1D<value_type>::Iterator KIter;
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
            The upper left and lower right corners of the kernel must be passed.
            A comma-separated initializer list is given after the assignment operator.
            This function is used like this:

            \code
            // define horizontal Sobel filter
            vigra::Kernel2D<float> sobel;

            sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =
            0.125, 0.0, -0.125,
            0.25,  0.0, -0.25,
            0.125, 0.0, -0.125;
            \endcode

            The norm is set to the sum of the initialzer values. If the wrong number of
            values is given, a run-time error results. It is, however, possible to give
            just one initializer. This creates an averaging filter with the given constant:

            \code
            vigra::Kernel2D<float> average3x3;

            average3x3.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 1.0/9.0;
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
    Kernel2D & initExplicitly(Diff2D upperleft, Diff2D lowerright)
    {
        vigra_precondition(upperleft.x <= 0 && upperleft.y <= 0,
                           "Kernel2D::initExplicitly(): left borders must be <= 0.");
        vigra_precondition(lowerright.x >= 0 && lowerright.y >= 0,
                           "Kernel2D::initExplicitly(): right borders must be >= 0.");

        left_ = Point2D(upperleft);
        right_ = Point2D(lowerright);

        int w = right_.x - left_.x + 1;
        int h = right_.y - left_.y + 1;
        kernel_.resize(w, h);

        return *this;
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

        /** Access kernel entry at given position.
         */
    value_type & operator()(int x, int y)
    { return kernel_[Diff2D(x,y) - left_]; }

        /** Read kernel entry at given position.
         */
    value_type operator()(int x, int y) const
    { return kernel_[Diff2D(x,y) - left_]; }

        /** Norm of the kernel (i.e. sum of its elements).
         */
    value_type norm() const { return norm_; }

        /** The kernels default accessor.
         */
    Accessor accessor() const { return Accessor(); }

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
tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor,
       Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(Kernel2D<T> & k)

{
    return
        tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor,
        Diff2D, Diff2D, BorderTreatmentMode>(
            k.center(),
            k.accessor(),
            k.upperLeft(), k.lowerRight(),
            k.borderTreatment());
}

template <class T>
inline
tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor,
       Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(Kernel2D<T> & k, BorderTreatmentMode border)

{
    return
        tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor,
        Diff2D, Diff2D, BorderTreatmentMode>(
            k.center(),
            k.accessor(),
            k.upperLeft(), k.lowerRight(),
            border);
}


} // namespace vigra

#endif // VIGRA_STDCONVOLUTION_HXX
