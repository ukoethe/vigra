/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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


#ifndef VIGRA_SEPARABLECONVOLUTION_HXX
#define VIGRA_SEPARABLECONVOLUTION_HXX

#include <cmath>
#include "utilities.hxx"
#include "numerictraits.hxx"
#include "imageiteratoradapter.hxx"
#include "bordertreatment.hxx"
#include "gaussians.hxx"
#include "array_vector.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                internalConvolveLineWrap              */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void internalConvolveLineWrap(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                              DestIterator id, DestAccessor da,
                              KernelIterator kernel, KernelAccessor ka,
                              int kleft, int kright)
{
  //    int w = iend - is;
    int w = std::distance( is, iend );

    typedef typename NumericTraits<typename
                      SrcAccessor::value_type>::RealPromote SumType;

    SrcIterator ibegin = is;

    for(int x=0; x<w; ++x, ++is, ++id)
    {
        KernelIterator ik = kernel + kright;
        SumType sum = NumericTraits<SumType>::zero();

        if(x < kright)
        {
            int x0 = x - kright;
            SrcIterator iss = iend + x0;

            for(; x0; ++x0, --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }

            iss = ibegin;
            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }
        else if(w-x <= -kleft)
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = iend;
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }

            int x0 = -kleft - w + x + 1;
            iss = ibegin;

            for(; x0; --x0, --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }
        else
        {
            SrcIterator iss = is - kright;
            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }

        da.set(NumericTraits<typename
                      DestAccessor::value_type>::fromRealPromote(sum), id);
    }
}

/********************************************************/
/*                                                      */
/*                internalConvolveLineClip              */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor,
          class Norm>
void internalConvolveLineClip(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                              DestIterator id, DestAccessor da,
                              KernelIterator kernel, KernelAccessor ka,
                              int kleft, int kright, Norm norm)
{
  //    int w = iend - is;
    int w = std::distance( is, iend );

    typedef typename NumericTraits<typename
                      SrcAccessor::value_type>::RealPromote SumType;

    SrcIterator ibegin = is;

    for(int x=0; x<w; ++x, ++is, ++id)
    {
        KernelIterator ik = kernel + kright;
        SumType sum = NumericTraits<SumType>::zero();

        if(x < kright)
        {
            int x0 = x - kright;
            Norm clipped = NumericTraits<Norm>::zero();

            for(; x0; ++x0, --ik)
            {
                clipped += ka(ik);
            }

            SrcIterator iss = ibegin;
            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }

            sum = norm / (norm - clipped) * sum;
        }
        else if(w-x <= -kleft)
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = iend;
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }

            Norm clipped = NumericTraits<Norm>::zero();

            int x0 = -kleft - w + x + 1;

            for(; x0; --x0, --ik)
            {
                clipped += ka(ik);
            }

            sum = norm / (norm - clipped) * sum;
        }
        else
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }

        da.set(NumericTraits<typename
                      DestAccessor::value_type>::fromRealPromote(sum), id);
    }
}

/********************************************************/
/*                                                      */
/*             internalConvolveLineReflect              */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void internalConvolveLineReflect(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                              DestIterator id, DestAccessor da,
                              KernelIterator kernel, KernelAccessor ka,
                              int kleft, int kright)
{
  //    int w = iend - is;
    int w = std::distance( is, iend );

    typedef typename NumericTraits<typename
                      SrcAccessor::value_type>::RealPromote SumType;

    SrcIterator ibegin = is;

    for(int x=0; x<w; ++x, ++is, ++id)
    {
        KernelIterator ik = kernel + kright;
        SumType sum = NumericTraits<SumType>::zero();

        if(x < kright)
        {
            int x0 = x - kright;
            SrcIterator iss = ibegin - x0;

            for(; x0; ++x0, --ik, --iss)
            {
                sum += ka(ik) * sa(iss);
            }

            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }
        else if(w-x <= -kleft)
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = iend;
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }

            int x0 = -kleft - w + x + 1;
            iss = iend - 2;

            for(; x0; --x0, --ik, --iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }
        else
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }

        da.set(NumericTraits<typename
                      DestAccessor::value_type>::fromRealPromote(sum), id);
    }
}

/********************************************************/
/*                                                      */
/*             internalConvolveLineRepeat               */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void internalConvolveLineRepeat(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                              DestIterator id, DestAccessor da,
                              KernelIterator kernel, KernelAccessor ka,
                              int kleft, int kright)
{
  //    int w = iend - is;
    int w = std::distance( is, iend );

    typedef typename NumericTraits<typename
                      SrcAccessor::value_type>::RealPromote SumType;

    SrcIterator ibegin = is;

    for(int x=0; x<w; ++x, ++is, ++id)
    {
        KernelIterator ik = kernel + kright;
        SumType sum = NumericTraits<SumType>::zero();

        if(x < kright)
        {
            int x0 = x - kright;
            SrcIterator iss = ibegin;

            for(; x0; ++x0, --ik)
            {
                sum += ka(ik) * sa(iss);
            }

            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }
        else if(w-x <= -kleft)
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = iend;
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }

            int x0 = -kleft - w + x + 1;
            iss = iend - 1;

            for(; x0; --x0, --ik)
            {
                sum += ka(ik) * sa(iss);
            }
        }
        else
        {
            SrcIterator iss = is + (-kright);
            SrcIterator isend = is + (1 - kleft);
            for(; iss != isend ; --ik, ++iss)
            {
                sum += ka(ik) * sa(iss);
            }
        }

        da.set(NumericTraits<typename
                      DestAccessor::value_type>::fromRealPromote(sum), id);
    }
}

/********************************************************/
/*                                                      */
/*              internalConvolveLineAvoid               */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void internalConvolveLineAvoid(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                              DestIterator id, DestAccessor da,
                              KernelIterator kernel, KernelAccessor ka,
                              int kleft, int kright)
{
  //    int w = iend - is;
    int w = std::distance( is, iend );

    typedef typename NumericTraits<typename
                      SrcAccessor::value_type>::RealPromote SumType;

    is += kright;
    id += kright;

    for(int x=kright; x<w+kleft; ++x, ++is, ++id)
    {
        KernelIterator ik = kernel + kright;
        SumType sum = NumericTraits<SumType>::zero();

        SrcIterator iss = is + (-kright);
        SrcIterator isend = is + (1 - kleft);
        for(; iss != isend ; --ik, ++iss)
        {
            sum += ka(ik) * sa(iss);
        }

        da.set(NumericTraits<typename
                      DestAccessor::value_type>::fromRealPromote(sum), id);
    }
}

/********************************************************/
/*                                                      */
/*         Separable convolution functions              */
/*                                                      */
/********************************************************/

/** \addtogroup SeparableConvolution One-dimensional and separable convolution functions

    Perform 1D convolution and separable filtering in 2 dimensions.

    These generic convolution functions implement
    the standard convolution operation for a wide range of images and
    signals that fit into the required interface. They need a suitable
    kernel to operate.
*/
//@{

/** \brief Performs a 1 dimensional convolution of the source signal using the given
    kernel.

    The KernelIterator must point to the center iterator, and
    the kernel's size is given by its left (kleft <= 0) and right
    (kright >= 0) borders. The signal must always be larger than the kernel.
    At those positions where the kernel does not completely fit
    into the signal's range, the specified \ref BorderTreatmentMode is
    applied.

    The signal's value_type (SrcAccessor::value_type) must be a
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
        void convolveLine(SrcIterator is, SrcIterator isend, SrcAccessor sa,
                          DestIterator id, DestAccessor da,
                          KernelIterator ik, KernelAccessor ka,
                          int kleft, int kright, BorderTreatmentMode border)
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void convolveLine(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          tuple5<KernelIterator, KernelAccessor,
                                 int, int, BorderTreatmentMode> kernel)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="separableconvolution_8hxx-source.html">vigra/separableconvolution.hxx</a>"


    \code
    std::vector<float> src, dest;
    ...

    // define binomial filter of size 5
    static float kernel[] =
           { 1.0/16.0, 4.0/16.0, 6.0/16.0, 4.0/16.0, 1.0/16.0};

    typedef vigra::StandardAccessor<float> FAccessor;
    typedef vigra::StandardAccessor<float> KernelAccessor;


    vigra::convolveLine(src.begin(), src.end(), FAccessor(), dest.begin(), FAccessor(),
             kernel+2, KernelAccessor(), -2, 2, BORDER_TREATMENT_REFLECT);
    //       ^^^^^^^^  this is the center of the kernel

    \endcode

    <b> Required Interface:</b>

    \code
    RandomAccessIterator is, isend;
    RandomAccessIterator id;
    RandomAccessIterator ik;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    KernelAccessor kernel_accessor;

    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(is);

    s = s + s;
    s = kernel_accessor(ik) * s;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), id);

    \endcode

    If border == BORDER_TREATMENT_CLIP:

    \code
    NumericTraits<KernelAccessor::value_type>::RealPromote k = kernel_accessor(ik);

    k = k + k;
    k = k - k;
    k = k * k;
    k = k / k;

    \endcode

    <b> Preconditions:</b>

    \code
    kleft <= 0
    kright >= 0
    iend - is >= kright + kleft + 1
    \endcode

    If border == BORDER_TREATMENT_CLIP: Sum of kernel elements must be
    != 0.

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void convolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  KernelIterator ik, KernelAccessor ka,
                  int kleft, int kright, BorderTreatmentMode border)
{
    typedef typename KernelAccessor::value_type KernelValue;

    vigra_precondition(kleft <= 0,
                 "convolveLine(): kleft must be <= 0.\n");
    vigra_precondition(kright >= 0,
                 "convolveLine(): kright must be >= 0.\n");

    //    int w = iend - is;
    int w = std::distance( is, iend );

    vigra_precondition(w >= kright - kleft + 1,
                 "convolveLine(): kernel longer than line\n");

    switch(border)
    {
      case BORDER_TREATMENT_WRAP:
      {
        internalConvolveLineWrap(is, iend, sa, id, da, ik, ka, kleft, kright);
        break;
      }
      case BORDER_TREATMENT_AVOID:
      {
        internalConvolveLineAvoid(is, iend, sa, id, da, ik, ka, kleft, kright);
        break;
      }
      case BORDER_TREATMENT_REFLECT:
      {
        internalConvolveLineReflect(is, iend, sa, id, da, ik, ka, kleft, kright);
        break;
      }
      case BORDER_TREATMENT_REPEAT:
      {
        internalConvolveLineRepeat(is, iend, sa, id, da, ik, ka, kleft, kright);
        break;
      }
      case BORDER_TREATMENT_CLIP:
      {
        // find norm of kernel
        typedef typename KernelAccessor::value_type KT;
        KT norm = NumericTraits<KT>::zero();
        KernelIterator iik = ik + kleft;
        for(int i=kleft; i<=kright; ++i, ++iik) norm += ka(iik);

        vigra_precondition(norm != NumericTraits<KT>::zero(),
                     "convolveLine(): Norm of kernel must be != 0"
                     " in mode BORDER_TREATMENT_CLIP.\n");

        internalConvolveLineClip(is, iend, sa, id, da, ik, ka, kleft, kright, norm);
        break;
      }
      default:
      {
        vigra_precondition(0,
                     "convolveLine(): Unknown border treatment mode.\n");
      }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
inline
void convolveLine(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  tuple5<KernelIterator, KernelAccessor,
                         int, int, BorderTreatmentMode> kernel)
{
    convolveLine(src.first, src.second, src.third,
                 dest.first, dest.second,
                 kernel.first, kernel.second,
                 kernel.third, kernel.fourth, kernel.fifth);
}

/********************************************************/
/*                                                      */
/*                      separableConvolveX              */
/*                                                      */
/********************************************************/

/** \brief Performs a 1 dimensional convolution in x direction.

    It calls \link SeparableConvolution#convolveLine convolveLine\endlink() for every row of the
    image. See \link SeparableConvolution#convolveLine convolveLine\endlink() for more information about required interfaces
    and vigra_preconditions.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void separableConvolveX(SrcImageIterator supperleft,
                                SrcImageIterator slowerright, SrcAccessor sa,
                                DestImageIterator dupperleft, DestAccessor da,
                                KernelIterator ik, KernelAccessor ka,
                                int kleft, int kright, BorderTreatmentMode border)
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void separableConvolveX(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                pair<DestImageIterator, DestAccessor> dest,
                                tuple5<KernelIterator, KernelAccessor,
                                             int, int, BorderTreatmentMode> kernel)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="separableconvolution_8hxx-source.html">vigra/separableconvolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // define Gaussian kernel with std. deviation 3.0
    vigra::Kernel1D<double> kernel;
    kernel.initGaussian(3.0);

    vigra::separableConvolveX(srcImageRange(src), destImage(dest), kernel1d(kernel));

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void separableConvolveX(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                        KernelIterator ik, KernelAccessor ka,
                        int kleft, int kright, BorderTreatmentMode border)
{
    typedef typename KernelAccessor::value_type KernelValue;

    vigra_precondition(kleft <= 0,
                 "separableConvolveX(): kleft must be <= 0.\n");
    vigra_precondition(kright >= 0,
                 "separableConvolveX(): kright must be >= 0.\n");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;

    vigra_precondition(w >= kright - kleft + 1,
                 "separableConvolveX(): kernel longer than line\n");

    int y;

    for(y=0; y<h; ++y, ++supperleft.y, ++dupperleft.y)
    {
        typename SrcIterator::row_iterator rs = supperleft.rowIterator();
        typename DestIterator::row_iterator rd = dupperleft.rowIterator();

        convolveLine(rs, rs+w, sa, rd, da,
                     ik, ka, kleft, kright, border);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
inline void
separableConvolveX(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  tuple5<KernelIterator, KernelAccessor,
                         int, int, BorderTreatmentMode> kernel)
{
    separableConvolveX(src.first, src.second, src.third,
                 dest.first, dest.second,
                 kernel.first, kernel.second,
                 kernel.third, kernel.fourth, kernel.fifth);
}



/********************************************************/
/*                                                      */
/*                      separableConvolveY              */
/*                                                      */
/********************************************************/

/** \brief Performs a 1 dimensional convolution in y direction.

    It calls \link SeparableConvolution#convolveLine convolveLine\endlink() for every column of the
    image. See \link SeparableConvolution#convolveLine convolveLine\endlink() for more information about required interfaces
    and vigra_preconditions.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void separableConvolveY(SrcImageIterator supperleft,
                                SrcImageIterator slowerright, SrcAccessor sa,
                                DestImageIterator dupperleft, DestAccessor da,
                                KernelIterator ik, KernelAccessor ka,
                                int kleft, int kright, BorderTreatmentMode border)
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void separableConvolveY(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                pair<DestImageIterator, DestAccessor> dest,
                                tuple5<KernelIterator, KernelAccessor,
                                             int, int, BorderTreatmentMode> kernel)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="separableconvolution_8hxx-source.html">vigra/separableconvolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // define Gaussian kernel with std. deviation 3.0
    vigra::Kernel1D kernel;
    kernel.initGaussian(3.0);

    vigra::separableConvolveY(srcImageRange(src), destImage(dest), kernel1d(kernel));

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void separableConvolveY(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                        KernelIterator ik, KernelAccessor ka,
                        int kleft, int kright, BorderTreatmentMode border)
{
    typedef typename KernelAccessor::value_type KernelValue;

    vigra_precondition(kleft <= 0,
                 "separableConvolveY(): kleft must be <= 0.\n");
    vigra_precondition(kright >= 0,
                 "separableConvolveY(): kright must be >= 0.\n");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;

    vigra_precondition(h >= kright - kleft + 1,
                 "separableConvolveY(): kernel longer than line\n");

    int x;

    for(x=0; x<w; ++x, ++supperleft.x, ++dupperleft.x)
    {
        typename SrcIterator::column_iterator cs = supperleft.columnIterator();
        typename DestIterator::column_iterator cd = dupperleft.columnIterator();

        convolveLine(cs, cs+h, sa, cd, da,
                     ik, ka, kleft, kright, border);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
inline void
separableConvolveY(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  tuple5<KernelIterator, KernelAccessor,
                         int, int, BorderTreatmentMode> kernel)
{
    separableConvolveY(src.first, src.second, src.third,
                 dest.first, dest.second,
                 kernel.first, kernel.second,
                 kernel.third, kernel.fourth, kernel.fifth);
}

//@}

/********************************************************/
/*                                                      */
/*                      Kernel1D                        */
/*                                                      */
/********************************************************/

/** \brief Generic 1 dimensional convolution kernel.

    This kernel may be used for convolution of 1 dimensional signals or for
    separable convolution of multidimensional signals.

    Convlution functions access the kernel via a 1 dimensional random access
    iterator which they get by calling \ref center(). This iterator
    points to the center of the kernel. The kernel's size is given by its left() (<=0)
    and right() (>= 0) methods. The desired border treatment mode is
    returned by borderTreatment().

    The different init functions create a kernel with the specified
    properties. The kernel's value_type must be a linear space, i.e. it
    must define multiplication with doubles and NumericTraits.


    The kernel defines a factory function kernel1d() to create an argument object
    (see \ref KernelArgumentObjectFactories).

    <b> Usage:</b>

    <b>\#include</b> "<a href="stdconvolution_8hxx-source.html">vigra/stdconvolution.hxx</a>"

    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // define Gaussian kernel with std. deviation 3.0
    vigra::Kernel1D kernel;
    kernel.initGaussian(3.0);

    vigra::separableConvolveX(srcImageRange(src), destImage(dest), kernel1d(kernel));
    \endcode

    <b> Required Interface:</b>

    \code
    value_type v = vigra::NumericTraits<value_type>::one(); // if norm is not
                                                            // given explicitly
    double d;

    v = d * v;
    \endcode
*/

template <class ARITHTYPE>
class Kernel1D
{
  public:
    typedef ArrayVector<ARITHTYPE> InternalVector;

        /** the kernel's value type
        */
    typedef typename InternalVector::value_type value_type;

        /** the kernel's reference type
        */
    typedef typename InternalVector::reference reference;

        /** the kernel's const reference type
        */
    typedef typename InternalVector::const_reference const_reference;

        /** deprecated -- use Kernel1D::iterator
        */
    typedef typename InternalVector::iterator Iterator;

        /** 1D random access iterator over the kernel's values
        */
    typedef typename InternalVector::iterator iterator;

        /** const 1D random access iterator over the kernel's values
        */
    typedef typename InternalVector::const_iterator const_iterator;

        /** the kernel's accessor
        */
    typedef StandardAccessor<ARITHTYPE> Accessor;

        /** the kernel's const accessor
        */
    typedef StandardConstAccessor<ARITHTYPE> ConstAccessor;

    struct InitProxy
    {
        InitProxy(Iterator i, int count, value_type & norm)
        : iter_(i), base_(i),
          count_(count), sum_(count),
          norm_(norm)
        {}

        ~InitProxy()
        {
            vigra_precondition(count_ == 1 || count_ == sum_,
                  "Kernel1D::initExplicitly(): "
                  "Wrong number of init values.");
        }

        InitProxy & operator,(value_type const & v)
        {
            if(sum_ == count_) 
                norm_ = *iter_;

            norm_ += v;

            --count_;

            if(count_ > 0)
            {
                ++iter_;
                *iter_ = v;
            }

            return *this;
        }

        Iterator iter_, base_;
        int count_, sum_;
        value_type & norm_;
    };

    static value_type one() { return NumericTraits<value_type>::one(); }

        /** Default constructor.
            Creates a kernel of size 1 which would copy the signal
            unchanged.
        */
    Kernel1D()
    : kernel_(),
      left_(0),
      right_(0),
      border_treatment_(BORDER_TREATMENT_REFLECT),
      norm_(one())
    {
        kernel_.push_back(norm_);
    }

        /** Copy constructor.
        */
    Kernel1D(Kernel1D const & k)
    : kernel_(k.kernel_),
      left_(k.left_),
      right_(k.right_),
      border_treatment_(k.border_treatment_),
      norm_(k.norm_)
    {}

        /** Copy assignment.
        */
    Kernel1D & operator=(Kernel1D const & k)
    {
        if(this != &k)
        {
            left_ = k.left_;
            right_ = k.right_;
            border_treatment_ = k.border_treatment_;
            norm_ = k.norm_;
            kernel_ = k.kernel_;
        }
        return *this;
    }

        /** Initialization.
            This initializes the kernel with the given constant. The norm becomes
            v*size().

            Instead of a single value an initializer list of length size()
            can be used like this:

            \code
            vigra::Kernel1D<float> roberts_gradient_x;

            roberts_gradient_x.initExplicitly(0, 1) = 1.0, -1.0;
            \endcode

            In this case, the norm will be set to the sum of the init values.
            An initializer list of wrong length will result in a run-time error.
        */
    InitProxy operator=(value_type const & v)
    {
        int size = right_ - left_ + 1;
        for(unsigned int i=0; i<kernel_.size(); ++i) kernel_[i] = v;
        norm_ = (double)size*v;

        return InitProxy(kernel_.begin(), size, norm_);
    }

        /** Destructor.
        */
    ~Kernel1D()
    {}

        /**
            Init as a sampled Gaussian function. The radius of the kernel is
            always 3*std_dev. '<tt>norm</tt>' denotes the sum of all bins of the kernel
            (i.e. the kernel is corrected for the normalization error introduced
             by windowing the Gaussian to a finite interval). However,
            if <tt>norm</tt> is 0.0, the kernel is normalized to 1 by the analytic
            expression for the Gaussian, and <b>no</b> correction for the windowing
            error is performed.

            Precondition:
            \code
            std_dev >= 0.0
            \endcode

            Postconditions:
            \code
            1. left()  == -(int)(3.0*std_dev + 0.5)
            2. right() ==  (int)(3.0*std_dev + 0.5)
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == norm
            \endcode
        */
    void initGaussian(double std_dev, value_type norm);

        /** Init as a Gaussian function with norm 1.
         */
    void initGaussian(double std_dev)
    {
        initGaussian(std_dev, one());
    }


        /**
            Init as Lindeberg's discrete analog of the Gaussian function. The radius of the kernel is
            always 3*std_dev. 'norm' denotes the sum of all bins of the kernel.

            Precondition:
            \code
            std_dev >= 0.0
            \endcode

            Postconditions:
            \code
            1. left()  == -(int)(3.0*std_dev + 0.5)
            2. right() ==  (int)(3.0*std_dev + 0.5)
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == norm
            \endcode
        */
    void initDiscreteGaussian(double std_dev, value_type norm);

        /** Init as a LOineberg's discrete analog of the Gaussian function
            with norm 1.
         */
    void initDiscreteGaussian(double std_dev)
    {
        initDiscreteGaussian(std_dev, one());
    }

        /**
            Init as a Gaussian derivative of order '<tt>order</tt>'.
            The radius of the kernel is always <tt>3*std_dev + 0.5*order</tt>.
            '<tt>norm</tt>' denotes the norm of the kernel so that the
            following condition is fulfilled:

            \f[ \sum_{i=left()}^{right()}
                         \frac{(-i)^{order}kernel[i]}{order!} = norm
            \f]

            Thus, the kernel will be corrected for the error introduced
            by windowing the Gaussian to a finite interval. However,
            if <tt>norm</tt> is 0.0, the kernel is normalized to 1 by the analytic
            expression for the Gaussian derivative, and <b>no</b> correction for the
            windowing error is performed.

            Preconditions:
            \code
            1. std_dev >= 0.0
            2. order   >= 1
            \endcode

            Postconditions:
            \code
            1. left()  == -(int)(3.0*std_dev + 0.5*order + 0.5)
            2. right() ==  (int)(3.0*std_dev + 0.5*order + 0.5)
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == norm
            \endcode
        */
    void initGaussianDerivative(double std_dev, int order, value_type norm);

        /** Init as a Gaussian derivative with norm 1.
         */
    void initGaussianDerivative(double std_dev, int order)
    {
        initGaussianDerivative(std_dev, order, one());
    }

        /**
            Init an optimal 3-tap smoothing filter.
            The filter values are 
            
            \code
            [0.216, 0.568, 0.216]
            \endcode
            
            These values are optimal in the sense that the 3x3 filter obtained by separable application
            of this filter is the best possible 3x3 approximation to a Gaussian filter.
            The equivalent Gaussian has sigma = 0.680.
 
            Postconditions:
            \code
            1. left()  == -1
            2. right() ==  1
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalSmoothing3()
    {
        this->initExplicitly(-1, 1) = 0.216, 0.568, 0.216;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init an optimal 3-tap smoothing filter to be used in the context of first derivative computation.
            This filter must be used in conjunction with the symmetric difference filter (see initSymmetricDifference()), 
            such that the difference filter is applied along one dimension, and the smoothing filter along the other.
            The filter values are 
            
            \code
            [0.224365, 0.55127, 0.224365]
            \endcode
            
            These values are optimal in the sense that the 3x3 filter obtained by combining 
            this filter with the symmetric difference is the best possible 3x3 approximation to a 
            Gaussian first derivative filter. The equivalent Gaussian has sigma = 0.675.
 
            Postconditions:
            \code
            1. left()  == -1
            2. right() ==  1
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalFirstDerivativeSmoothing3()
    {
        this->initExplicitly(-1, 1) = 0.224365, 0.55127, 0.224365;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init an optimal 3-tap smoothing filter to be used in the context of second derivative computation.
            This filter must be used in conjunction with the 3-tap second difference filter (see initSecondDifference3()), 
            such that the difference filter is applied along one dimension, and the smoothing filter along the other.
            The filter values are 
            
            \code
            [0.13, 0.74, 0.13]
            \endcode
            
            These values are optimal in the sense that the 3x3 filter obtained by combining 
            this filter with the 3-tap second difference is the best possible 3x3 approximation to a 
            Gaussian second derivative filter. The equivalent Gaussian has sigma = 0.433.
 
            Postconditions:
            \code
            1. left()  == -1
            2. right() ==  1
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalSecondDerivativeSmoothing3()
    {
        this->initExplicitly(-1, 1) = 0.13, 0.74, 0.13;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init an optimal 5-tap smoothing filter.
            The filter values are 
            
            \code
            [0.03134, 0.24, 0.45732, 0.24, 0.03134]
            \endcode
            
            These values are optimal in the sense that the 5x5 filter obtained by separable application
            of this filter is the best possible 5x5 approximation to a Gaussian filter.
            The equivalent Gaussian has sigma = 0.867.
 
            Postconditions:
            \code
            1. left()  == -2
            2. right() ==  2
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalSmoothing5()
    {
        this->initExplicitly(-2, 2) = 0.03134, 0.24, 0.45732, 0.24, 0.03134;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init an optimal 5-tap smoothing filter to be used in the context of first derivative computation.
           This filter must be used in conjunction with the optimal 5-tap first derivative filter 
           (see initOptimalFirstDerivative5()),  such that the derivative filter is applied along one dimension, 
           and the smoothing filter along the other. The filter values are 
            
            \code
            [0.04255, 0.241, 0.4329, 0.241, 0.04255]
            \endcode
            
            These values are optimal in the sense that the 5x5 filter obtained by combining 
            this filter with the optimal 5-tap first derivative is the best possible 5x5 approximation to a 
            Gaussian first derivative filter. The equivalent Gaussian has sigma = 0.906.
 
            Postconditions:
            \code
            1. left()  == -2
            2. right() ==  2
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalFirstDerivativeSmoothing5()
    {
        this->initExplicitly(-2, 2) = 0.04255, 0.241, 0.4329, 0.241, 0.04255;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init an optimal 5-tap smoothing filter to be used in the context of second derivative computation.
           This filter must be used in conjunction with the optimal 5-tap second derivative filter 
           (see initOptimalSecondDerivative5()), such that the derivative filter is applied along one dimension, 
           and the smoothing filter along the other. The filter values are 
            
            \code
            [0.0243, 0.23556, 0.48028, 0.23556, 0.0243]
            \endcode
            
            These values are optimal in the sense that the 5x5 filter obtained by combining 
            this filter with the optimal 5-tap second derivative is the best possible 5x5 approximation to a 
            Gaussian second derivative filter. The equivalent Gaussian has sigma = 0.817.
 
            Postconditions:
            \code
            1. left()  == -2
            2. right() ==  2
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalSecondDerivativeSmoothing5()
    {
        this->initExplicitly(-2, 2) = 0.0243, 0.23556, 0.48028, 0.23556, 0.0243;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init a 5-tap filter as defined by Peter Burt in the context of pyramid creation.
            The filter values are
            
            \code
            [a, 0.25, 0.5-2*a, 0.25, a]
            \endcode
            
            The default <tt>a = 0.04785</tt> is optimal in the sense that it minimizes the difference
            to a true Gaussian filter (which would have sigma = 0.975). For other values of <tt>a</tt>, the scale 
            of the most similar Gaussian can be approximated by
            
            \code
            sigma = 5.1 * a + 0.731
            \endcode
 
            Preconditions:
            \code
            0 <= a <= 0.125
            \endcode

            Postconditions:
            \code
            1. left()  == -2
            2. right() ==  2
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initBurtFilter(double a = 0.04785)
    {
        vigra_precondition(a >= 0.0 && a <= 0.125,
            "Kernel1D::initBurtFilter(): 0 <= a <= 0.125 required.");
        this->initExplicitly(-2, 2) = a, 0.25, 0.5 - 2.0*a, 0.25, a;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /**
            Init as a Binomial filter. 'norm' denotes the sum of all bins
            of the kernel.

            Precondition:
            \code
            radius   >= 0
            \endcode

            Postconditions:
            \code
            1. left()  == -radius
            2. right() ==  radius
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == norm
            \endcode
        */
    void initBinomial(int radius, value_type norm);

        /** Init as a Binomial filter with norm 1.
         */
    void initBinomial(int radius)
    {
        initBinomial(radius, one());
    }

        /**
            Init as an Averaging filter. 'norm' denotes the sum of all bins
            of the kernel. The window size is (2*radius+1) * (2*radius+1)

            Precondition:
            \code
            radius   >= 0
            \endcode

            Postconditions:
            \code
            1. left()  == -radius
            2. right() ==  radius
            3. borderTreatment() == BORDER_TREATMENT_CLIP
            4. norm() == norm
            \endcode
        */
    void initAveraging(int radius, value_type norm);

        /** Init as an Averaging filter with norm 1.
         */
    void initAveraging(int radius)
    {
        initAveraging(radius, one());
    }

        /**
           Init as a symmetric gradient filter of the form
           <TT>[ 0.5 * norm, 0.0 * norm, -0.5 * norm]</TT>
           
           <b>Deprecated</b>. Use initSymmetricDifference() instead.

            Postconditions:
            \code
            1. left()  == -1
            2. right() ==  1
            3. borderTreatment() == BORDER_TREATMENT_REPEAT
            4. norm() == norm
            \endcode
        */
    void
    initSymmetricGradient(value_type norm )
    {
        initSymmetricDifference(norm);
        setBorderTreatment(BORDER_TREATMENT_REPEAT);
    }

        /** Init as a symmetric gradient filter with norm 1.
           
           <b>Deprecated</b>. Use initSymmetricDifference() instead.
         */
    void initSymmetricGradient()
    {
        initSymmetricGradient(one());
    }

    void
    initSymmetricDifference(value_type norm );

        /** Init as the 3-tap symmetric difference filter
             The filter values are
             
             \code
             [0.5, 0, -0.5]
             \endcode

            Postconditions:
            \code
            1. left()  == -1
            2. right() ==  1
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
          */
    void initSymmetricDifference()
    {
        initSymmetricDifference(one());
    }

        /**
            Init the 3-tap second difference filter.
            The filter values are
             
             \code
             [1, -2, 1]
             \endcode

            Postconditions:
            \code
            1. left()  == -1
            2. right() ==  1
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1
            \endcode
        */
    void
    initSecondDifference3()
    {
        this->initExplicitly(-1, 1) = 1.0, -2.0, 1.0;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }
    
        /**
            Init an optimal 5-tap first derivative filter.
           This filter must be used in conjunction with the corresponding 5-tap smoothing filter 
           (see initOptimalFirstDerivativeSmoothing5()), such that the derivative filter is applied along one dimension,
            and the smoothing filter along the other.
            The filter values are 
            
            \code
            [0.1, 0.3, 0.0, -0.3, -0.1]
            \endcode
            
            These values are optimal in the sense that the 5x5 filter obtained by combining 
            this filter with the corresponding 5-tap smoothing filter is the best possible 5x5 approximation to a 
            Gaussian first derivative filter. The equivalent Gaussian has sigma = 0.906.
            
            If the filter is instead separably combined with itself, an almost optimal approximation of the
            mixed second Gaussian derivative at scale sigma = 0.899 results.
 
            Postconditions:
            \code
            1. left()  == -2
            2. right() ==  2
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalFirstDerivative5()
    {
        this->initExplicitly(-2, 2) = 0.1, 0.3, 0.0, -0.3, -0.1;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }
    
        /**
            Init an optimal 5-tap second derivative filter.
           This filter must be used in conjunction with the corresponding 5-tap smoothing filter 
           (see initOptimalSecondDerivativeSmoothing5()), such that the derivative filter is applied along one dimension,
            and the smoothing filter along the other.
            The filter values are 
            
            \code
            [0.22075, 0.117, -0.6755, 0.117, 0.22075]
            \endcode
            
            These values are optimal in the sense that the 5x5 filter obtained by combining 
            this filter with the corresponding 5-tap smoothing filter is the best possible 5x5 approximation to a 
            Gaussian second derivative filter. The equivalent Gaussian has sigma = 0.817.
 
            Postconditions:
            \code
            1. left()  == -2
            2. right() ==  2
            3. borderTreatment() == BORDER_TREATMENT_REFLECT
            4. norm() == 1.0
            \endcode
        */
    void initOptimalSecondDerivative5()
    {
        this->initExplicitly(-2, 2) = 0.22075, 0.117, -0.6755, 0.117, 0.22075;
        this->setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

        /** Init the kernel by an explicit initializer list.
            The left and right boundaries of the kernel must be passed.
            A comma-separated initializer list is given after the assignment
            operator. This function is used like this:

            \code
            // define horizontal Roberts filter
            vigra::Kernel1D<float> roberts_gradient_x;

            roberts_gradient_x.initExplicitly(0, 1) = 1.0, -1.0;
            \endcode

            The norm is set to the sum of the initialzer values. If the wrong number of
            values is given, a run-time error results. It is, however, possible to give
            just one initializer. This creates an averaging filter with the given constant:

            \code
            vigra::Kernel1D<float> average5x1;

            average5x1.initExplicitly(-2, 2) = 1.0/5.0;
            \endcode

            Here, the norm is set to value*size().

            <b> Preconditions:</b>

            \code

            1. left <= 0
            2. right >= 0
            3. the number of values in the initializer list
               is 1 or equals the size of the kernel.
            \endcode
        */
    Kernel1D & initExplicitly(int left, int right)
    {
        vigra_precondition(left <= 0,
                     "Kernel1D::initExplicitly(): left border must be <= 0.");
        vigra_precondition(right >= 0,
                     "Kernel1D::initExplicitly(): right border must be >= 0.");

        right_ = right;
        left_ = left;

        kernel_.resize(right - left + 1);

        return *this;
    }

        /** Get iterator to center of kernel

            Postconditions:
            \code

            center()[left()] ... center()[right()] are valid kernel positions
            \endcode
        */
    iterator center()
    {
        return kernel_.begin() - left();
    }

    const_iterator center() const
    {
        return kernel_.begin() - left();
    }

        /** Access kernel value at specified location.

            Preconditions:
            \code

            left() <= location <= right()
            \endcode
        */
    reference operator[](int location)
    {
        return kernel_[location - left()];
    }

    const_reference operator[](int location) const
    {
        return kernel_[location - left()];
    }

        /** left border of kernel (inclusive), always <= 0
        */
    int left() const { return left_; }

        /** right border of kernel (inclusive), always >= 0
        */
    int right() const { return right_; }

        /** size of kernel (right() - left() + 1)
        */
    int size() const { return right_ - left_ + 1; }

        /** current border treatment mode
        */
    BorderTreatmentMode borderTreatment() const
    { return border_treatment_; }

        /** Set border treatment mode.
        */
    void setBorderTreatment( BorderTreatmentMode new_mode)
    { border_treatment_ = new_mode; }

        /** norm of kernel
        */
    value_type norm() const { return norm_; }

        /** set a new norm and normalize kernel, use the normalization formula
            for the given </tt>derivativeOrder</tt>.
        */
    void
    normalize(value_type norm, unsigned int derivativeOrder = 0, double offset = 0.0);

        /** normalize kernel to norm 1.
        */
    void
    normalize()
    {
        normalize(one());
    }

        /** get a const accessor
        */
    ConstAccessor accessor() const { return ConstAccessor(); }

        /** get an accessor
        */
    Accessor accessor() { return Accessor(); }

  private:
    InternalVector kernel_;
    int left_, right_;
    BorderTreatmentMode border_treatment_;
    value_type norm_;
};

template <class ARITHTYPE>
void Kernel1D<ARITHTYPE>::normalize(value_type norm,
                          unsigned int derivativeOrder,
                          double offset)
{
    typedef typename NumericTraits<value_type>::RealPromote TmpType;

    // find kernel sum
    Iterator k = kernel_.begin();
    TmpType sum = NumericTraits<TmpType>::zero();

    if(derivativeOrder == 0)
    {
        for(; k < kernel_.end(); ++k)
        {
            sum += *k;
        }
    }
    else
    {
        unsigned int faculty = 1;
        for(unsigned int i = 2; i <= derivativeOrder; ++i)
            faculty *= i;
        for(double x = left() + offset; k < kernel_.end(); ++x, ++k)
        {
            sum += *k * VIGRA_CSTD::pow(-x, int(derivativeOrder)) / faculty;
        }
    }

    vigra_precondition(sum != NumericTraits<value_type>::zero(),
                    "Kernel1D<ARITHTYPE>::normalize(): "
                    "Cannot normalize a kernel with sum = 0");
    // normalize
    sum = norm / sum;
    k = kernel_.begin();
    for(; k != kernel_.end(); ++k)
    {
        *k = *k * sum;
    }

    norm_ = norm;
}

/***********************************************************************/

template <class ARITHTYPE>
void Kernel1D<ARITHTYPE>::initGaussian(double std_dev,
                                       value_type norm)
{
    vigra_precondition(std_dev >= 0.0,
              "Kernel1D::initGaussian(): Standard deviation must be >= 0.");

    if(std_dev > 0.0)
    {
        Gaussian<ARITHTYPE> gauss(std_dev);

        // first calculate required kernel sizes
        int radius = (int)(3.0 * std_dev + 0.5);
        if(radius == 0)
            radius = 1;

        // allocate the kernel
        kernel_.erase(kernel_.begin(), kernel_.end());
        kernel_.reserve(radius*2+1);

        for(ARITHTYPE x = -radius; x <= radius; ++x)
        {
            kernel_.push_back(gauss(x));
        }
        left_ = -radius;
        right_ = radius;
    }
    else
    {
        kernel_.erase(kernel_.begin(), kernel_.end());
        kernel_.push_back(1.0);
        left_ = 0;
        right_ = 0;
    }

    if(norm != 0.0)
        normalize(norm);
    else
        norm_ = 1.0;

    // best border treatment for Gaussians is BORDER_TREATMENT_REFLECT
    border_treatment_ = BORDER_TREATMENT_REFLECT;
}

/***********************************************************************/

template <class ARITHTYPE>
void Kernel1D<ARITHTYPE>::initDiscreteGaussian(double std_dev,
                                       value_type norm)
{
    vigra_precondition(std_dev >= 0.0,
              "Kernel1D::initDiscreteGaussian(): Standard deviation must be >= 0.");

    if(std_dev > 0.0)
    {
        // first calculate required kernel sizes
        int radius = (int)(3.0*std_dev + 0.5);
        if(radius == 0)
            radius = 1;

        double f = 2.0 / std_dev / std_dev;

        // allocate the working array
        int maxIndex = (int)(2.0 * (radius + 5.0 * VIGRA_CSTD::sqrt((double)radius)) + 0.5);
        InternalVector warray(maxIndex+1);
        warray[maxIndex] = 0.0;
        warray[maxIndex-1] = 1.0;

        for(int i = maxIndex-2; i >= radius; --i)
        {
            warray[i] = warray[i+2] + f * (i+1) * warray[i+1];
            if(warray[i] > 1.0e40)
            {
                warray[i+1] /= warray[i];
                warray[i] = 1.0;
            }
        }

        // the following rescaling ensures that the numbers stay in a sensible range
        // during the rest of the iteration, so no other rescaling is needed
        double er = VIGRA_CSTD::exp(-radius*radius / (2.0*std_dev*std_dev));
        warray[radius+1] = er * warray[radius+1] / warray[radius];
        warray[radius] = er;

        for(int i = radius-1; i >= 0; --i)
        {
            warray[i] = warray[i+2] + f * (i+1) * warray[i+1];
            er += warray[i];
        }

        double scale = norm / (2*er - warray[0]);

        initExplicitly(-radius, radius);
        iterator c = center();

        for(int i=0; i<=radius; ++i)
        {
            c[i] = c[-i] = warray[i] * scale;
        }
    }
    else
    {
        kernel_.erase(kernel_.begin(), kernel_.end());
        kernel_.push_back(norm);
        left_ = 0;
        right_ = 0;
    }

    norm_ = norm;

    // best border treatment for Gaussians is BORDER_TREATMENT_REFLECT
    border_treatment_ = BORDER_TREATMENT_REFLECT;
}

/***********************************************************************/

template <class ARITHTYPE>
void
Kernel1D<ARITHTYPE>::initGaussianDerivative(double std_dev,
                    int order,
                    value_type norm)
{
    vigra_precondition(order >= 0,
              "Kernel1D::initGaussianDerivative(): Order must be >= 0.");

    if(order == 0)
    {
        initGaussian(std_dev, norm);
        return;
    }

    vigra_precondition(std_dev > 0.0,
              "Kernel1D::initGaussianDerivative(): "
              "Standard deviation must be > 0.");

    Gaussian<ARITHTYPE> gauss(std_dev, order);

    // first calculate required kernel sizes
    int radius = (int)(3.0 * std_dev + 0.5 * order + 0.5);
    if(radius == 0)
        radius = 1;

    // allocate the kernels
    kernel_.clear();
    kernel_.reserve(radius*2+1);

    // fill the kernel and calculate the DC component
    // introduced by truncation of the Gaussian
    ARITHTYPE dc = 0.0;
    for(ARITHTYPE x = -radius; x <= radius; ++x)
    {
        kernel_.push_back(gauss(x));
        dc += kernel_[kernel_.size()-1];
    }
    dc /= (2.0*radius + 1.0);

    // remove DC, but only if kernel correction is permitted by a non-zero
    // value for norm
    if(norm != 0.0)
    {
        for(unsigned int i=0; i < kernel_.size(); ++i)
        {
            kernel_[i] -= dc;
        }
    }

    left_ = -radius;
    right_ = radius;

    if(norm != 0.0)
        normalize(norm, order);
    else
        norm_ = 1.0;

    // best border treatment for Gaussian derivatives is
    // BORDER_TREATMENT_REFLECT
    border_treatment_ = BORDER_TREATMENT_REFLECT;
}

/***********************************************************************/

template <class ARITHTYPE>
void
Kernel1D<ARITHTYPE>::initBinomial(int radius,
                                  value_type norm)
{
    vigra_precondition(radius > 0,
              "Kernel1D::initBinomial(): Radius must be > 0.");

    // allocate the kernel
    InternalVector kernel(radius*2+1);

    int i,j;
    for(i=0; i<radius*2+1; ++i) kernel[i] = 0;

    // fill kernel
    typename InternalVector::iterator x = kernel.begin() + radius;
    x[radius] = 1.0;

    for(j=radius-1; j>=-radius; --j)
    {
        for(i=j; i<radius; ++i)
        {
            x[i] = (x[i] + x[i+1]) / 2.0;
        }
        x[radius] /= 2.0;
    }

    // normalize
    kernel_.erase(kernel_.begin(), kernel_.end());
    kernel_.reserve(radius*2+1);

    for(i=0; i<=radius*2+1; ++i)
    {
        kernel_.push_back(kernel[i] * norm);
    }

    left_ = -radius;
    right_ = radius;
    norm_ = norm;

    // best border treatment for Binomial is BORDER_TREATMENT_REFLECT
    border_treatment_ = BORDER_TREATMENT_REFLECT;
}

/***********************************************************************/

template <class ARITHTYPE>
void Kernel1D<ARITHTYPE>::initAveraging(int radius,
                                        value_type norm)
{
    vigra_precondition(radius > 0,
              "Kernel1D::initAveraging(): Radius must be > 0.");

    // calculate scaling
    double scale = 1.0 / (radius * 2 + 1);

    // normalize
    kernel_.erase(kernel_.begin(), kernel_.end());
    kernel_.reserve(radius*2+1);

    for(int i=0; i<=radius*2+1; ++i)
    {
        kernel_.push_back(scale * norm);
    }

    left_ = -radius;
    right_ = radius;
    norm_ = norm;

    // best border treatment for Averaging is BORDER_TREATMENT_CLIP
    border_treatment_ = BORDER_TREATMENT_CLIP;
}

/***********************************************************************/

template <class ARITHTYPE>
void
Kernel1D<ARITHTYPE>::initSymmetricDifference(value_type norm)
{
    kernel_.erase(kernel_.begin(), kernel_.end());
    kernel_.reserve(3);

    kernel_.push_back(0.5 * norm);
    kernel_.push_back(0.0 * norm);
    kernel_.push_back(-0.5 * norm);

    left_ = -1;
    right_ = 1;
    norm_ = norm;

    // best border treatment for symmetric difference is
    // BORDER_TREATMENT_REFLECT
    border_treatment_ = BORDER_TREATMENT_REFLECT;
}

/**************************************************************/
/*                                                            */
/*         Argument object factories for Kernel1D             */
/*                                                            */
/*     (documentation: see vigra/convolution.hxx)             */
/*                                                            */
/**************************************************************/

template <class KernelIterator, class KernelAccessor>
inline
tuple5<KernelIterator, KernelAccessor, int, int, BorderTreatmentMode>
kernel1d(KernelIterator ik, KernelAccessor ka,
       int kleft, int kright, BorderTreatmentMode border)
{
    return
      tuple5<KernelIterator, KernelAccessor, int, int, BorderTreatmentMode>(
                                                ik, ka, kleft, kright, border);
}

template <class T>
inline
tuple5<typename Kernel1D<T>::const_iterator, typename Kernel1D<T>::ConstAccessor,
       int, int, BorderTreatmentMode>
kernel1d(Kernel1D<T> const & k)

{
    return
        tuple5<typename Kernel1D<T>::const_iterator, typename Kernel1D<T>::ConstAccessor,
               int, int, BorderTreatmentMode>(
                                     k.center(),
                                     k.accessor(),
                                     k.left(), k.right(),
                                     k.borderTreatment());
}

template <class T>
inline
tuple5<typename Kernel1D<T>::const_iterator, typename Kernel1D<T>::ConstAccessor,
       int, int, BorderTreatmentMode>
kernel1d(Kernel1D<T> const & k, BorderTreatmentMode border)

{
    return
        tuple5<typename Kernel1D<T>::const_iterator, typename Kernel1D<T>::ConstAccessor,
               int, int, BorderTreatmentMode>(
                                     k.center(),
                                     k.accessor(),
                                     k.left(), k.right(),
                                     border);
}


} // namespace vigra

#endif // VIGRA_SEPARABLECONVOLUTION_HXX
