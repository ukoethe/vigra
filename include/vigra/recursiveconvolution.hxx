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
 
 
#ifndef VIGRA_RECURSIVECONVOLUTION_HXX
#define VIGRA_RECURSIVECONVOLUTION_HXX

#include <cmath>
#include <vector>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/imageiteratoradapter.hxx"
#include "vigra/bordertreatment.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*         Recursive convolution functions              */
/*                                                      */
/********************************************************/

/** \addtogroup RecursiveConvolution Recursive convolution functions
    
    First order recursive filters and their specialization for 
    the exponential filter and its derivatives (1D and separable 2D).
    These filters are very fast, and the speed does not depend on the 
    filter size. 
*/
//@{

/********************************************************/
/*                                                      */
/*                   recursiveFilterLine                */
/*                                                      */
/********************************************************/

/** \brief Performs a 1-dimensional recursive convolution of the source signal.

    The function performs a causal and an anti-causal recursive filtering
    with the given filter parameter <TT>b</TT> and border treatment 
    <TT>border</TT>. Thus, the result is always a filtering with linear phase.
    \f[
        \begin{array}{rcl}
        a_{i, causal} & = & source_i + b * a_{i-1, causal} \\
        a_{i, anticausal} & = & source_i + b * a_{i+1, anticausal} \\
        dest_i & = & \frac{1 - b}{1 + b}(a_{i, causal} + a_{i, anticausal})
        \end{array}
    \f]
   
    The signal's value_type (SrcAccessor::value_type) must be a
    linear space over <TT>double</TT>,
    i.e. addition of source values, multiplication with <TT>double</TT>,
    and <TT>NumericTraits</TT> must be defined.     
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
        void recursiveFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                     DestIterator id, DestAccessor ad, double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    
    \code
    vector<float> src, dest;    
    ...
    
    vigra::DefaultAccessor<vector<float>::iterator, float> FAccessor;
    
    
    vigra::recursiveFilterLine(src.begin(), src.end(), FAccessor(), 
                               dest.begin(), FAccessor(), 
                               0.5, BORDER_TREATMENT_REFLECT);
    \endcode

    <b> Required Interface:</b>
    
    \code
    RandomAccessIterator is, isend;
    RandomAccessIterator id;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(is);
    double d;
    
    s = s + s;
    s = d * s;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), id);

    \endcode

    <b> Preconditions:</b>
    
    \code
    -1 < b < 1
    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void recursiveFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                         DestIterator id, DestAccessor ad, double b, BorderTreatmentMode border)
{
    int w = isend - is;
    SrcIterator istart = is;
    
    int x;
    
    vigra_precondition(-1.0 < b && b < 1.0,
                 "recursiveFilterLine(): -1 < factor < 1 required.\n");
                 
    if(b == 0.0)
    {
        for(; is != isend; ++is, ++id)
        {
            ad.set(as(is), id);
        }
        return;
    }

    double eps = 0.00001;
    int kernelw = std::min(w-1, (int)(VIGRA_CSTD::log(eps)/VIGRA_CSTD::log(VIGRA_CSTD::fabs(b))));
    
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TempType;
    typedef NumericTraits<typename DestAccessor::value_type> DestTraits;
    
    // speichert den Ergebnis der linkseitigen Filterung.
    std::vector<TempType> vline(w);
    typename std::vector<TempType>::iterator line = vline.begin();
    
    double norm = (1.0 - b) / (1.0 + b);

    TempType old;
    
    if(border == BORDER_TREATMENT_REPEAT ||
       border == BORDER_TREATMENT_AVOID)
    {
         old = (1.0 / (1.0 - b)) * as(is);
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        is += kernelw;
        old = (1.0 / (1.0 - b)) * as(is);
        for(x = 0; x < kernelw; ++x, --is)
            old = as(is) + b * old;
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        is = isend - (kernelw + 1); 
        old = (1.0 / (1.0 - b)) * as(is);
	for(x = 0; x < kernelw; ++x, ++is)
	    old = as(is) + b * old;
    }
    else if(border == BORDER_TREATMENT_CLIP)
    {
        old = NumericTraits<TempType>::zero();
    }
    else
        vigra_fail("recursiveFilterLine(): Unknown border treatment mode.\n");

    // left side of filter
    for(x=0, is = istart; x < w; ++x, ++is)
    {
        old = as(is) + b * old;
        line[x] = old;
    }

    // right side of the filter
    if(border == BORDER_TREATMENT_REPEAT ||
       border == BORDER_TREATMENT_AVOID)
    {
        is = isend - 1;
        old = (1.0 / (1.0 - b)) * as(is);
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        is = isend - (kernelw + 1);
        old = (1.0 / (1.0 - b)) * as(is); 
        for(x = 0; x < kernelw; ++x, ++is)
	   old = as(is) + b * old;
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
      is = istart + (kernelw);
      old = (1.0 / (1.0 - b)) * as(is);; 
      for(x = 0; x < kernelw; ++x, --is)
          old = as(is) + b * old;
    }
    else if(border == BORDER_TREATMENT_CLIP)
    {
        old = NumericTraits<TempType>::zero();
    }
    
    is = isend - 1;
    id += w - 1;
    if(border == BORDER_TREATMENT_CLIP)
    {    
      //Korrekturfaktoren für b
        double bright = b;
        double bleft = VIGRA_CSTD::pow(b, w);// b^w

        for(x=w-1; x>=0; --x, --is, --id)
        {    
            TempType f = b * old;
            old = as(is) + f;
            double norm = (1.0 - b) / (1.0 + b - bleft - bright);
            bleft /= b;
            bright *= b;
            ad.set(norm * (line[x] + f), id);
        }
    }
    else if(border == BORDER_TREATMENT_AVOID)
    {
        for(x=w-1; x >= kernelw; --x, --is, --id)
        {    
            TempType f = b * old;
            old = as(is) + f;
            if(x < w - kernelw)
                ad.set(DestTraits::fromRealPromote(norm * (line[x] + f)), id);
        }
    }
    else
    {
        for(x=w-1; x>=0; --x, --is, --id)
        {    
            TempType f = b * old;
            old = as(is) + f;
            ad.set(DestTraits::fromRealPromote(norm * (line[x] + f)), id);
        }
    }
}
            
/********************************************************/
/*                                                      */
/*                    recursiveSmoothLine               */
/*                                                      */
/********************************************************/

/** \brief Convolves the image with a 1-dimensional exponential filter.

    This function calls \ref recursiveFilterLine() with <TT>b = exp(-1.0/scale)</TT>
    and <TT>border = BORDER_TREATMENT_REPEAT</TT>. See 
    \ref recursiveFilterLine() for more documentation.
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
        void recursiveSmoothLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                     DestIterator id, DestAccessor ad, double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    
    \code
    vector<float> src, dest;    
    ...
    
    vigra::DefaultAccessor<vector<float>::iterator, float> FAccessor;
    
    
    vigra::recursiveSmoothLine(src.begin(), src.end(), FAccessor(), 
                        dest.begin(), FAccessor(), 3.0);
    \endcode

    <b> Required Interface:</b>
    
    \code
    RandomAccessIterator is, isend;
    RandomAccessIterator id;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(is);
    double d;
    
    s = s + s;
    s = d * s;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), id);

    \endcode

    <b> Preconditions:</b>
    
    \code
    scale > 0
    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void recursiveSmoothLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                         DestIterator id, DestAccessor ad, double scale)
{
    vigra_precondition(scale >= 0,
                 "recursiveSmoothLine(): scale must be >= 0.\n");
                 
    double b = (scale == 0.0) ? 
                    0.0 :
                    VIGRA_CSTD::exp(-1.0/scale);
    
    recursiveFilterLine(is, isend, as, id, ad, b, BORDER_TREATMENT_REPEAT);
}
            
/********************************************************/
/*                                                      */
/*             recursiveFirstDerivativeLine             */
/*                                                      */
/********************************************************/

/** \brief Performs a 1 dimensional recursive convolution of the source signal.

    It uses the first derivative an exponential  <TT>d/dx exp(-abs(x)/scale)</TT> as 
    a kernel. The signal's value_type (SrcAccessor::value_type) must be a
    linear space over <TT>double</TT>,
    i.e. addition and subtraction of source values, multiplication with 
    <TT>double</TT>, and <TT>NumericTraits</TT> must be defined. Border 
    treatment is always <TT>BORDER_TREATMENT_REPEAT</TT>.
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
        void recursiveFirstDerivativeLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                     DestIterator id, DestAccessor ad, double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    
    \code
    vector<float> src, dest;    
    ...
    
    vigra::DefaultAccessor<vector<float>::iterator, float> FAccessor;
    
    
    vigra::recursiveFirstDerivativeLine(src.begin(), src.end(), FAccessor(), 
                        dest.begin(), FAccessor(), 3.0);
    \endcode

    <b> Required Interface:</b>
    
    \code
    RandomAccessIterator is, isend;
    RandomAccessIterator id;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(is);
    double d;
    
    s = s + s;
    s = -s;
    s = d * s;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), id);

    \endcode

    <b> Preconditions:</b>
    
    \code
    scale > 0
    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void recursiveFirstDerivativeLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                         DestIterator id, DestAccessor ad, double scale)
{
    vigra_precondition(scale > 0,
                 "recursiveFirstDerivativeLine(): scale must be > 0.\n");

    int w = isend -is;
    
    int x;
    
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote 
    TempType;
    typedef NumericTraits<typename DestAccessor::value_type> DestTraits;

    std::vector<TempType> vline(w);
    typename std::vector<TempType>::iterator line = vline.begin();
    
    double b = VIGRA_CSTD::exp(-1.0/scale);
    double norm = (1.0 - b) * (1.0 - b) / 2.0 / b;
    TempType old = (1.0 / (1.0 - b)) * as(is);

    // left side of filter
    for(x=0; x<w; ++x, ++is)
    {
        old = as(is) + b * old;
        line[x] = -old;
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

        old = as(is) + b * old;

        ad.set(DestTraits::fromRealPromote(norm * (line[x] + old)), id);
    }
}
            
/********************************************************/
/*                                                      */
/*            recursiveSecondDerivativeLine             */
/*                                                      */
/********************************************************/

/** \brief Performs a 1 dimensional recursive convolution of the source signal.

    It uses the second derivative an exponential  <TT>d2/dx2 exp(-abs(x)/scale)</TT> as 
    a kernel. The signal's value_type (SrcAccessor::value_type) must be a
    linear space over <TT>double</TT>,
    i.e. addition and subtraction of source values, multiplication with 
    <TT>double</TT>, and <TT>NumericTraits</TT> must be defined. Border 
    treatment is always <TT>BORDER_TREATMENT_REPEAT</TT>.
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
        void recursiveSecondDerivativeLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                     DestIterator id, DestAccessor ad, double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    
    \code
    vector<float> src, dest;    
    ...
    
    vigra::DefaultAccessor<vector<float>::iterator, float> FAccessor;
    
    
    vigra::recursiveSecondDerivativeLine(src.begin(), src.end(), FAccessor(), 
                        dest.begin(), FAccessor(), 3.0);
    \endcode

    <b> Required Interface:</b>
    
    \code
    RandomAccessIterator is, isend;
    RandomAccessIterator id;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(is);
    double d;
    
    s = s + s;
    s = s - s;
    s = d * s;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), id);

    \endcode

    <b> Preconditions:</b>
    
    \code
    scale > 0
    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void recursiveSecondDerivativeLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                         DestIterator id, DestAccessor ad, double scale)
{
    vigra_precondition(scale > 0,
                 "recursiveSecondDerivativeLine(): scale must be > 0.\n");

    int w = isend -is;
    
    int x;
    
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote 
    TempType;
    typedef NumericTraits<typename DestAccessor::value_type> DestTraits;
    
    std::vector<TempType> vline(w);
    typename std::vector<TempType>::iterator line = vline.begin();
        
    double b = VIGRA_CSTD::exp(-1.0/scale);
    double a = -2.0 / (1.0 - b);
    double norm = (1.0 - b) * (1.0 - b) * (1.0 - b) / (1.0 + b);
    TempType old = (1.0 / (1.0 - b)) * as(is);

    // left side of filter
    for(x=0; x<w; ++x, ++is)
    {
        line[x] = old;
        old = as(is) + b * old;
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

        TempType f = old + a * as(is);
        old = as(is) + b * old;
        ad.set(DestTraits::fromRealPromote(norm * (line[x] + f)), id);
    }
}
            
/********************************************************/
/*                                                      */
/*                   recursiveFilterX                   */
/*                                                      */
/********************************************************/

/** \brief Performs 1 dimensional recursive smoothing in x direction.

    It calls \ref recursiveFilterLine() for every row of the
    image. See \ref recursiveFilterLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterX(SrcImageIterator supperleft, 
                               SrcImageIterator slowerright, SrcAccessor as,
                               DestImageIterator dupperleft, DestAccessor ad, 
                               double b, BorderTreatmentMode border);
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterX(
                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double b, BorderTreatmentMode border);
            }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothX(srcImageRange(src), destImage(dest), 
           0.5, BORDER_TREATMENT_REFLECT);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveFilterX(SrcImageIterator supperleft, 
                       SrcImageIterator slowerright, SrcAccessor as,
                       DestImageIterator dupperleft, DestAccessor ad, 
                       double b, BorderTreatmentMode border)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int y;
    
    for(y=0; y<h; ++y, ++supperleft.y, ++dupperleft.y)
    {
        typename SrcImageIterator::row_iterator rs = supperleft.rowIterator();
        typename DestImageIterator::row_iterator rd = dupperleft.rowIterator();

        recursiveFilterLine(rs, rs+w, as, 
                             rd, ad, 
                             b, border);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveFilterX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double b, BorderTreatmentMode border)
{
    recursiveFilterX(src.first, src.second, src.third,
                      dest.first, dest.second, b, border);
}
            
/********************************************************/
/*                                                      */
/*                    recursiveSmoothX                  */
/*                                                      */
/********************************************************/

/** \brief Performs 1 dimensional recursive smoothing in x direction.

    It calls \ref recursiveSmoothLine() for every row of the
    image. See \ref recursiveSmoothLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSmoothX(SrcImageIterator supperleft, 
                  SrcImageIterator slowerright, SrcAccessor as,
                  DestImageIterator dupperleft, DestAccessor ad, 
                  double scale)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSmoothX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveSmoothX(SrcImageIterator supperleft, 
                      SrcImageIterator slowerright, SrcAccessor as,
                      DestImageIterator dupperleft, DestAccessor ad, 
              double scale)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int y;
    
    for(y=0; y<h; ++y, ++supperleft.y, ++dupperleft.y)
    {
        typename SrcImageIterator::row_iterator rs = supperleft.rowIterator();
        typename DestImageIterator::row_iterator rd = dupperleft.rowIterator();

        recursiveSmoothLine(rs, rs+w, as, 
                            rd, ad, 
                            scale);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveSmoothX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
        double scale)
{
    recursiveSmoothX(src.first, src.second, src.third,
                     dest. first, dest.second, scale);
}
            
/********************************************************/
/*                                                      */
/*                     recursiveFilterY                 */
/*                                                      */
/********************************************************/

/** \brief Performs 1 dimensional recursive smoothing in y direction.

    It calls \ref recursiveFilterLine() for every column of the
    image. See \ref recursiveFilterLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterY(SrcImageIterator supperleft, 
                              SrcImageIterator slowerright, SrcAccessor as,
                              DestImageIterator dupperleft, DestAccessor ad, 
                              double b, BorderTreatmentMode border);
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterY(
                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double b, BorderTreatmentMode border);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveFilterY(SrcImageIterator supperleft, 
                       SrcImageIterator slowerright, SrcAccessor as,
                       DestImageIterator dupperleft, DestAccessor ad, 
                       double b, BorderTreatmentMode border)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int x;
    
    for(x=0; x<w; ++x, ++supperleft.x, ++dupperleft.x)
    {
        typename SrcImageIterator::column_iterator cs = supperleft.columnIterator();
        typename DestImageIterator::column_iterator cd = dupperleft.columnIterator();

        recursiveFilterLine(cs, cs+h, as, 
                            cd, ad, 
                            b, border);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveFilterY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double b, BorderTreatmentMode border)
{
    recursiveFilterY(src.first, src.second, src.third,
                      dest.first, dest.second, b, border);
}

/********************************************************/
/*                                                      */
/*                     recursiveSmoothY                 */
/*                                                      */
/********************************************************/

/** \brief Performs 1 dimensional recursive smoothing in y direction.

    It calls \ref recursiveSmoothLine() for every column of the
    image. See \ref recursiveSmoothLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSmoothY(SrcImageIterator supperleft, 
                  SrcImageIterator slowerright, SrcAccessor as,
                  DestImageIterator dupperleft, DestAccessor ad, 
                  double scale)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSmoothY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveSmoothY(SrcImageIterator supperleft, 
                      SrcImageIterator slowerright, SrcAccessor as,
                      DestImageIterator dupperleft, DestAccessor ad, 
              double scale)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int x;
    
    for(x=0; x<w; ++x, ++supperleft.x, ++dupperleft.x)
    {
        typename SrcImageIterator::column_iterator cs = supperleft.columnIterator();
        typename DestImageIterator::column_iterator cd = dupperleft.columnIterator();

        recursiveSmoothLine(cs, cs+h, as, 
                            cd, ad, 
                            scale);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveSmoothY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
{
    recursiveSmoothY(src.first, src.second, src.third,
                     dest. first, dest.second, scale);
}
            
/********************************************************/
/*                                                      */
/*              recursiveFirstDerivativeX               */
/*                                                      */
/********************************************************/

/** \brief Recursively calculates the 1 dimensional first derivative in x 
    direction.
    
    It calls \ref recursiveFirstDerivativeLine() for every 
    row of the image. See \ref recursiveFirstDerivativeLine() for more 
    information about required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveFirstDerivativeX(SrcImageIterator supperleft, 
                  SrcImageIterator slowerright, SrcAccessor as,
                  DestImageIterator dupperleft, DestAccessor ad, 
                  double scale)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveFirstDerivativeX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveFirstDerivativeX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveFirstDerivativeX(SrcImageIterator supperleft, 
                      SrcImageIterator slowerright, SrcAccessor as,
                      DestImageIterator dupperleft, DestAccessor ad, 
              double scale)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int y;
    
    for(y=0; y<h; ++y, ++supperleft.y, ++dupperleft.y)
    {
        typename SrcImageIterator::row_iterator rs = supperleft.rowIterator();
        typename DestImageIterator::row_iterator rd = dupperleft.rowIterator();

        recursiveFirstDerivativeLine(rs, rs+w, as, 
                                     rd, ad, 
                                     scale);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveFirstDerivativeX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
        double scale)
{
    recursiveFirstDerivativeX(src.first, src.second, src.third,
                          dest. first, dest.second, scale);
}
            
/********************************************************/
/*                                                      */
/*              recursiveFirstDerivativeY               */
/*                                                      */
/********************************************************/

/** \brief Recursively calculates the 1 dimensional first derivative in y 
    direction.
    
    It calls \ref recursiveFirstDerivativeLine() for every 
    column of the image. See \ref recursiveFirstDerivativeLine() for more 
    information about required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveFirstDerivativeY(SrcImageIterator supperleft, 
                  SrcImageIterator slowerright, SrcAccessor as,
                  DestImageIterator dupperleft, DestAccessor ad, 
                  double scale)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveFirstDerivativeY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveFirstDerivativeY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveFirstDerivativeY(SrcImageIterator supperleft, 
                      SrcImageIterator slowerright, SrcAccessor as,
                      DestImageIterator dupperleft, DestAccessor ad, 
              double scale)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int x;
    
    for(x=0; x<w; ++x, ++supperleft.x, ++dupperleft.x)
    {
        typename SrcImageIterator::column_iterator cs = supperleft.columnIterator();
        typename DestImageIterator::column_iterator cd = dupperleft.columnIterator();

        recursiveFirstDerivativeLine(cs, cs+h, as, 
                                     cd, ad, 
                                     scale);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveFirstDerivativeY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
        double scale)
{
    recursiveFirstDerivativeY(src.first, src.second, src.third,
                          dest. first, dest.second, scale);
}
            
/********************************************************/
/*                                                      */
/*             recursiveSecondDerivativeX               */
/*                                                      */
/********************************************************/

/** \brief Recursively calculates the 1 dimensional second derivative in x 
    direction.
    
    It calls \ref recursiveSecondDerivativeLine() for every 
    row of the image. See \ref recursiveSecondDerivativeLine() for more 
    information about required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSecondDerivativeX(SrcImageIterator supperleft, 
                  SrcImageIterator slowerright, SrcAccessor as,
                  DestImageIterator dupperleft, DestAccessor ad, 
                  double scale)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSecondDerivativeX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSecondDerivativeX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveSecondDerivativeX(SrcImageIterator supperleft, 
                      SrcImageIterator slowerright, SrcAccessor as,
                      DestImageIterator dupperleft, DestAccessor ad, 
              double scale)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int y;
    
    for(y=0; y<h; ++y, ++supperleft.y, ++dupperleft.y)
    {
        typename SrcImageIterator::row_iterator rs = supperleft.rowIterator();
        typename DestImageIterator::row_iterator rd = dupperleft.rowIterator();

        recursiveSecondDerivativeLine(rs, rs+w, as, 
                                      rd, ad, 
                                      scale);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveSecondDerivativeX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
        double scale)
{
    recursiveSecondDerivativeX(src.first, src.second, src.third,
                          dest. first, dest.second, scale);
}
            
/********************************************************/
/*                                                      */
/*             recursiveSecondDerivativeY               */
/*                                                      */
/********************************************************/

/** \brief Recursively calculates the 1 dimensional second derivative in y 
    direction.
    
    It calls \ref recursiveSecondDerivativeLine() for every 
    column of the image. See \ref recursiveSecondDerivativeLine() for more 
    information about required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSecondDerivativeY(SrcImageIterator supperleft, 
                  SrcImageIterator slowerright, SrcAccessor as,
                  DestImageIterator dupperleft, DestAccessor ad, 
                  double scale)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void recursiveSecondDerivativeY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
            double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSecondDerivativeY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveSecondDerivativeY(SrcImageIterator supperleft, 
                      SrcImageIterator slowerright, SrcAccessor as,
                      DestImageIterator dupperleft, DestAccessor ad, 
              double scale)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int x;
    
    for(x=0; x<w; ++x, ++supperleft.x, ++dupperleft.x)
    {
        typename SrcImageIterator::column_iterator cs = supperleft.columnIterator();
        typename DestImageIterator::column_iterator cd = dupperleft.columnIterator();

        recursiveSecondDerivativeLine(cs, cs+h, as, 
                                      cd, ad, 
                                      scale);
    }
}
            
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveSecondDerivativeY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
        double scale)
{
    recursiveSecondDerivativeY(src.first, src.second, src.third,
                          dest. first, dest.second, scale);
}
            
//@}

} // namespace vigra

#endif // VIGRA_RECURSIVECONVOLUTION_HXX
