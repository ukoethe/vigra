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
 
 
#ifndef VIGRA_RECURSIVECONVOLUTION_HXX
#define VIGRA_RECURSIVECONVOLUTION_HXX

#include <cmath>
#include <vector>
#include "utilities.hxx"
#include "numerictraits.hxx"
#include "imageiteratoradapter.hxx"
#include "bordertreatment.hxx"
#include "array_vector.hxx"

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

    The function performs a causal and an anti-causal first or second order 
    recursive filtering with the given filter parameter <TT>b1</TT> and 
    border treatment <TT>border</TT> (first order filter, <TT>b2 = 0</TT>) or parameters 
    <TT>b1, b2</TT> and <TT>BORDER_TREATMENT_REFLECT</TT> (second order filter). Thus, 
    the result is always a filtering with linear phase.
    \f[
        \begin{array}{rcl}
        a_{i, causal} & = & source_i + b1 * a_{i-1, causal} + b2 * a_{i-2, causal} \\
        a_{i, anticausal} & = & source_i + b1 * a_{i+1, anticausal} + b2 * a_{i+2, anticausal} \\
        dest_i & = & \frac{1 - b1 - b2}{1 + b1 + b2}(a_{i, causal} + a_{i, anticausal} - source_i)
        \end{array}
    \f]
   
    The signal's value_type (SrcAccessor::value_type) must be a
    linear space over <TT>double</TT>,
    i.e. addition of source values, multiplication with <TT>double</TT>,
    and <TT>NumericTraits</TT> must be defined.     
    
    <b> Declaration:</b>
    
    <b>First order recursive filter:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
        void recursiveFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                     DestIterator id, DestAccessor ad, 
                     double b1, BorderTreatmentMode border)
    }
    \endcode
    
    <b>Second order recursive filter:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
        void recursiveFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                     DestIterator id, DestAccessor ad, 
                     double b1, double b2)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
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
    -1 < b  < 1
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveFilterLine)

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
                 
    // trivial case: b == 0.0 is an identity filter => simply copy the data and return
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
    typedef typename DestTraits::RealPromote RealPromote;
    
    // store result of causal filtering
    std::vector<TempType> vline(w);
    typename std::vector<TempType>::iterator line = vline.begin();
    
    double norm = (1.0 - b) / (1.0 + b);

    TempType old;
    
    if(border == BORDER_TREATMENT_REPEAT ||
       border == BORDER_TREATMENT_AVOID)
    {
         old = TempType((1.0 / (1.0 - b)) * as(is));
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        is += kernelw;
        old = TempType((1.0 / (1.0 - b)) * as(is));
        for(x = 0; x < kernelw; ++x, --is)
            old = TempType(as(is) + b * old);
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        is = isend - kernelw; 
        old = TempType((1.0 / (1.0 - b)) * as(is));
        for(x = 0; x < kernelw; ++x, ++is)
            old = TempType(as(is) + b * old);
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
        old = TempType(as(is) + b * old);
        line[x] = old;
    }

    // right side of the filter
    if(border == BORDER_TREATMENT_REPEAT ||
       border == BORDER_TREATMENT_AVOID)
    {
        is = isend - 1;
        old = TempType((1.0 / (1.0 - b)) * as(is));
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        old = line[w-2];
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
      is = istart + kernelw - 1;
      old = TempType((1.0 / (1.0 - b)) * as(is));
      for(x = 0; x < kernelw; ++x, --is)
          old = TempType(as(is) + b * old);
    }
    else if(border == BORDER_TREATMENT_CLIP)
    {
        old = NumericTraits<TempType>::zero();
    }
    
    is = isend - 1;
    id += w - 1;
    if(border == BORDER_TREATMENT_CLIP)
    {    
       // correction factors for b
        double bright = b;
        double bleft = VIGRA_CSTD::pow(b, w);

        for(x=w-1; x>=0; --x, --is, --id)
        {    
            TempType f = TempType(b * old);
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
            TempType f = TempType(b * old);
            old = as(is) + f;
            if(x < w - kernelw)
                ad.set(DestTraits::fromRealPromote(RealPromote(norm * (line[x] + f))), id);
        }
    }
    else
    {
        for(x=w-1; x>=0; --x, --is, --id)
        {    
            TempType f = TempType(b * old);
            old = as(is) + f;
            ad.set(DestTraits::fromRealPromote(RealPromote(norm * (line[x] + f))), id);
        }
    }
}
            
/********************************************************/
/*                                                      */
/*            recursiveFilterLine (2nd order)           */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void recursiveFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                         DestIterator id, DestAccessor ad, double b1, double b2)
{
    int w = isend - is;
    SrcIterator istart = is;
    
    int x;
    
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TempType;
    typedef NumericTraits<typename DestAccessor::value_type> DestTraits;
    
    // speichert den Ergebnis der linkseitigen Filterung.
    std::vector<TempType> vline(w+1);
    typename std::vector<TempType>::iterator line = vline.begin();
    
    double norm  = 1.0 - b1 - b2;
    double norm1 = (1.0 - b1 - b2) / (1.0 + b1 + b2);
    double norm2 = norm * norm;
    

    // init left side of filter
    int kernelw = std::min(w-1, std::max(8, (int)(1.0 / norm + 0.5)));  
    is += (kernelw - 2);
    line[kernelw] = as(is);
    line[kernelw-1] = as(is);
    for(x = kernelw - 2; x > 0; --x, --is)
    {
        line[x] = detail::RequiresExplicitCast<TempType>::cast(as(is) + b1 * line[x+1] + b2 * line[x+2]);
    }
    line[0] = detail::RequiresExplicitCast<TempType>::cast(as(is) + b1 * line[1] + b2 * line[2]);
    ++is;
    line[1] = detail::RequiresExplicitCast<TempType>::cast(as(is) + b1 * line[0] + b2 * line[1]);
    ++is;
    for(x=2; x < w; ++x, ++is)
    {
        line[x] = detail::RequiresExplicitCast<TempType>::cast(as(is) + b1 * line[x-1] + b2 * line[x-2]);
    }
    line[w] = line[w-1];

    line[w-1] = detail::RequiresExplicitCast<TempType>::cast(norm1 * (line[w-1] + b1 * line[w-2] + b2 * line[w-3]));
    line[w-2] = detail::RequiresExplicitCast<TempType>::cast(norm1 * (line[w-2] + b1 * line[w] + b2 * line[w-2]));
    id += w-1;
    ad.set(line[w-1], id);
    --id;
    ad.set(line[w-2], id);
    --id;
    for(x=w-3; x>=0; --x, --id, --is)
    {    
        line[x] = detail::RequiresExplicitCast<TempType>::cast(norm2 * line[x] + b1 * line[x+1] + b2 * line[x+2]);
        ad.set(line[x], id);
    }
}
            
/********************************************************/
/*                                                      */
/*              recursiveGaussianFilterLine             */
/*                                                      */
/********************************************************/

// AUTHOR: Sebastian Boppel

/** \brief Compute a 1-dimensional recursive approximation of Gaussian smoothing.

    The function applies a causal and an anti-causal third order recursive filter 
    which optimally approximates the Gaussian filter, as proposed in
    
    I. Young, L. van Vliet: <i>Recursive implementation of the Gaussian filter</i><br>
    Signal Processing 44:139-151, 1995
    
    The formulas for transforming the given scale parameter <tt>sigma</tt> into the actual filter coefficients
    are taken from Luigi Rosa's Matlab implementation.
   
    The signal's value_type (SrcAccessor::value_type) must be a
    linear space over <TT>double</TT>, i.e. addition of source values, multiplication with <TT>double</TT>,
    and <TT>NumericTraits</TT> must be defined.     
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        recursiveGaussianFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                                    DestIterator id, DestAccessor ad, 
                                    double sigma);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    
    \code
    vector<float> src, dest;    
    ...
    
    vigra::DefaultAccessor<vector<float>::iterator, float> FAccessor;
    double sigma = 2.5;
    
    vigra::recursiveGaussianFilterLine(src.begin(), src.end(), FAccessor(), 
                                       dest.begin(), FAccessor(), 
                                       sigma);
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
    0 <= sigma (absolute values are used for negative sigma)
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveGaussianFilterLine)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
recursiveGaussianFilterLine(SrcIterator is, SrcIterator isend, SrcAccessor as,
                            DestIterator id, DestAccessor ad, 
                            double sigma)
{
    //coefficients taken out Luigi Rosa's implementation for Matlab
    double q = 1.31564 * (std::sqrt(1.0 + 0.490811 * sigma*sigma) - 1.0);
    double qq = q*q;
    double qqq = qq*q;
    double b0 = 1.0/(1.57825 + 2.44413*q + 1.4281*qq + 0.422205*qqq);
    double b1 = (2.44413*q + 2.85619*qq + 1.26661*qqq)*b0;
    double b2 = (-1.4281*qq - 1.26661*qqq)*b0;
    double b3 = 0.422205*qqq*b0;
    double B = 1.0 - (b1 + b2 + b3);
    
    int w = isend - is;
    vigra_precondition(w >= 4,
        "recursiveGaussianFilterLine(): line must have at least length 4.");
        
    int kernelw = std::min(w-4, (int)(4.0*sigma));
 
    int x;
    
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TempType;
    typedef NumericTraits<typename DestAccessor::value_type> DestTraits;
    
    // speichert das Ergebnis der linkseitigen Filterung.
    std::vector<TempType> yforward(w);
    
    std::vector<TempType> ybackward(w, 0.0);
    
    // initialise the filter for reflective boundary conditions
    for(x=kernelw; x>=0; --x)
    {
        ybackward[x] = detail::RequiresExplicitCast<TempType>::cast(B*as(is, x) + (b1*ybackward[x+1]+b2*ybackward[x+2]+b3*ybackward[x+3]));
    }

	//from left to right - causal - forward
    yforward[0] = detail::RequiresExplicitCast<TempType>::cast(B*as(is) + (b1*ybackward[1]+b2*ybackward[2]+b3*ybackward[3]));

    ++is;    
    yforward[1] = detail::RequiresExplicitCast<TempType>::cast(B*as(is) + (b1*yforward[0]+b2*ybackward[1]+b3*ybackward[2]));

    ++is;
    yforward[2] = detail::RequiresExplicitCast<TempType>::cast(B*as(is) + (b1*yforward[1]+b2*yforward[0]+b3*ybackward[1]));

    ++is;
    for(x=3; x < w; ++x, ++is)
    {
        yforward[x] = detail::RequiresExplicitCast<TempType>::cast(B*as(is) + (b1*yforward[x-1]+b2*yforward[x-2]+b3*yforward[x-3]));
    }
    
    //from right to left - anticausal - backward
    ybackward[w-1] = detail::RequiresExplicitCast<TempType>::cast(B*yforward[w-1] + (b1*yforward[w-2]+b2*yforward[w-3]+b3*yforward[w-4]));
        
    ybackward[w-2] = detail::RequiresExplicitCast<TempType>::cast(B*yforward[w-2] + (b1*ybackward[w-1]+b2*yforward[w-2]+b3*yforward[w-3]));
    
    ybackward[w-3] = detail::RequiresExplicitCast<TempType>::cast(B*yforward[w-3] + (b1*ybackward[w-2]+b2*ybackward[w-1]+b3*yforward[w-2]));
    
    for(x=w-4; x>=0; --x)
    {
        ybackward[x] = detail::RequiresExplicitCast<TempType>::cast(B*yforward[x]+(b1*ybackward[x+1]+b2*ybackward[x+2]+b3*ybackward[x+3]));
    }

	// output
    for(x=0; x < w; ++x, ++id)
    {
        ad.set(ybackward[x], id);
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void recursiveSmoothLine)

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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void recursiveFirstDerivativeLine)

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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
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
doxygen_overloaded_function(template <...> void recursiveSecondDerivativeLine)

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
    TempType old = detail::RequiresExplicitCast<TempType>::cast((1.0 / (1.0 - b)) * as(is));

    // left side of filter
    for(x=0; x<w; ++x, ++is)
    {
        line[x] = old;
        old = detail::RequiresExplicitCast<TempType>::cast(as(is) + b * old);
    }
    
    // right side of the filter
    --is;
    old = detail::RequiresExplicitCast<TempType>::cast((1.0 / (1.0 - b)) * as(is));
    id += w;
    ++is;
    
    for(x=w-1; x>=0; --x)
    {    
        --is;
        --id;

        TempType f = detail::RequiresExplicitCast<TempType>::cast(old + a * as(is));
        old = detail::RequiresExplicitCast<TempType>::cast(as(is) + b * old);
        ad.set(DestTraits::fromRealPromote(detail::RequiresExplicitCast<TempType>::cast(norm * (line[x] + f))), id);
    }
}
            
/********************************************************/
/*                                                      */
/*                   recursiveFilterX                   */
/*                                                      */
/********************************************************/

/** \brief Performs 1 dimensional recursive filtering (1st and 2nd order) in x direction.

    It calls \ref recursiveFilterLine() for every row of the
    image. See \ref recursiveFilterLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        // first order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterX(SrcImageIterator supperleft, 
                               SrcImageIterator slowerright, SrcAccessor as,
                               DestImageIterator dupperleft, DestAccessor ad, 
                               double b, BorderTreatmentMode border);

        // second order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterX(SrcImageIterator supperleft, 
                               SrcImageIterator slowerright, SrcAccessor as,
                               DestImageIterator dupperleft, DestAccessor ad, 
                               double b1, double b2);
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // first order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterX(
                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double b, BorderTreatmentMode border);

        // second order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterX(
                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double b1, double b2);
            }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothX(srcImageRange(src), destImage(dest), 
           0.5, BORDER_TREATMENT_REFLECT);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveFilterX)

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
/*            recursiveFilterX (2nd order)              */
/*                                                      */
/********************************************************/

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveFilterX(SrcImageIterator supperleft, 
                       SrcImageIterator slowerright, SrcAccessor as,
                       DestImageIterator dupperleft, DestAccessor ad, 
                       double b1, double b2)
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
                             b1, b2);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveFilterX(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
                       double b1, double b2)
{
    recursiveFilterX(src.first, src.second, src.third,
                      dest.first, dest.second, b1, b2);
}
            


/********************************************************/
/*                                                      */
/*               recursiveGaussianFilterX               */
/*                                                      */
/********************************************************/

// AUTHOR: Sebastian Boppel

/** \brief Compute 1 dimensional recursive approximation of Gaussian smoothing in y direction.

    It calls \ref recursiveGaussianFilterLine() for every column of the
    image. See \ref recursiveGaussianFilterLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void 
        recursiveGaussianFilterX(SrcImageIterator supperleft, SrcImageIterator slowerright, SrcAccessor as,
                                 DestImageIterator dupperleft, DestAccessor ad, 
                                 double sigma);
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void 
        recursiveGaussianFilterX(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                 pair<DestImageIterator, DestAccessor> dest, 
                                 double sigma);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveGaussianFilterX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveGaussianFilterX)

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void 
recursiveGaussianFilterX(SrcImageIterator supperleft, SrcImageIterator slowerright, SrcAccessor as,
                         DestImageIterator dupperleft, DestAccessor ad, 
                         double sigma)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int y;
    
    for(y=0; y<h; ++y, ++supperleft.y, ++dupperleft.y)
    {
        typename SrcImageIterator::row_iterator rs = supperleft.rowIterator();
        typename DestImageIterator::row_iterator rd = dupperleft.rowIterator();

        recursiveGaussianFilterLine(rs, rs+w, as, 
                                    rd, ad, 
                                    sigma);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void 
recursiveGaussianFilterX(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                         pair<DestImageIterator, DestAccessor> dest, 
                         double sigma)
{
    recursiveGaussianFilterX(src.first, src.second, src.third,
                             dest.first, dest.second, sigma);
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
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveSmoothX)

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

/** \brief Performs 1 dimensional recursive filtering (1st and 2nd order) in y direction.

    It calls \ref recursiveFilterLine() for every column of the
    image. See \ref recursiveFilterLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        // first order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterY(SrcImageIterator supperleft, 
                               SrcImageIterator slowerright, SrcAccessor as,
                               DestImageIterator dupperleft, DestAccessor ad, 
                               double b, BorderTreatmentMode border);

        // second order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterY(SrcImageIterator supperleft, 
                               SrcImageIterator slowerright, SrcAccessor as,
                               DestImageIterator dupperleft, DestAccessor ad, 
                               double b1, double b2);
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // first order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterY(
                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double b, BorderTreatmentMode border);

        // second order filter
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void recursiveFilterY(
                    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double b1, double b2);
            }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveFilterY(srcImageRange(src), destImage(dest), -0.6, -0.06);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveFilterY)

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
/*            recursiveFilterY (2nd order)              */
/*                                                      */
/********************************************************/

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void recursiveFilterY(SrcImageIterator supperleft, 
                       SrcImageIterator slowerright, SrcAccessor as,
                       DestImageIterator dupperleft, DestAccessor ad, 
                       double b1, double b2)
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
                            b1, b2);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void recursiveFilterY(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest, 
                       double b1, double b2)
{
    recursiveFilterY(src.first, src.second, src.third,
                      dest.first, dest.second, b1, b2);
}
            

/********************************************************/
/*                                                      */
/*               recursiveGaussianFilterY               */
/*                                                      */
/********************************************************/

// AUTHOR: Sebastian Boppel

/** \brief Compute 1 dimensional recursive approximation of Gaussian smoothing in y direction.

    It calls \ref recursiveGaussianFilterLine() for every column of the
    image. See \ref recursiveGaussianFilterLine() for more information about 
    required interfaces and vigra_preconditions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void 
        recursiveGaussianFilterY(SrcImageIterator supperleft, SrcImageIterator slowerright, SrcAccessor as,
                                 DestImageIterator dupperleft, DestAccessor ad, 
                                 double sigma);
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void 
        recursiveGaussianFilterY(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                 pair<DestImageIterator, DestAccessor> dest, 
                                 double sigma);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveGaussianFilterY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveGaussianFilterY)

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void 
recursiveGaussianFilterY(SrcImageIterator supperleft, SrcImageIterator slowerright, SrcAccessor as,
                         DestImageIterator dupperleft, DestAccessor ad, 
                         double sigma)
{
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    int x;
    
    for(x=0; x<w; ++x, ++supperleft.x, ++dupperleft.x)
    {
        typename SrcImageIterator::column_iterator cs = supperleft.columnIterator();
        typename DestImageIterator::column_iterator cd = dupperleft.columnIterator();

        recursiveGaussianFilterLine(cs, cs+h, as, 
                                    cd, ad, 
                                    sigma);
    } 
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void 
recursiveGaussianFilterY(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                         pair<DestImageIterator, DestAccessor> dest, 
                         double sigma)
{
    recursiveGaussianFilterY(src.first, src.second, src.third,
                             dest.first, dest.second, sigma);
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
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSmoothY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveSmoothY)

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
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveFirstDerivativeX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveFirstDerivativeX)

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
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveFirstDerivativeY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveFirstDerivativeY)

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
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSecondDerivativeX(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveSecondDerivativeX)

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
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    
    <b>\#include</b> \<<a href="recursiveconvolution_8hxx-source.html">vigra/recursiveconvolution.hxx</a>\><br>
    Namespace: vigra
    
    \code
    vigra::FImage src(w,h), dest(w,h);    
    ...
    
    vigra::recursiveSecondDerivativeY(srcImageRange(src), destImage(dest), 3.0);
    
    \endcode

*/
doxygen_overloaded_function(template <...> void recursiveSecondDerivativeY)

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
