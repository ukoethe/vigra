/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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


#ifndef VIGRA_POLARFILTERS_SPATIAL_HXX
#define VIGRA_POLARFILTERS_SPATIAL_HXX

#include <cmath>
#include <vector>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/imageiteratoradapter.hxx"
#include "vigra/bordertreatment.hxx"
#include "vigra/separableconvolution.hxx"

namespace vigra {

namespace detail {

/***********************************************************************/

typedef std::vector<Kernel1D<double> > KernelArray;

void
initGaussianPolarFilters1(double std_dev, KernelArray & k)
{
    typedef Kernel1D<double>::iterator iterator;
    
    vigra_precondition(std_dev >= 0.0,
              "initGaussianPolarFilter1(): "
              "Standard deviation must be >= 0.");
              
    k.resize(4);
                            
    int radius = (int)(4.0*std_dev + 0.5);
    double f = 1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / std_dev;  // norm
    iterator c;
    int ix;

    double a =  2.*0.16116715428 / VIGRA_CSTD::pow(std_dev, 5);
    double b = -2.*0.68931741459 / VIGRA_CSTD::pow(std_dev, 3);
    
    std_dev *= 1.08179074376;
    double sigma22 = -0.5 / std_dev / std_dev;
    
    for(unsigned int i=0; i<k.size(); ++i)
    {
        k[i].initExplicitly(-radius, radius);
        k[i].setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }
    
    c = k[0].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * VIGRA_CSTD::exp(sigma22 * x * x);
    }
    
    c = k[1].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }
    
    c = k[2].center();
    double b2 = b / 3.0;
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * (b2 + a * x * x) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
    
    c = k[3].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * (b + a * x * x) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
}

void
initGaussianPolarFilters2(double std_dev, KernelArray & k)
{
    typedef Kernel1D<double>::iterator iterator;
    
    vigra_precondition(std_dev >= 0.0,
              "initGaussianPolarFilter2(): "
              "Standard deviation must be >= 0.");
              
    k.resize(3);
                            
    int radius = (int)(4.0*std_dev + 0.5);
    double f = 1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / std_dev;  // norm
    iterator c;
    int ix;

    for(unsigned int i=0; i<k.size(); ++i)
    {
        k[i].initExplicitly(-radius, radius);
        k[i].setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }
        
    double sigma2 = std_dev*std_dev;   
    double sigma22 = -0.5 / sigma2;

    c = k[0].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * VIGRA_CSTD::exp(sigma22 * x * x);
    }
    
    c = k[1].center();
    double f1 = f / sigma2;
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f1 * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }
    
    c = k[2].center();
    double f2 = f / (sigma2 * sigma2);
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f2 * (x * x - sigma2) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
}

void
initGaussianPolarFilters3(double std_dev, KernelArray & k)
{
    typedef Kernel1D<double>::iterator iterator;
    
    vigra_precondition(std_dev >= 0.0,
              "initGaussianPolarFilter3(): "
              "Standard deviation must be >= 0.");
              
    k.resize(4);
                            
    int radius = (int)(4.0*std_dev + 0.5);
    double f = 1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / std_dev;  // norm
    iterator c;
    int ix;

    for(unsigned int i=0; i<k.size(); ++i)
    {
        k[i].initExplicitly(-radius, radius);
        k[i].setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }
        
    double a =  0.3229320833 / VIGRA_CSTD::pow(std_dev, 5);
    double b = -1.3786348292 / VIGRA_CSTD::pow(std_dev, 3);
    
    std_dev *= 1.15470053838;
    double sigma22 = -0.5 / std_dev / std_dev;
        
    c = k[0].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * VIGRA_CSTD::exp(sigma22 * x * x);
    }
        
    c = k[1].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }
        
    c = k[2].center();
    double b2 = b / 3.0;
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * (b2 + a * x * x) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
        
    c = k[3].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * (b + a * x * x) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven>
void evenPolarFilters(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                  DestIteratorEven dupperleft_even, DestAccessorEven even,
                  double scale)
{
    vigra_precondition(even.size(dupperleft_even) == 3,
                       "polarFilters2(): image for even output must have 3 bands.");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    typedef typename 
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 3> > TmpImage;    
    typedef typename TmpImage::traverser TmpTraverser;
    TmpImage t(w, h);
    TmpTraverser tul(t.upperLeft());
    TmpTraverser tlr(t.lowerRight());
    
    KernelArray k2;
    initGaussianPolarFilters2(scale, k2);
    
    // calculate filter responses for even filters  
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k2[2], k2[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k2[1], k2[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k2[0], k2[2]);

    // create even tensor from filter responses  
    for(; tul.y != tlr.y; ++tul.y, ++dupperleft_even.y)
    {
        typename TmpTraverser::row_iterator tr = tul.rowIterator();
        typename TmpTraverser::row_iterator trend = tr + w;
        typename DestIteratorEven::row_iterator e = dupperleft_even.rowIterator();
        for(; tr != trend; ++tr, ++e)
        {
            even.setComponent(sq((*tr)[0]) + sq((*tr)[1]), e, 0);
            even.setComponent(-(*tr)[1] * ((*tr)[0] + (*tr)[2]), e, 1);
            even.setComponent(sq((*tr)[1]) + sq((*tr)[2]), e, 2);
        }      
    }
}

} // namespace detail

/** \addtogroup CommonConvolutionFilters Common Filters
*/
//@{

/********************************************************/
/*                                                      */
/*                       polarFilters                   */
/*                                                      */
/********************************************************/

/** \brief Calculate polar 2-jet filter results for an image.

    These functions calculates a spatial domain approximation of
    the polar separable filters as described in
    
    U. Köthe: <it>"Integrated Edge and Junction Detection with the Boundary Tensor"</it>, 
     in: ICCV 03, Proc. of 9th Intl. Conf. on Computer Vision, Nice 2003, vol. 1, 
     pp. 424-431, Los Alamitos: IEEE Computer Society, 2003
     
    The output images must have 3 bands. The first one represents the even filter's
    tensor, the second one contains the odd filter's tensor. The sum of these tensors
    gives the boundary tensor.
    
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorEven, class DestAccessorEven,
                  class DestIteratorOdd, class DestAccessorOdd>
        void polarFilters2(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                          DestIteratorEven dupperleft_even, DestAccessorEven even,
                          DestIteratorOdd dupperleft_odd, DestAccessorOdd odd,
                          double scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorEven, class DestAccessorEven,
                  class DestIteratorOdd, class DestAccessorOdd>
        inline void
        polarFilters2(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorEven, DestAccessorEven> desteven,
                          pair<DestIteratorOdd, DestAccessorOdd> destodd,
                          double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
void polarFilters2(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                  DestIteratorEven dupperleft_even, DestAccessorEven even,
                  DestIteratorOdd dupperleft_odd, DestAccessorOdd odd,
                  double scale)
{
    vigra_precondition(even.size(dupperleft_even) == 3,
                       "polarFilters2(): image for even output must have 3 bands.");
    vigra_precondition(odd.size(dupperleft_odd) == 3,
                       "polarFilters2(): image for odd output must have 3 bands.");

    detail::evenPolarFilters(supperleft, slowerright, sa, 
                             dupperleft_even, even, scale);
    
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    typedef typename 
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 4> > TmpImage;    
    typedef typename TmpImage::traverser TmpTraverser;
    TmpImage t(w, h);
    TmpTraverser tul(t.upperLeft());
    TmpTraverser tlr(t.lowerRight());
    
    detail::KernelArray k1;
    detail::initGaussianPolarFilters1(scale, k1);
    
    // calculate filter responses for odd filters  
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[3], k1[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[2], k1[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[1], k1[2]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[0], k1[3]);
                  
    // create odd tensor from filter responses  
    for(; tul.y != tlr.y; ++tul.y, ++dupperleft_odd.y)
    {
        typename TmpTraverser::row_iterator tr = tul.rowIterator();
        typename TmpTraverser::row_iterator trend = tr + w;
        typename DestIteratorOdd::row_iterator o = dupperleft_odd.rowIterator();
        for(; tr != trend; ++tr, ++o)
        {
            TmpType d0 = (*tr)[0] + (*tr)[2];
            TmpType d1 = -(*tr)[1] - (*tr)[3];
            
            odd.setComponent(sq(d0), o, 0);
            odd.setComponent(d0 * d1, o, 1);
            odd.setComponent(sq(d1), o, 2);
        }      
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
inline 
void polarFilters2(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorEven, DestAccessorEven> even,
                  pair<DestIteratorOdd, DestAccessorOdd> odd,
                  double scale)
{
    polarFilters2(src.first, src.second, src.third,
                 even.first, even.second, odd.first, odd.second, scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
void polarFilters1(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                  DestIteratorEven dupperleft_even, DestAccessorEven even,
                  DestIteratorOdd dupperleft_odd, DestAccessorOdd odd,
                  double scale)
{
    vigra_precondition(even.size(dupperleft_even) == 3,
                       "polarFilters2(): image for even output must have 3 bands.");
    vigra_precondition(odd.size(dupperleft_odd) == 3,
                       "polarFilters2(): image for odd output must have 3 bands.");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    typedef typename 
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 4> > TmpImage;    
    typedef typename TmpImage::traverser TmpTraverser;
    TmpImage t(w, h);
    TmpTraverser tul(t.upperLeft());
    TmpTraverser tlr(t.lowerRight());
 
    detail::KernelArray k2;
    detail::initGaussianPolarFilters2(scale, k2);
    
    // calculate filter responses for even filters  
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k2[2], k2[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k2[0], k2[2]);

    // create even tensor from filter responses  
    for(; tul.y != tlr.y; ++tul.y, ++dupperleft_even.y)
    {
        typename TmpTraverser::row_iterator tr = tul.rowIterator();
        typename TmpTraverser::row_iterator trend = tr + w;
        typename DestIteratorEven::row_iterator e = dupperleft_even.rowIterator();
        for(; tr != trend; ++tr, ++e)
        {
            even.setComponent(0.5*sq((*tr)[0]-(*tr)[1]), e, 0);
            even.setComponent(0, e, 1);
            even.setComponent(0.5*sq((*tr)[0]-(*tr)[1]), e, 2);
        }      
    }
   
    detail::KernelArray k1;
    detail::initGaussianPolarFilters1(scale, k1);
    
    // calculate filter responses for odd filters  
    tmpBand.setIndex(0);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[3], k1[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[2], k1[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[1], k1[2]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t, tmpBand), k1[0], k1[3]);
                  
    // create odd tensor from filter responses  
    tul = t.upperLeft();
    for(; tul.y != tlr.y; ++tul.y, ++dupperleft_odd.y)
    {
        typename TmpTraverser::row_iterator tr = tul.rowIterator();
        typename TmpTraverser::row_iterator trend = tr + w;
        typename DestIteratorOdd::row_iterator o = dupperleft_odd.rowIterator();
        for(; tr != trend; ++tr, ++o)
        {
            TmpType d0 = (*tr)[0] + (*tr)[2];
            TmpType d1 = -(*tr)[1] - (*tr)[3];
            
            odd.setComponent(sq(d0), o, 0);
            odd.setComponent(d0 * d1, o, 1);
            odd.setComponent(sq(d1), o, 2);
        }      
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
inline 
void polarFilters1(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorEven, DestAccessorEven> even,
                  pair<DestIteratorOdd, DestAccessorOdd> odd,
                  double scale)
{
    polarFilters1(src.first, src.second, src.third,
                 even.first, even.second, odd.first, odd.second, scale);
}

/** \brief Calculate polar 3-jet filter results for an image.

    These functions operates similarly to polarFilters2(), but they 
    include 3rd order filter responses into the odd tensor. 
    
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorEven, class DestAccessorEven,
                  class DestIteratorOdd, class DestAccessorOdd>
        void polarFilters3(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                          DestIteratorEven dupperleft_even, DestAccessorEven even,
                          DestIteratorOdd dupperleft_odd, DestAccessorOdd odd,
                          double scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorEven, class DestAccessorEven,
                  class DestIteratorOdd, class DestAccessorOdd>
        inline void
        polarFilters3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorEven, DestAccessorEven> desteven,
                          pair<DestIteratorOdd, DestAccessorOdd> destodd,
                          double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
void polarFilters3(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                  DestIteratorEven dupperleft_even, DestAccessorEven even,
                  DestIteratorOdd dupperleft_odd, DestAccessorOdd odd,
                  double scale)
{
    vigra_precondition(even.size(dupperleft_even) == 3,
                       "polarFilters3(): image for even output must have 3 bands.");
    vigra_precondition(odd.size(dupperleft_odd) == 3,
                       "polarFilters3(): image for odd output must have 3 bands.");

    detail::evenPolarFilters(supperleft, slowerright, sa, 
                             dupperleft_even, even, scale);
    
    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;
    
    typedef typename 
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 4> > TmpImage;    
    typedef typename TmpImage::traverser TmpTraverser;
    TmpImage t1(w, h), t2(w, h);
    TmpTraverser tul1(t1.upperLeft());
    TmpTraverser tlr1(t1.lowerRight());
    TmpTraverser tul2(t2.upperLeft());
    
    detail::KernelArray k1, k3;
    detail::initGaussianPolarFilters1(scale, k1);
    detail::initGaussianPolarFilters3(scale, k3);
    
    // calculate filter responses for odd filters
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t1.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k1[3], k1[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k1[1], k1[2]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k3[3], k3[0]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k3[1], k3[2]);
                  
    tmpBand.setIndex(0);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k1[0], k1[3]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k1[2], k1[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k3[0], k3[3]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k3[2], k3[1]);
                  
    // create odd tensor from filter responses  
    for(; tul1.y != tlr1.y; ++tul1.y, ++tul2.y, ++dupperleft_odd.y)
    {
        typename TmpTraverser::row_iterator tr1 = tul1.rowIterator();
        typename TmpTraverser::row_iterator trend1 = tr1 + w;
        typename TmpTraverser::row_iterator tr2 = tul2.rowIterator();
        typename DestIteratorOdd::row_iterator o = dupperleft_odd.rowIterator();
        for(; tr1 != trend1; ++tr1, ++tr2, ++o)
        {
            TmpType d111 = 0.75*((*tr1)[0]+(*tr1)[1]-(*tr1)[3])+0.25*(*tr1)[2];
            TmpType d122 = 0.25*((*tr1)[0]+(*tr1)[1]-(*tr1)[2])+0.75*(*tr1)[3];
            TmpType d112 = 0.25*((*tr2)[0]+(*tr2)[1]-(*tr2)[2])+0.75*(*tr2)[3];
            TmpType d222 = 0.75*((*tr2)[0]+(*tr2)[1]-(*tr2)[3])+0.25*(*tr2)[2];
            TmpType d2 = sq(d112);
            TmpType d3 = sq(d122);
            
            odd.setComponent(sq(d111) + 2.0*d2 + d3, o, 0);
            odd.setComponent(d111*d112 + 2.0*d112*d122 + d122*d222, o, 1);
            odd.setComponent(d2 + 2.0*d3 + sq(d222), o, 2);
        }      
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
inline 
void polarFilters3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorEven, DestAccessorEven> even,
                  pair<DestIteratorOdd, DestAccessorOdd> odd,
                  double scale)
{
    polarFilters3(src.first, src.second, src.third,
                 even.first, even.second, odd.first, odd.second, scale);
}

//@}

} // namespace vigra

#endif // VIGRA_POLARFILTERS_SPATIAL_HXX
