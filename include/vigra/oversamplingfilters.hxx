/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2003 by Ullrich Koethe                  */
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

#ifndef VIGRA_OVERSAMPLINGFILTER_HXX
#define VIGRA_OVERSAMPLINGFILTER_HXX

#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/array_vector.hxx"

namespace vigra {

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class KernelIter>
void 
oversamplingConvolutionLine(SrcIter s, SrcIter send, SrcAcc src,
                            DestIter d, DestAcc dest,
                            KernelIter intKernel, int intRadius,
                            KernelIter halfintKernel, int halfintRadius)
{
    typedef typename 
        NumericTraits<typename SrcAcc::value_type>::RealPromote
        TmpType;
        
    int wo = send - s;
    int wn = 2*(wo-1)+1;
    int wo2 = 2*wo - 2;
    
    vigra_precondition(2*intRadius+1 <= wo && 2*halfintRadius <= wo,
       "oversamplingConvolutionLine(): kernel larger than line");
    
    int i;
    for(i=0; i<wn; ++i, ++d)
    {
        TmpType sum = NumericTraits<TmpType>::zero();
        KernelIter k;
        int lbound, hbound;
        
        if(i%2 == 0)
        {
            // integer convolution
            lbound = i/2 - intRadius;
            hbound = i/2 + intRadius;
            k = intKernel + 2*intRadius;
        }
        else
        {
            // half-integer convolution
            lbound = i/2 - halfintRadius + 1;
            hbound = i/2 + halfintRadius;
            k = halfintKernel + 2*halfintRadius - 1;
        }
            
        if(lbound < 0 || hbound >= wo)
        {    
            for(int m=lbound; m <= hbound; ++m, --k)
            {
                int mm = (m < 0) ?
                            -m :
                            (m >= wo) ?
                                wo2 - m :
                                m;
                sum += *k * src(s, mm);
            }
        }
        else
        {
            SrcIter ss = s + lbound;
            SrcIter ssend = s + hbound;
            
            for(; ss <= ssend; ++ss, --k)
            {
                sum += *k * src(ss);
            }
        }
        
        dest.set(sum, d);
    }
    
}

template <class Vector>
int makeSmoothingKernels(Vector & intKernel, Vector & halfintKernel, double scale)
{
    int intRadius = int(VIGRA_CSTD::ceil(3.0*scale));
    int halfintRadius = intRadius + 1;
    
    double s2 = (scale > 0.0) ? 
                    scale*scale*2.0 :
                    1.0;
    
    intKernel.resize(2*intRadius+1);
    double sum = VIGRA_CSTD::exp(0.0);
    for(int i = 1; i <= intRadius; ++i)
    {
        sum += 2.0*VIGRA_CSTD::exp(-i*i/s2);
    }
    for(int i = -intRadius; i <= intRadius; ++i)
    {
        intKernel[i + intRadius] = VIGRA_CSTD::exp(-i*i/s2) / sum;
    }
    
    halfintKernel.resize(2*halfintRadius);
    sum = 0.0;
    for(int i = 0; i < halfintRadius; ++i)
    {
        sum += 2.0*VIGRA_CSTD::exp(-(i+0.5)*(i+0.5)/s2);
    }
    for(int i = -halfintRadius; i < halfintRadius; ++i)
    {
        halfintKernel[i + halfintRadius] = VIGRA_CSTD::exp(-(i+0.5)*(i+0.5)/s2) / sum;
    }
    
    return intRadius;
}

template <class Vector>
int makeDerivativeKernels(Vector & intKernel, Vector & halfintKernel, double scale)
{
    int intRadius = (scale > 0.0) ?
                         int(VIGRA_CSTD::ceil(3.0*scale)) :
                         1;
    int halfintRadius = intRadius;
    intKernel.resize(2*intRadius+1);
    halfintKernel.resize(2*halfintRadius);
    
    if (scale == 0.0)
    {
        intKernel[0] = 0.5;
        intKernel[1] = 0.0;
        intKernel[2] = -0.5;
        halfintKernel[0] = 1.0;
        halfintKernel[1] = -1.0;
    }
    else
    {
        double s2 = scale*scale*2.0;

        double sum = 0.0;
        for(int i = 1; i <= intRadius; ++i)
        {
            double x = i;
            sum += 2.0*VIGRA_CSTD::exp(-x*x/s2)*x*x;
        }
        for(int i = -intRadius; i <= intRadius; ++i)
        {
            double x = i;
            intKernel[i + intRadius] = -x*VIGRA_CSTD::exp(-x*x/s2) / sum;
        }

        halfintKernel.resize(2*halfintRadius);
        sum = 0.0;
        for(int i = 1; i <= halfintRadius; ++i)
        {
            double x = i - 0.5;
            sum += 2.0*VIGRA_CSTD::exp(-x*x/s2)*x*x;
        }
        for(int i = -halfintRadius; i < halfintRadius; ++i)
        {
            double x = i+0.5;
            halfintKernel[i + halfintRadius] = -x*VIGRA_CSTD::exp(-x*x/s2) / sum;
        }
    }
    return intRadius;
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc>
void 
oversamplingConvolution(SrcIter sul, SrcIter slr, SrcAcc src,
                        DestIter dul, DestAcc dest,
                        double scale, int xorder, int yorder)
{
    typedef typename 
        NumericTraits<typename SrcAcc::value_type>::RealPromote
        TmpType;
    typedef BasicImage<TmpType> TmpImage;
    typedef typename TmpImage::traverser TmpIter;
        
    vigra_precondition(scale >= 0.0,
       "oversamplingConvolution(): scale must be >= 0");

    int wo = slr.x - sul.x;
    int ho = slr.y - sul.y;
    //int wn = 2*(wo-1)+1;
    int hn = 2*(ho-1)+1;

    ArrayVector<double> intKernel, halfintKernel;
    int intRadius = makeSmoothingKernels(intKernel, halfintKernel, scale);
    int halfintRadius = intRadius + 1;
    
    ArrayVector<double> intDerivativeKernel, halfintDerivativeKernel;
    int intDerivativeRadius = makeDerivativeKernels(intDerivativeKernel, halfintDerivativeKernel, scale);
    int halfintDerivativeRadius = intDerivativeRadius;
    
    vigra_precondition(2*intDerivativeRadius+1 <= wo && 2*halfintDerivativeRadius <= wo &&
                       2*intDerivativeRadius+1 <= ho && 2*halfintDerivativeRadius <= ho,
       "oversamplingConvolution(): kernel larger than image");
    
    TmpImage tmp(wo, hn);
    
    SrcIter s = sul;
    TmpIter t = tmp.upperLeft();
    
    for(; s.x < slr.x; ++s.x, ++t.x)
    {
        typename SrcIter::column_iterator sc = s.columnIterator();
        typename TmpIter::column_iterator tc = t.columnIterator();
        if(yorder)
        {
            oversamplingConvolutionLine(sc, sc+ho, src, tc, tmp.accessor(),
                                  intDerivativeKernel.begin(), intDerivativeRadius, 
                                  halfintDerivativeKernel.begin(), halfintDerivativeRadius);
        }
        else
        {
            oversamplingConvolutionLine(sc, sc+ho, src, tc, tmp.accessor(),
                                  intKernel.begin(), intRadius, 
                                  halfintKernel.begin(), halfintRadius);
        }
    }
    
    t = tmp.upperLeft();
    TmpIter tlr = tmp.lowerRight();
    
    for(; t.y < tlr.y; ++t.y, ++dul.y)
    {
        typename TmpIter::row_iterator tr = t.rowIterator();
        typename DestIter::row_iterator dr = dul.rowIterator();
        if(xorder)
        {
            oversamplingConvolutionLine(tr, tr+wo, tmp.accessor(), dr, dest,
                                  intDerivativeKernel.begin(), intDerivativeRadius, 
                                  halfintDerivativeKernel.begin(), halfintDerivativeRadius);
        }
        else
        {
            oversamplingConvolutionLine(tr, tr+wo, tmp.accessor(), dr, dest,
                                  intKernel.begin(), intRadius, 
                                  halfintKernel.begin(), halfintRadius);
        }
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc>
inline void 
oversamplingConvolution(triple<SrcIter, SrcIter, SrcAcc> src,
                        pair<DestIter, DestAcc> dest,
                        double scale, int xorder, int yorder)
{
    oversamplingConvolution(src.first, src.second, src.third,
                            dest.first, dest.second, scale, xorder, yorder);
}

} // namespace vigra 


#endif /* VIGRA_OVERSAMPLINGFILTER_HXX */
