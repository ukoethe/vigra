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

#ifndef VIGRA_RESAMPLING_CONVOLUTION_HXX
#define VIGRA_RESAMPLING_CONVOLUTION_HXX

#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/rational.hxx"

namespace vigra {

namespace resampling_detail
{

struct MapTargetToSourceCoordinate
{
    MapTargetToSourceCoordinate(Rational<int> const & samplingRatio, 
                                Rational<int> const & offset)
    : samplingRatio_(samplingRatio),
      offset_(offset)
    {}
    
    int operator()(int i) const
    {
        Rational<int> r = i / samplingRatio_ + offset_;
        return r.numerator() / r.denominator();
    }
    
    double toDouble(int i) const
    {
        Rational<int> r = i / samplingRatio_ + offset_;
        return double(r.numerator()) / double(r.denominator());
    }
    
    Rational<int> samplingRatio_, offset_;
};

} // namespace resampling_detail

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class KernelArray,
          class Functor>
void 
resamplingConvolveLine(SrcIter s, SrcIter send, SrcAcc src,
                       DestIter d, DestIter dend, DestAcc dest,
                       KernelArray const & kernels,
                       Functor mapTargetToSourceCoordinate)
{
    typedef typename 
        NumericTraits<typename SrcAcc::value_type>::RealPromote
        TmpType;
    typedef typename KernelArray::value_type Kernel;
    typedef typename Kernel::const_iterator KernelIter;
    
    int wo = send - s;
    int wn = dend - d;
    int wo2 = 2*wo - 2;
    
    int i;
    typename KernelArray::const_iterator kernel = kernels.begin();
    for(i=0; i<wn; ++i, ++d, ++kernel)
    {
        // use the kernels periodically
        if(kernel == kernels.end())
            kernel = kernels.begin();
        
        // calculate current target point into source location
        int is = mapTargetToSourceCoordinate(i);
        
        TmpType sum = NumericTraits<TmpType>::zero();

        int lbound = is - kernel->right(), 
            hbound = is - kernel->left();
                    
        KernelIter k = kernel->center() + kernel->right();
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

template <class Kernel, class MapCoordinate, class KernelArray>
void
createResamplingKernels(Kernel const & kernel, 
             MapCoordinate const & mapCoordinate, KernelArray & kernels)
{
    for(unsigned int idest = 0; idest < kernels.size(); ++idest)
    {
        int isrc = mapCoordinate(idest);
        double idsrc = mapCoordinate.toDouble(idest);
        double offset = isrc - idsrc;
        int left = VIGRA_CSTD::floor(offset - kernel.radius());
        int right = VIGRA_CSTD::ceil(offset + kernel.radius());
        kernels[idest].initExplicitly(left, right);
        
        double x = left - offset;
        for(int i = left; i <= right; ++i, ++x)
            kernels[idest][i] = kernel(x);
        // normalize!???
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
void 
resamplingConvolveX(SrcIter sul, SrcIter slr, SrcAcc src,
                    DestIter dul, DestIter dlr, DestAcc dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    int wold = slr.x - sul.x;
    int hold = slr.y - sul.y;
    int wnew = dlr.x - dul.x;
    
    int period = lcm(samplingRatio.numerator(), samplingRatio.denominator());
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(samplingRatio, offset);
    
    ArrayVector<Kernel1D<double> > kernels(period);
    
    createResamplingKernels(kernel, mapCoordinate, kernels);

    for(; sul.y < slr.y; ++sul.y, ++dul.y)
    {
        typename SrcIter::row_iterator sr = sul.rowIterator();
        typename DestIter::row_iterator dr = dul.rowIterator();
        resamplingConvolveLine(sr, sr+wold, src, dr, dr+wnew, dest,
                               kernels, mapCoordinate);
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
inline void 
resamplingConvolveX(triple<SrcIter, SrcIter, SrcAcc> src,
                    triple<DestIter, DestIter, DestAcc> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveX(src.first, src.second, src.third,
                        dest.first, dest.second, dest.third,
                        kernel, samplingRatio, offset);
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
void 
resamplingConvolveY(SrcIter sul, SrcIter slr, SrcAcc src,
                    DestIter dul, DestIter dlr, DestAcc dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    int wold = slr.x - sul.x;
    int hold = slr.y - sul.y;
    int hnew = dlr.y - dul.y;
    
    int period = lcm(samplingRatio.numerator(), samplingRatio.denominator());
    
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(samplingRatio, offset);
    
    ArrayVector<Kernel1D<double> > kernels(period);
    
    createResamplingKernels(kernel, mapCoordinate, kernels);

    for(; sul.x < slr.x; ++sul.x, ++dul.x)
    {
        typename SrcIter::column_iterator sc = sul.columnIterator();
        typename DestIter::column_iterator dc = dul.columnIterator();
        resamplingConvolveLine(sc, sc+hold, src, dc, dc+hnew, dest,
                               kernels, mapCoordinate);
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
inline void 
resamplingConvolveY(triple<SrcIter, SrcIter, SrcAcc> src,
                    triple<DestIter, DestIter, DestAcc> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveY(src.first, src.second, src.third,
                        dest.first, dest.second, dest.third,
                        kernel, samplingRatio, offset);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Kernel>
void resamplingConvolveImage(SrcIterator sul,SrcIterator slr, SrcAccessor src,
                   DestIterator dul, DestIterator dlr, DestAccessor dest,
                   Kernel const & kx, 
                   Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                   Kernel const & ky, 
                   Rational<int> const & samplingRatioY, Rational<int> const & offsetY)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    
    BasicImage<TmpType> tmp(dlr.x - dul.x, slr.y - sul.y);

    resamplingConvolveX(srcIterRange(sul, slr, src),
                        destImageRange(tmp), 
                        kx, samplingRatioX, offsetX);
    resamplingConvolveY(srcImageRange(tmp),
                        destIterRange(dul, dlr, dest), 
                        ky, samplingRatioY, offsetY);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Kernel>
inline void
resamplingConvolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   triple<DestIterator, DestIterator, DestAccessor> dest,
                   Kernel const & kx, 
                   Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                   Kernel const & ky, 
                   Rational<int> const & samplingRatioY, Rational<int> const & offsetY)
{
    resamplingConvolveImage(src.first, src.second, src.third,
                            dest.first, dest.second, dest.third,
                            kx, samplingRatioX, offsetX,
                            ky, samplingRatioY, offsetY);
}

} // namespace vigra 


#endif /* VIGRA_RESAMPLING_CONVOLUTION_HXX */
