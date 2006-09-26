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

#ifndef VIGRA_RFFTW_HXX
#define VIGRA_RFFTW_HXX

#include "array_vector.hxx"
#include "fftw.hxx"
#include <rfftw.h>

namespace vigra {

namespace detail {

struct FFTWSinCosConfig
{
    ArrayVector<fftw_real> twiddles;
    ArrayVector<fftw_real> fftwInput;
    ArrayVector<fftw_complex> fftwTmpResult;
    fftw_real norm;
    rfftwnd_plan fftwPlan;
};

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Config>
void 
cosineTransformLineImpl(SrcIterator s, SrcIterator send, SrcAccessor src, 
                        DestIterator d, DestAccessor dest,
                        Config & config)
{
    int size = send - s;
    int ns2 = size / 2;
    int nm1 = size - 1;
    int modn = size % 2;

    if(size <= 0)
        return;
    
    fftw_real const * twiddles = config.twiddles.begin();
    fftw_real * fftwInput = config.fftwInput.begin();
    fftw_complex * fftwTmpResult = config.fftwTmpResult.begin();
    fftw_real norm = config.norm;
    rfftwnd_plan fftwPlan = config.fftwPlan;

    switch(size)
    {
      case 1:
      {
        dest.set(src(s) / norm, d);
        break;
      }
      case 2:
      {
        dest.set((src(s) + src(s, 1)) / norm, d);
        dest.set((src(s) - src(s, 1)) / norm, d, 1);
        break;
      }
      case 3:
      {
        fftw_real x1p3 = src(s) + src(s, 2);
        fftw_real tx2 =  2.0 * src(s, 1);

        dest.set((x1p3 + tx2) / norm, d, 0);
        dest.set((src(s) - src(s, 2)) / norm, d, 1);
        dest.set((x1p3 - tx2) / norm, d, 2);
        break;
      }
      default:
      {
        fftw_real c1 = src(s) - src(s, nm1);
        fftwInput[0] = src(s) + src(s, nm1);
        for(int k=1; k<ns2; ++k)
        {
            int kc = nm1 - k;
            fftw_real t1 = src(s, k) + src(s, kc);
            fftw_real t2 = src(s, k) - src(s, kc);
            c1 = c1 + twiddles[kc] * t2;
            t2 = twiddles[k] * t2;
            fftwInput[k] = t1 - t2;
            fftwInput[kc] = t1 + t2;
        }

        if (modn != 0)
        {
            fftwInput[ns2] = 2.0*src(s, ns2);
        }
        rfftwnd_one_real_to_complex(fftwPlan, fftwInput, fftwTmpResult);
        dest.set(fftwTmpResult[0].re / norm, d, 0);
        for(int k=1; k<(size+1)/2; ++k)
        {
            dest.set(fftwTmpResult[k].re, d, 2*k-1);
            dest.set(fftwTmpResult[k].im, d, 2*k);
        }
        fftw_real xim2 = dest(d, 1);
        dest.set(c1 / norm, d, 1);
        for(int k=3; k<size; k+=2)
        {
            fftw_real xi = dest(d, k);
            dest.set(dest(d, k-2) - dest(d, k-1) / norm, d, k);
            dest.set(xim2 / norm, d, k-1);
            xim2 = xi;
        }
        if (modn != 0)
        {
            dest.set(xim2 / norm, d, size-1);
        }
      }
    }
}

} // namespace detail

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
void cosineTransformX(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                      DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    detail::FFTWSinCosConfig config;

    // horizontal transformation
    int ns2 = w / 2;
    int nm1 = w - 1;
    int modn = w % 2;
    
    config.twiddles.resize(w+1);
    config.fftwInput.resize(w+1);
    config.fftwTmpResult.resize(w+1);
    config.norm = norm;
    config.fftwPlan = rfftw2d_create_plan(1, nm1, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );

    fftw_real dt = M_PI / nm1;
    for(int k=1; k<ns2; ++k)
    {
        fftw_real f = dt * k;
        config.twiddles[k] = 2.0*VIGRA_CSTD::sin(f);
        config.twiddles[nm1-k] = 2.0*VIGRA_CSTD::cos(f);
    }

    for(; sul.y != slr.y; ++sul.y, ++dul.y)
    {
        typename SrcTraverser::row_iterator s = sul.rowIterator();
        typename SrcTraverser::row_iterator send = s + w;
        typename DestTraverser::row_iterator d = dul.rowIterator();
        cosineTransformLineImpl(s, send, src, d, dest, config);
    }

    rfftwnd_destroy_plan(config.fftwPlan);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
void cosineTransformY(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                      DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    detail::FFTWSinCosConfig config;

    // horizontal transformation
    int ns2 = h / 2;
    int nm1 = h - 1;
    int modn = h % 2;
    
    config.twiddles.resize(h + 1);
    config.fftwInput.resize(h + 1);
    config.fftwTmpResult.resize(h + 1);
    config.norm = norm;
    config.fftwPlan = rfftw2d_create_plan(1, nm1, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );

    fftw_real dt = M_PI / nm1;
    for(int k=1; k<ns2; ++k)
    {
        fftw_real f = dt * k;
        config.twiddles[k] = 2.0*VIGRA_CSTD::sin(f);
        config.twiddles[nm1-k] = 2.0*VIGRA_CSTD::cos(f);
    }

    for(; sul.x != slr.x; ++sul.x, ++dul.x)
    {
        typename SrcTraverser::column_iterator s = sul.columnIterator();
        typename SrcTraverser::column_iterator send = s + h;
        typename DestTraverser::column_iterator d = dul.columnIterator();
        cosineTransformLineImpl(s, send, src, d, dest, config);
    }

    rfftwnd_destroy_plan(config.fftwPlan);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void 
realFourierTransformXEvenYEven(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                      DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    BasicImage<fftw_real> tmp(slr - sul);
    cosineTransformX(sul, slr, src, tmp.upperLeft(), tmp.accessor(), 1.0);
    cosineTransformY(tmp.upperLeft(), tmp.lowerRight(), tmp.accessor(), dul, dest, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void 
realFourierTransformXEvenYEven(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                      pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    realFourierTransformXEvenYEven(src.first, src.second, src.third, dest.first, dest.second, norm);
}

} // namespace vigra

#endif // VIGRA_RFFTW_HXX
