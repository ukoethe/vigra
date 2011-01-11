/************************************************************************/
/*                                                                      */
/*               Copyright 2009-2010 by Ullrich Koethe                  */
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

#ifndef VIGRA_MULTI_FFT_HXX
#define VIGRA_MULTI_FFT_HXX

#include "fftw3.hxx"
#include "multi_array.hxx"
#include "navigator.hxx"
#include "copyimage.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                    Fourier Transform                 */
/*                                                      */
/********************************************************/

/** \addtogroup FourierTransform Fast Fourier Transform

    This documentation describes the VIGRA interface to FFTW version 3. The interface
    to the old FFTW version 2 (file "vigra/fftw.hxx") is deprecated.

    VIGRA uses the <a href="http://www.fftw.org/">FFTW Fast Fourier
    Transform</a> package to perform Fourier transformations. VIGRA
    provides a wrapper for FFTW's complex number type (FFTWComplex),
    but FFTW's functions are used verbatim. If the image is stored as
    a FFTWComplexImage, the simplest call to an FFT function is like this:

    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // create a plan with estimated performance optimization
    fftw_plan forwardPlan = fftw_plan_dft_2d(height, width,
                                (fftw_complex *)spatial.begin(), (fftw_complex *)fourier.begin(),
                                FFTW_FORWARD, FFTW_ESTIMATE );
    // calculate FFT (this can be repeated as often as needed,
    //                with fresh data written into the source array)
    fftw_execute(forwardPlan);

    // release the plan memory
    fftw_destroy_plan(forwardPlan);

    // likewise for the inverse transform
    fftw_plan backwardPlan = fftw_plan_dft_2d(height, width,
                                 (fftw_complex *)fourier.begin(), (fftw_complex *)spatial.begin(),
                                 FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(backwardPlan);
    fftw_destroy_plan(backwardPlan);

    // do not forget to normalize the result according to the image size
    transformImage(srcImageRange(spatial), destImage(spatial),
                   std::bind1st(std::multiplies<FFTWComplex>(), 1.0 / width / height));
    \endcode

    Note that in the creation of a plan, the height must be given
    first. Note also that <TT>spatial.begin()</TT> may only be passed
    to <TT>fftw_plan_dft_2d</TT> if the transform shall be applied to the
    entire image. When you want to restrict operation to an ROI, you
    can create a copy of the ROI in an image of appropriate size, or
    you may use the Guru interface to FFTW.

    More information on using FFTW can be found <a href="http://www.fftw.org/doc/">here</a>.

    FFTW produces fourier images that have the DC component (the
    origin of the Fourier space) in the upper left corner. Often, one
    wants the origin in the center of the image, so that frequencies
    always increase towards the border of the image. This can be
    achieved by calling \ref moveDCToCenter(). The inverse
    transformation is done by \ref moveDCToUpperLeft().

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
    Namespace: vigra
*/

/** \addtogroup FourierTransform
*/
//@{

/********************************************************/
/*                                                      */
/*                     moveDCToCenter                   */
/*                                                      */
/********************************************************/

/** \brief Rearrange the quadrants of a Fourier image so that the origin is
          in the image center.

    FFTW produces fourier images where the DC component (origin of
    fourier space) is located in the upper left corner of the
    image. The quadrants are placed like this (using a 4x4 image for
    example):

    \code
            DC 4 3 3
             4 4 3 3
             1 1 2 2
             1 1 2 2
    \endcode

    After applying the function, the quadrants are at their usual places:

    \code
            2 2  1 1
            2 2  1 1
            3 3 DC 4
            3 3  4 4
    \endcode

    This transformation can be reversed by \ref moveDCToUpperLeft().
    Note that the transformation must not be executed in place - input
    and output images must be different.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
        void moveDCToCenter(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void moveDCToCenter(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest);
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<vigra/fftw3.hxx\><br>
        Namespace: vigra

    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // create a plan with estimated performance optimization
    fftw_plan forwardPlan = fftw_plan_dft_2d(height, width,
                                (fftw_complex *)spatial.begin(), (fftw_complex *)fourier.begin(),
                                FFTW_FORWARD, FFTW_ESTIMATE );
    // calculate FFT
    fftw_execute(forwardPlan);

    vigra::FFTWComplexImage rearrangedFourier(width, height);
    moveDCToCenter(srcImageRange(fourier), destImage(rearrangedFourier));

    // delete the plan
    fftw_destroy_plan(forwardPlan);
    \endcode
*/
template <unsigned int N, class T, class C>
void moveDCToCenter(MultiArrayView<N, T, C> a)
{
    typedef typename MultiArrayView<N, T, C>::traverser Traverser;
    typedef MultiArrayNavigator<Traverser, N> Navigator;
    typedef typename Navigator::iterator Iterator;
    
    for(unsigned int d = 0; d < N; ++d)
    {
        Navigator nav(a.traverser_begin(), a.shape(), d);

        for( ; nav.hasMore(); nav++ )
        {
            Iterator i = nav.begin();
            int s  = nav.end() - i;
            int s2 = s/2;
                
            if(even(s))
            {
                for(int k=0; k<s2; ++k)
                {
                    std::swap(i[k], i[k+s2]);
                }
            }
            else            
            {
                T v = i[0];
                for(int k=0; k<s2; ++k)
                {
                    i[k] = i[k+s2+1];
                    i[k+s2+1] = i[k+1];
                }
                i[s2] = v;
            }
        }
    }
}

template <unsigned int N, class T, class C>
void moveDCToUpperLeft(MultiArrayView<N, T, C> a)
{
    typedef typename MultiArrayView<N, T, C>::traverser Traverser;
    typedef MultiArrayNavigator<Traverser, N> Navigator;
    typedef typename Navigator::iterator Iterator;
    
    for(unsigned int d = 0; d < N; ++d)
    {
        Navigator nav(a.traverser_begin(), a.shape(), d);

        for( ; nav.hasMore(); nav++ )
        {
            Iterator i = nav.begin();
            int s  = nav.end() - i;
            int s2 = s/2;
                
            if(even(s))
            {
                for(int k=0; k<s2; ++k)
                {
                    std::swap(i[k], i[k+s2]);
                }
            }
            else            
            {
                T v = i[s2];
                for(int k=s2; k>0; --k)
                {
                    i[k] = i[k+s2];
                    i[k+s2] = i[k-1];
                }
                i[0] = v;
            }
        }
    }
}

namespace detail
{

inline fftw_plan 
fftwPlanCreate(unsigned int N, int* shape, 
               fftw_complex * in,  int* instrides,  int instep,
               fftw_complex * out, int* outstrides, int outstep,
               int sign, unsigned int planner_flags)
{
    return fftw_plan_many_dft(N, shape, 1,
                              in, instrides, instep, 0,
                              out, outstrides, outstep, 0,
                              sign, planner_flags);
}

inline fftwf_plan 
fftwPlanCreate(unsigned int N, int* shape, 
               fftwf_complex * in,  int* instrides,  int instep,
               fftwf_complex * out, int* outstrides, int outstep,
               int sign, unsigned int planner_flags)
{
    return fftwf_plan_many_dft(N, shape, 1,
                               in, instrides, instep, 0,
                               out, outstrides, outstep, 0,
                               sign, planner_flags);
}

inline fftwl_plan 
fftwPlanCreate(unsigned int N, int* shape, 
               fftwl_complex * in,  int* instrides,  int instep,
               fftwl_complex * out, int* outstrides, int outstep,
               int sign, unsigned int planner_flags)
{
    return fftwl_plan_many_dft(N, shape, 1,
                               in, instrides, instep, 0,
                               out, outstrides, outstep, 0,
                               sign, planner_flags);
}

inline void fftwPlanDestroy(fftw_plan plan)
{
    if(plan != 0)
        fftw_destroy_plan(plan);
}

inline void fftwPlanDestroy(fftwf_plan plan)
{
    if(plan != 0)
        fftwf_destroy_plan(plan);
}

inline void fftwPlanDestroy(fftwl_plan plan)
{
    if(plan != 0)
        fftwl_destroy_plan(plan);
}

inline void 
fftwPlanExecute(fftw_plan plan) 
{
    fftw_execute(plan);
}

inline void 
fftwPlanExecute(fftwf_plan plan) 
{
    fftwf_execute(plan);
}

inline void 
fftwPlanExecute(fftwl_plan plan) 
{
    fftwl_execute(plan);
}

inline void 
fftwPlanExecute(fftw_plan plan, fftw_complex * in,  fftw_complex * out) 
{
    fftw_execute_dft(plan, in, out);
}

inline void 
fftwPlanExecute(fftwf_plan plan, fftwf_complex * in,  fftwf_complex * out) 
{
    fftwf_execute_dft(plan, in, out);
}

inline void 
fftwPlanExecute(fftwl_plan plan, fftwl_complex * in,  fftwl_complex * out) 
{
    fftwl_execute_dft(plan, in, out);
}

} // namespace detail

template <class Real>
class FFTWPlan
{
    typedef ArrayVector<int> Shape;
    typedef typename FFTWReal2Complex<Real>::plan_type PlanType;
    typedef typename FFTWComplex<Real>::complex_type Complex;
    
    PlanType plan;
    Shape shape, instrides, outstrides;
    
  public:
    FFTWPlan()
    : plan(0)
    {}
    
    template <unsigned int N, class C1, class C2>
    FFTWPlan(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
             MultiArrayView<N, FFTWComplex<Real>, C2> out,
             int SIGN, unsigned int planner_flags)
    : plan(0)
    {
        init(in, out, SIGN, planner_flags);
    }
    

    FFTWPlan(FFTWPlan const & other)
    : plan(other.plan)
    {
        FFTWPlan & o = const_cast<FFTWPlan &>(other);
		shape.swap(o.shape);
        instrides.swap(o.instrides);
        outstrides.swap(o.outstrides);
        o.plan = 0; // act like std::auto_ptr
    }
    
    FFTWPlan & operator=(FFTWPlan const & other)
    {
        if(this != &other)
        {
            FFTWPlan & o = const_cast<FFTWPlan &>(other);
			plan = o.plan;
            shape.swap(o.shape);
            instrides.swap(o.instrides);
            outstrides.swap(o.outstrides);
            o.plan = 0; // act like std::auto_ptr
        }
        return *this;
    }

    ~FFTWPlan()
    {
        detail::fftwPlanDestroy(plan);
    }

    template <unsigned int N, class C1, class C2>
    void init(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
              MultiArrayView<N, FFTWComplex<Real>, C2> out,
              int SIGN, unsigned int planner_flags);
        
    template <unsigned int N, class C1, class C2>
    void execute(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                 MultiArrayView<N, FFTWComplex<Real>, C2> out) const;

    void execute() const
	{
		vigra_precondition(plan != 0, "FFTWPlan::execute(): plan is NULL.");
        detail::fftwPlanExecute(plan);
	}
    
  public:
    static int bestPaddingSize(int s)
    {
        // Image sizes where FFTW is fast. The list contains all
        // numbers less than 100000 whose prime decomposition is of the form
        // 2^a*3^b*5^c*7^d*11^e*13^f, where e+f is either 0 or 1, and the 
        // other exponents are arbitrary
        static const int size = 1330;
        static int goodSizes[size] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
            18, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 35, 36, 39, 40, 42, 44, 45, 48, 
            49, 50, 52, 54, 55, 56, 60, 63, 64, 65, 66, 70, 72, 75, 77, 78, 80, 81, 
            84, 88, 90, 91, 96, 98, 99, 100, 104, 105, 108, 110, 112, 117, 120, 125, 
            126, 128, 130, 132, 135, 140, 144, 147, 150, 154, 156, 160, 162, 165, 
            168, 175, 176, 180, 182, 189, 192, 195, 196, 198, 200, 208, 210, 216, 
            220, 224, 225, 231, 234, 240, 243, 245, 250, 252, 256, 260, 264, 270, 
            273, 275, 280, 288, 294, 297, 300, 308, 312, 315, 320, 324, 325, 330, 
            336, 343, 350, 351, 352, 360, 364, 375, 378, 384, 385, 390, 392, 396, 
            400, 405, 416, 420, 432, 440, 441, 448, 450, 455, 462, 468, 480, 486, 
            490, 495, 500, 504, 512, 520, 525, 528, 539, 540, 546, 550, 560, 567, 
            576, 585, 588, 594, 600, 616, 624, 625, 630, 637, 640, 648, 650, 660, 
            672, 675, 686, 693, 700, 702, 704, 720, 728, 729, 735, 750, 756, 768, 
            770, 780, 784, 792, 800, 810, 819, 825, 832, 840, 864, 875, 880, 882, 
            891, 896, 900, 910, 924, 936, 945, 960, 972, 975, 980, 990, 1000, 1008, 
            1024, 1029, 1040, 1050, 1053, 1056, 1078, 1080, 1092, 1100, 1120, 1125, 
            1134, 1152, 1155, 1170, 1176, 1188, 1200, 1215, 1225, 1232, 1248, 1250, 
            1260, 1274, 1280, 1296, 1300, 1320, 1323, 1344, 1350, 1365, 1372, 1375, 
            1386, 1400, 1404, 1408, 1440, 1456, 1458, 1470, 1485, 1500, 1512, 1536, 
            1540, 1560, 1568, 1575, 1584, 1600, 1617, 1620, 1625, 1638, 1650, 1664, 
            1680, 1701, 1715, 1728, 1750, 1755, 1760, 1764, 1782, 1792, 1800, 1820, 
            1848, 1872, 1875, 1890, 1911, 1920, 1925, 1944, 1950, 1960, 1980, 2000, 
            2016, 2025, 2048, 2058, 2079, 2080, 2100, 2106, 2112, 2156, 2160, 2184, 
            2187, 2200, 2205, 2240, 2250, 2268, 2275, 2304, 2310, 2340, 2352, 2376, 
            2400, 2401, 2430, 2450, 2457, 2464, 2475, 2496, 2500, 2520, 2548, 2560, 
            2592, 2600, 2625, 2640, 2646, 2673, 2688, 2695, 2700, 2730, 2744, 2750, 
            2772, 2800, 2808, 2816, 2835, 2880, 2912, 2916, 2925, 2940, 2970, 3000, 
            3024, 3072, 3080, 3087, 3120, 3125, 3136, 3150, 3159, 3168, 3185, 3200, 
            3234, 3240, 3250, 3276, 3300, 3328, 3360, 3375, 3402, 3430, 3456, 3465, 
            3500, 3510, 3520, 3528, 3564, 3584, 3600, 3640, 3645, 3675, 3696, 3744, 
            3750, 3773, 3780, 3822, 3840, 3850, 3888, 3900, 3920, 3960, 3969, 4000, 
            4032, 4050, 4095, 4096, 4116, 4125, 4158, 4160, 4200, 4212, 4224, 4312, 
            4320, 4368, 4374, 4375, 4400, 4410, 4455, 4459, 4480, 4500, 4536, 4550, 
            4608, 4620, 4680, 4704, 4725, 4752, 4800, 4802, 4851, 4860, 4875, 4900, 
            4914, 4928, 4950, 4992, 5000, 5040, 5096, 5103, 5120, 5145, 5184, 5200, 
            5250, 5265, 5280, 5292, 5346, 5376, 5390, 5400, 5460, 5488, 5500, 5544, 
            5600, 5616, 5625, 5632, 5670, 5733, 5760, 5775, 5824, 5832, 5850, 5880, 
            5940, 6000, 6048, 6075, 6125, 6144, 6160, 6174, 6237, 6240, 6250, 6272, 
            6300, 6318, 6336, 6370, 6400, 6468, 6480, 6500, 6552, 6561, 6600, 6615, 
            6656, 6720, 6750, 6804, 6825, 6860, 6875, 6912, 6930, 7000, 7020, 7040, 
            7056, 7128, 7168, 7200, 7203, 7280, 7290, 7350, 7371, 7392, 7425, 7488, 
            7500, 7546, 7560, 7644, 7680, 7700, 7776, 7800, 7840, 7875, 7920, 7938, 
            8000, 8019, 8064, 8085, 8100, 8125, 8190, 8192, 8232, 8250, 8316, 8320, 
            8400, 8424, 8448, 8505, 8575, 8624, 8640, 8736, 8748, 8750, 8775, 8800, 
            8820, 8910, 8918, 8960, 9000, 9072, 9100, 9216, 9240, 9261, 9360, 9375, 
            9408, 9450, 9477, 9504, 9555, 9600, 9604, 9625, 9702, 9720, 9750, 9800, 
            9828, 9856, 9900, 9984, 10000, 10080, 10125, 10192, 10206, 10240, 10290, 
            10368, 10395, 10400, 10500, 10530, 10560, 10584, 10692, 10752, 10780, 
            10800, 10920, 10935, 10976, 11000, 11025, 11088, 11200, 11232, 11250, 
            11264, 11319, 11340, 11375, 11466, 11520, 11550, 11648, 11664, 11700, 
            11760, 11880, 11907, 12000, 12005, 12096, 12150, 12250, 12285, 12288, 
            12320, 12348, 12375, 12474, 12480, 12500, 12544, 12600, 12636, 12672, 
            12740, 12800, 12936, 12960, 13000, 13104, 13122, 13125, 13200, 13230, 
            13312, 13365, 13377, 13440, 13475, 13500, 13608, 13650, 13720, 13750, 
            13824, 13860, 14000, 14040, 14080, 14112, 14175, 14256, 14336, 14400, 
            14406, 14553, 14560, 14580, 14625, 14700, 14742, 14784, 14850, 14976, 
            15000, 15092, 15120, 15288, 15309, 15360, 15400, 15435, 15552, 15600, 
            15625, 15680, 15750, 15795, 15840, 15876, 15925, 16000, 16038, 16128, 
            16170, 16200, 16250, 16380, 16384, 16464, 16500, 16632, 16640, 16800, 
            16807, 16848, 16875, 16896, 17010, 17150, 17199, 17248, 17280, 17325, 
            17472, 17496, 17500, 17550, 17600, 17640, 17820, 17836, 17920, 18000, 
            18144, 18200, 18225, 18375, 18432, 18480, 18522, 18711, 18720, 18750, 
            18816, 18865, 18900, 18954, 19008, 19110, 19200, 19208, 19250, 19404, 
            19440, 19500, 19600, 19656, 19683, 19712, 19800, 19845, 19968, 20000, 
            20160, 20250, 20384, 20412, 20475, 20480, 20580, 20625, 20736, 20790, 
            20800, 21000, 21060, 21120, 21168, 21384, 21504, 21560, 21600, 21609, 
            21840, 21870, 21875, 21952, 22000, 22050, 22113, 22176, 22275, 22295, 
            22400, 22464, 22500, 22528, 22638, 22680, 22750, 22932, 23040, 23100, 
            23296, 23328, 23400, 23520, 23625, 23760, 23814, 24000, 24010, 24057, 
            24192, 24255, 24300, 24375, 24500, 24570, 24576, 24640, 24696, 24750, 
            24948, 24960, 25000, 25088, 25200, 25272, 25344, 25480, 25515, 25600, 
            25725, 25872, 25920, 26000, 26208, 26244, 26250, 26325, 26400, 26411, 
            26460, 26624, 26730, 26754, 26880, 26950, 27000, 27216, 27300, 27440, 
            27500, 27648, 27720, 27783, 28000, 28080, 28125, 28160, 28224, 28350, 
            28431, 28512, 28665, 28672, 28800, 28812, 28875, 29106, 29120, 29160, 
            29250, 29400, 29484, 29568, 29700, 29952, 30000, 30184, 30240, 30375, 
            30576, 30618, 30625, 30720, 30800, 30870, 31104, 31185, 31200, 31213, 
            31250, 31360, 31500, 31590, 31680, 31752, 31850, 32000, 32076, 32256, 
            32340, 32400, 32500, 32760, 32768, 32805, 32928, 33000, 33075, 33264, 
            33280, 33600, 33614, 33696, 33750, 33792, 33957, 34020, 34125, 34300, 
            34375, 34398, 34496, 34560, 34650, 34944, 34992, 35000, 35100, 35200, 
            35280, 35640, 35672, 35721, 35840, 36000, 36015, 36288, 36400, 36450, 
            36750, 36855, 36864, 36960, 37044, 37125, 37422, 37440, 37500, 37632, 
            37730, 37800, 37908, 38016, 38220, 38400, 38416, 38500, 38808, 38880, 
            39000, 39200, 39312, 39366, 39375, 39424, 39600, 39690, 39936, 40000, 
            40095, 40131, 40320, 40425, 40500, 40625, 40768, 40824, 40950, 40960, 
            41160, 41250, 41472, 41580, 41600, 42000, 42120, 42240, 42336, 42525, 
            42768, 42875, 43008, 43120, 43200, 43218, 43659, 43680, 43740, 43750, 
            43875, 43904, 44000, 44100, 44226, 44352, 44550, 44590, 44800, 44928, 
            45000, 45056, 45276, 45360, 45500, 45864, 45927, 46080, 46200, 46305, 
            46592, 46656, 46800, 46875, 47040, 47250, 47385, 47520, 47628, 47775, 
            48000, 48020, 48114, 48125, 48384, 48510, 48600, 48750, 49000, 49140, 
            49152, 49280, 49392, 49500, 49896, 49920, 50000, 50176, 50400, 50421, 
            50544, 50625, 50688, 50960, 51030, 51200, 51450, 51597, 51744, 51840, 
            51975, 52000, 52416, 52488, 52500, 52650, 52800, 52822, 52920, 53248, 
            53460, 53508, 53760, 53900, 54000, 54432, 54600, 54675, 54880, 55000, 
            55125, 55296, 55440, 55566, 56000, 56133, 56160, 56250, 56320, 56448, 
            56595, 56700, 56862, 56875, 57024, 57330, 57344, 57600, 57624, 57750, 
            58212, 58240, 58320, 58500, 58800, 58968, 59049, 59136, 59400, 59535, 
            59904, 60000, 60025, 60368, 60480, 60750, 61152, 61236, 61250, 61425, 
            61440, 61600, 61740, 61875, 62208, 62370, 62400, 62426, 62500, 62720, 
            63000, 63180, 63360, 63504, 63700, 64000, 64152, 64512, 64680, 64800, 
            64827, 65000, 65520, 65536, 65610, 65625, 65856, 66000, 66150, 66339, 
            66528, 66560, 66825, 66885, 67200, 67228, 67375, 67392, 67500, 67584, 
            67914, 68040, 68250, 68600, 68750, 68796, 68992, 69120, 69300, 69888, 
            69984, 70000, 70200, 70400, 70560, 70875, 71280, 71344, 71442, 71680, 
            72000, 72030, 72171, 72576, 72765, 72800, 72900, 73125, 73500, 73710, 
            73728, 73920, 74088, 74250, 74844, 74880, 75000, 75264, 75460, 75600, 
            75816, 76032, 76440, 76545, 76800, 76832, 77000, 77175, 77616, 77760, 
            78000, 78125, 78400, 78624, 78732, 78750, 78848, 78975, 79200, 79233, 
            79380, 79625, 79872, 80000, 80190, 80262, 80640, 80850, 81000, 81250, 
            81536, 81648, 81900, 81920, 82320, 82500, 82944, 83160, 83200, 83349, 
            84000, 84035, 84240, 84375, 84480, 84672, 85050, 85293, 85536, 85750, 
            85995, 86016, 86240, 86400, 86436, 86625, 87318, 87360, 87480, 87500, 
            87750, 87808, 88000, 88200, 88452, 88704, 89100, 89180, 89600, 89856, 
            90000, 90112, 90552, 90720, 91000, 91125, 91728, 91854, 91875, 92160, 
            92400, 92610, 93184, 93312, 93555, 93600, 93639, 93750, 94080, 94325, 
            94500, 94770, 95040, 95256, 95550, 96000, 96040, 96228, 96250, 96768, 
            97020, 97200, 97500, 98000, 98280, 98304, 98415, 98560, 98784, 99000, 
            99225, 99792, 99840 }; 

        if(s <= 0 || s > goodSizes[size-1])
            return s;
        return *std::upper_bound(goodSizes, goodSizes+size, s, std::less_equal<int>());
    }
        
};

template <class Real>
template <unsigned int N, class C1, class C2>
void
FFTWPlan<Real>::init(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                     MultiArrayView<N, FFTWComplex<Real>, C2> out,
                     int SIGN, unsigned int planner_flags)
{
    vigra_precondition(in.shape() == out.shape(),
        "FFTWPlan.init(): input and output must have the same shape.");
    vigra_precondition(in.strideOrdering() == out.strideOrdering(),
        "FFTWPlan.init(): input and output must have the same stride ordering.");

    MultiArrayView<N, FFTWComplex<Real>, StridedArrayTag> ins  = in.permuteStridesDescending();
    MultiArrayView<N, FFTWComplex<Real>, StridedArrayTag> outs = out.permuteStridesDescending();
    
	Shape newShape(ins.shape().begin(), ins.shape().end()),
          newIStrides(ins.stride().begin(), ins.stride().end()),
          newOStrides(outs.stride().begin(), outs.stride().end()),
          itotal(newShape), 
		  ototal(newShape);
    for(int j=1; j<N; ++j)
    {
        itotal[j] = ins.stride(j-1) / ins.stride(j);
        ototal[j] = outs.stride(j-1) / outs.stride(j);
    }
    
    PlanType newPlan = detail::fftwPlanCreate(N, newShape.begin(), 
                                 (Complex*)ins.data(), itotal.begin(), ins.stride(N-1),
                                 (Complex*)outs.data(), ototal.begin(), outs.stride(N-1),
                                  SIGN, planner_flags);
    detail::fftwPlanDestroy(plan);
    plan = newPlan;
    shape.swap(newShape);
    instrides.swap(newIStrides);
    outstrides.swap(newOStrides);
}

template <class Real>
template <unsigned int N, class C1, class C2>
void FFTWPlan<Real>::execute(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                             MultiArrayView<N, FFTWComplex<Real>, C2> out) const
{
    vigra_precondition(plan != 0, "FFTWPlan::execute(): plan is NULL.");
    
    vigra_precondition(shape.size() == N,
        "FFTWPlan::execute(): dimension mismatch between plan and data.");

    vigra_precondition(in.shape() == out.shape(),
        "FFTWPlan::execute(): input and output must have the same shape.");

    MultiArrayView<N, FFTWComplex<Real>, StridedArrayTag> ins  = in.permuteStridesDescending();
    MultiArrayView<N, FFTWComplex<Real>, StridedArrayTag> outs = out.permuteStridesDescending();
    
	vigra_precondition((ins.shape() == TinyVectorView<int, N>(shape.data())),
        "FFTWPlan::execute(): shape mismatch between plan and data.");
    vigra_precondition((ins.stride() == TinyVectorView<int, N>(instrides.data())),
        "FFTWPlan::execute(): strides mismatch between plan and input data.");
    vigra_precondition((outs.stride() == TinyVectorView<int, N>(outstrides.data())),
        "FFTWPlan::execute(): strides mismatch between plan and output data.");

	detail::fftwPlanExecute(plan, (Complex*)ins.data(), (Complex*)outs.data());
}

class FourierOptions
{
  public:
    
    ArrayVector<MultiArrayIndex> padded_shape;
    unsigned int planner_flags;
    bool normalize_inverse;
    
    FourierOptions()
    : planner_flags(FFTW_ESTIMATE),
      normalize_inverse(true)
    {}

	FourierOptions & plannerFlags(unsigned int flags)
	{
		planner_flags = flags;
		return *this;
	}

	FourierOptions & normalizeInverse(bool normalize = true)
	{
		normalize_inverse = normalize;
		return *this;
	}
};

template <unsigned int N, class Real, class C1, class C2>
inline void 
fourierTransform(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                 MultiArrayView<N, FFTWComplex<Real>, C2> out,
                 FourierOptions const & options = FourierOptions())
{
	FFTWPlan<Real>(in, out, FFTW_FORWARD, options.planner_flags).execute();
}

template <unsigned int N, class Real, class C1, class C2>
inline void 
fourierTransformInverse(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                        MultiArrayView<N, FFTWComplex<Real>, C2> out,
                        FourierOptions const & options = FourierOptions())
{
	FFTWPlan<Real>(in, out, FFTW_BACKWARD, options.planner_flags).execute();
    if(options.normalize_inverse)
    {
        Real norm = (Real)out.shape(0);
        for(int j=1; j<N; ++j)
        {
            norm *= (Real)out.shape(j);
        }
        out *= FFTWComplex<Real>(Real(1.0) / norm);
    }
}

template <unsigned int N, class Real, class C1, class C2>
void 
fourierTransform(MultiArrayView<N, Real, C1> in, 
                 MultiArrayView<N, FFTWComplex<Real>, C2> out,
                 FourierOptions const & options = FourierOptions())
{
    vigra_precondition(in.shape() == out.shape(),
        "fourierTransform(): input and output must have the same shape.");
    
    // copy the input array into the output and then perform an in-place FFT
    out = in;
	FFTWPlan<Real>(out, out, FFTW_FORWARD, options.planner_flags).execute();
}

template <unsigned int N, class Real, class C1, class C2>
void 
fourierTransformInverse2(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                         MultiArrayView<N, FFTWComplex<Real>, C2> out,
                         FourierOptions const & options = FourierOptions())
{
    using namespace vigra::functor;
    
    typedef FFTWComplex<Real> C;
    typedef MultiArrayView<1, C, UnstridedArrayTag> MV;
    typedef MultiArrayShape<1>::type S1;
    
    vigra_precondition(in.shape() == out.shape(),
        "fourierTransform2(): input and output must have the same shape.");
    
    MultiArrayIndex maxSize = 0;
    Real norm = (Real)1.0;
    for(unsigned int d=0; d<N; ++d)
    {
        maxSize = std::max(maxSize, in.shape(d));
        norm *= (Real)in.shape(d);
    }
    
    ArrayVector<C> tmp(maxSize);
    C * t = &tmp[0];
    
    ArrayVector<FFTWPlan<Real> > plans(N);
    for(unsigned int d=0; d<N; ++d)
        plans[d].init(MV(S1(in.shape(d)), t), MV(S1(in.shape(d)), t), 
                      FFTW_BACKWARD, options.planner_flags);
   
    typedef typename MultiArrayView<N, C, C2>::traverser Traverser;
    typedef MultiArrayNavigator<Traverser, N> Navigator;
    typedef typename Navigator::iterator Iterator;
    
    Navigator inav(in.traverser_begin(), in.shape(), 0);
    Navigator onav(out.traverser_begin(), out.shape(), 0);

    for( ; inav.hasMore(); inav++, onav++ )
    {
        if(options.normalize_inverse)
            std::transform(inav.begin(), inav.end(), t, Arg1() * Param(C(1.0 / norm)));
        else
            std::copy(inav.begin(), inav.end(), t);
        
        plans[0].execute();
        
        std::copy(t, t+out.shape(0), onav.begin());
    }
    
    for(unsigned int d = 1; d < N; ++d)
    {
        Navigator onav(out.traverser_begin(), out.shape(), d);

        for( ; onav.hasMore(); onav++ )
        {
            std::copy(onav.begin(), onav.end(), t);
            
            plans[d].execute();
            
            std::copy(t, t+out.shape(d), onav.begin());
        }
    }
}

#if 0
namespace detail {

template <class T>
void
fourierTransformImpl(FFTWComplexImage::const_traverser sul,
                     FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                     FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest, T sign)
{
    int w = int(slr.x - sul.x);
    int h = int(slr.y - sul.y);

    FFTWComplexImage sworkImage, dworkImage;

    fftw_complex * srcPtr = (fftw_complex *)(&*sul);
    fftw_complex * destPtr = (fftw_complex *)(&*dul);

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (h > 1 && &(*(sul + Diff2D(w, 0))) != &(*(sul + Diff2D(0, 1))))
    {
        sworkImage.resize(w, h);
        copyImage(srcIterRange(sul, slr, src), destImage(sworkImage));
        srcPtr = (fftw_complex *)(&(*sworkImage.upperLeft()));
    }
    if (h > 1 && &(*(dul + Diff2D(w, 0))) != &(*(dul + Diff2D(0, 1))))
    {
        dworkImage.resize(w, h);
        destPtr = (fftw_complex *)(&(*dworkImage.upperLeft()));
    }

    fftw_plan plan = fftw_plan_dft_2d(h, w, srcPtr, destPtr, sign, FFTW_ESTIMATE );
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if (h > 1 && &(*(dul + Diff2D(w, 0))) != &(*(dul + Diff2D(0, 1))))
    {
        copyImage(srcImageRange(dworkImage), destIter(dul, dest));
    }
}

} // namespace detail

/********************************************************/
/*                                                      */
/*                   fourierTransform                   */
/*                                                      */
/********************************************************/

/** \brief Compute forward and inverse Fourier transforms.

    In the forward direction, the input image may be scalar or complex, and the output image
    is always complex. In the inverse direction, both input and output must be complex.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor>
        void fourierTransform(SrcImageIterator srcUpperLeft,
                              SrcImageIterator srcLowerRight, SrcAccessor src,
                              FFTWComplexImage::traverser destUpperLeft, FFTWComplexImage::Accessor dest);

        void
        fourierTransformInverse(FFTWComplexImage::const_traverser sul,
                                FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                                FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor>
        void fourierTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                              pair<FFTWComplexImage::traverser, FFTWComplexImage::Accessor> dest);

        void
        fourierTransformInverse(triple<FFTWComplexImage::const_traverser,
                                       FFTWComplexImage::const_traverser, FFTWComplexImage::ConstAccessor> src,
                                pair<FFTWComplexImage::traverser, FFTWComplexImage::Accessor> dest);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
    Namespace: vigra

    \code
    // compute complex Fourier transform of a real image
    vigra::DImage src(w, h);
    vigra::FFTWComplexImage fourier(w, h);

    fourierTransform(srcImageRange(src), destImage(fourier));

    // compute inverse Fourier transform
    // note that both source and destination image must be of type vigra::FFTWComplexImage
    vigra::FFTWComplexImage inverseFourier(w, h);

    fourierTransform(srcImageRange(fourier), destImage(inverseFourier));
    \endcode
*/
doxygen_overloaded_function(template <...> void fourierTransform)

inline void
fourierTransform(FFTWComplexImage::const_traverser sul,
                 FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                 FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest)
{
    detail::fourierTransformImpl(sul, slr, src, dul, dest, FFTW_FORWARD);
}

template <class SrcImageIterator, class SrcAccessor>
void fourierTransform(SrcImageIterator srcUpperLeft,
                      SrcImageIterator srcLowerRight, SrcAccessor sa,
                      FFTWComplexImage::traverser destUpperLeft, FFTWComplexImage::Accessor da)
{
    // copy real input images into a complex one...
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    // ...and call the complex -> complex version of the algorithm
    FFTWComplexImage const & cworkImage = workImage;
    fourierTransform(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                     destUpperLeft, da);
}

template <class SrcImageIterator, class SrcAccessor>
inline
void fourierTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                      pair<FFTWComplexImage::traverser, FFTWComplexImage::Accessor> dest)
{
    fourierTransform(src.first, src.second, src.third, dest.first, dest.second);
}

/** \brief Compute inverse Fourier transforms.

    See \ref fourierTransform() for details.
*/
inline void
fourierTransformInverse(FFTWComplexImage::const_traverser sul,
                        FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                        FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest)
{
    detail::fourierTransformImpl(sul, slr, src, dul, dest, FFTW_BACKWARD);
}

template <class DestImageIterator, class DestAccessor>
void fourierTransformInverse(FFTWComplexImage::const_traverser sul,
                             FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                             DestImageIterator dul, DestAccessor dest)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    FFTWComplexImage workImage(w, h);
    fourierTransformInverse(sul, slr, src, workImage.upperLeft(), workImage.accessor());
    copyImage(srcImageRange(workImage), destIter(dul, dest));
}


template <class DestImageIterator, class DestAccessor>
inline void
fourierTransformInverse(triple<FFTWComplexImage::const_traverser,
                               FFTWComplexImage::const_traverser, FFTWComplexImage::ConstAccessor> src,
                        pair<DestImageIterator, DestAccessor> dest)
{
    fourierTransformInverse(src.first, src.second, src.third, dest.first, dest.second);
}



/********************************************************/
/*                                                      */
/*                   applyFourierFilter                 */
/*                                                      */
/********************************************************/

/** \brief Apply a filter (defined in the frequency domain) to an image.

    After transferring the image into the frequency domain, it is
    multiplied pixel-wise with the filter and transformed back. The
    result is put into the given destination image which must have the right size.
    The result will be normalized to compensate for the two FFTs.

    If the destination image is scalar, only the real part of the result image is
    retained. In this case, you are responsible for choosing a filter image
    which ensures a zero imaginary part of the result (e.g. use a real, even symmetric
    filter image, or a purely imaginary, odd symmetric on).

    The DC entry of the filter must be in the upper left, which is the
    position where FFTW expects it (see \ref moveDCToUpperLeft()).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class FilterImageIterator, class FilterAccessor,
                  class DestImageIterator, class DestAccessor>
        void applyFourierFilter(SrcImageIterator srcUpperLeft,
                                SrcImageIterator srcLowerRight, SrcAccessor sa,
                                FilterImageIterator filterUpperLeft, FilterAccessor fa,
                                DestImageIterator destUpperLeft, DestAccessor da);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class FilterImageIterator, class FilterAccessor,
                  class DestImageIterator, class DestAccessor>
        void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                pair<FilterImageIterator, FilterAccessor> filter,
                                pair<DestImageIterator, DestAccessor> dest);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
    Namespace: vigra

    \code
    // create a Gaussian filter in Fourier space
    vigra::FImage gaussFilter(w, h), filter(w, h);
    for(int y=0; y<h; ++y)
        for(int x=0; x<w; ++x)
        {
            xx = float(x - w / 2) / w;
            yy = float(y - h / 2) / h;

            gaussFilter(x,y) = std::exp(-(xx*xx + yy*yy) / 2.0 * scale);
        }

    // applyFourierFilter() expects the filter's DC in the upper left
    moveDCToUpperLeft(srcImageRange(gaussFilter), destImage(filter));

    vigra::FFTWComplexImage result(w, h);

    vigra::applyFourierFilter(srcImageRange(image), srcImage(filter), result);
    \endcode

    For inspection of the result, \ref FFTWMagnitudeAccessor might be
    useful. If you want to apply the same filter repeatedly, it may be more
    efficient to use the FFTW functions directly with FFTW plans optimized
    for good performance.
*/
doxygen_overloaded_function(template <...> void applyFourierFilter)

template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilter(SrcImageIterator srcUpperLeft,
                        SrcImageIterator srcLowerRight, SrcAccessor sa,
                        FilterImageIterator filterUpperLeft, FilterAccessor fa,
                        DestImageIterator destUpperLeft, DestAccessor da)
{
    // copy real input images into a complex one...
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    // ...and call the impl
    FFTWComplexImage const & cworkImage = workImage;
    applyFourierFilterImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                           filterUpperLeft, fa,
                           destUpperLeft, da);
}

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
inline
void applyFourierFilter(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da)
{
    int w = srcLowerRight.x - srcUpperLeft.x;
    int h = srcLowerRight.y - srcUpperLeft.y;

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (&(*(srcUpperLeft + Diff2D(w, 0))) == &(*(srcUpperLeft + Diff2D(0, 1))))
        applyFourierFilterImpl(srcUpperLeft, srcLowerRight, sa,
                               filterUpperLeft, fa,
                               destUpperLeft, da);
    else
    {
        FFTWComplexImage workImage(w, h);
        copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
                  destImage(workImage));

        FFTWComplexImage const & cworkImage = workImage;
        applyFourierFilterImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                               filterUpperLeft, fa,
                               destUpperLeft, da);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
inline
void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                        pair<FilterImageIterator, FilterAccessor> filter,
                        pair<DestImageIterator, DestAccessor> dest)
{
    applyFourierFilter(src.first, src.second, src.third,
                       filter.first, filter.second,
                       dest.first, dest.second);
}

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilterImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da)
{
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

    FFTWComplexImage complexResultImg(srcLowerRight - srcUpperLeft);

    // FFT from srcImage to complexResultImg
    fftw_plan forwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)&(*srcUpperLeft),
                               (fftw_complex *)complexResultImg.begin(),
                               FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(forwardPlan);
    fftw_destroy_plan(forwardPlan);

    // convolve in freq. domain (in complexResultImg)
    combineTwoImages(srcImageRange(complexResultImg), srcIter(filterUpperLeft, fa),
                     destImage(complexResultImg), std::multiplies<FFTWComplex<> >());

    // FFT back into spatial domain (inplace in complexResultImg)
    fftw_plan backwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)complexResultImg.begin(),
                               (fftw_complex *)complexResultImg.begin(),
                               FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(backwardPlan);
    fftw_destroy_plan(backwardPlan);

    typedef typename
        NumericTraits<typename DestAccessor::value_type>::isScalar
        isScalarResult;

    // normalization (after FFTs), maybe stripping imaginary part
    applyFourierFilterImplNormalization(complexResultImg, destUpperLeft, da,
                                        isScalarResult());
}

template <class DestImageIterator, class DestAccessor>
void applyFourierFilterImplNormalization(FFTWComplexImage const &srcImage,
                                         DestImageIterator destUpperLeft,
                                         DestAccessor da,
                                         VigraFalseType)
{
    double normFactor= 1.0/(srcImage.width() * srcImage.height());

    for(int y=0; y<srcImage.height(); y++, destUpperLeft.y++)
    {
        DestImageIterator dIt= destUpperLeft;
        for(int x= 0; x< srcImage.width(); x++, dIt.x++)
        {
            da.setComponent(srcImage(x, y).re()*normFactor, dIt, 0);
            da.setComponent(srcImage(x, y).im()*normFactor, dIt, 1);
        }
    }
}

inline
void applyFourierFilterImplNormalization(FFTWComplexImage const & srcImage,
        FFTWComplexImage::traverser destUpperLeft,
        FFTWComplexImage::Accessor da,
        VigraFalseType)
{
    transformImage(srcImageRange(srcImage), destIter(destUpperLeft, da),
                   linearIntensityTransform<FFTWComplex<> >(1.0/(srcImage.width() * srcImage.height())));
}

template <class DestImageIterator, class DestAccessor>
void applyFourierFilterImplNormalization(FFTWComplexImage const & srcImage,
                                         DestImageIterator destUpperLeft,
                                         DestAccessor da,
                                         VigraTrueType)
{
    double normFactor= 1.0/(srcImage.width() * srcImage.height());

    for(int y=0; y<srcImage.height(); y++, destUpperLeft.y++)
    {
        DestImageIterator dIt= destUpperLeft;
        for(int x= 0; x< srcImage.width(); x++, dIt.x++)
            da.set(srcImage(x, y).re()*normFactor, dIt);
    }
}

/**********************************************************/
/*                                                        */
/*                applyFourierFilterFamily                */
/*                                                        */
/**********************************************************/

/** \brief Apply an array of filters (defined in the frequency domain) to an image.

    This provides the same functionality as \ref applyFourierFilter(),
    but applying several filters at once allows to avoid
    repeated Fourier transforms of the source image.

    Filters and result images must be stored in \ref vigra::ImageArray data
    structures. In contrast to \ref applyFourierFilter(), this function adjusts
    the size of the result images and the the length of the array.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                                      SrcImageIterator srcLowerRight, SrcAccessor sa,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        void applyFourierFilterFamily(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
    Namespace: vigra

    \code
    // assuming the presence of a real-valued image named "image" to
    // be filtered in this example

    vigra::ImageArray<vigra::FImage> filters(16, image.size());

    for (int i=0; i<filters.size(); i++)
         // create some meaningful filters here
         createMyFilterOfScale(i, destImage(filters[i]));

    vigra::ImageArray<vigra::FFTWComplexImage> results();

    vigra::applyFourierFilterFamily(srcImageRange(image), filters, results);
    \endcode
*/
doxygen_overloaded_function(template <...> void applyFourierFilterFamily)

template <class SrcImageIterator, class SrcAccessor,
          class FilterType, class DestImage>
inline
void applyFourierFilterFamily(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                              const ImageArray<FilterType> &filters,
                              ImageArray<DestImage> &results)
{
    applyFourierFilterFamily(src.first, src.second, src.third,
                             filters, results);
}

template <class SrcImageIterator, class SrcAccessor,
          class FilterType, class DestImage>
void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                              SrcImageIterator srcLowerRight, SrcAccessor sa,
                              const ImageArray<FilterType> &filters,
                              ImageArray<DestImage> &results)
{
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    FFTWComplexImage const & cworkImage = workImage;
    applyFourierFilterFamilyImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                                 filters, results);
}

template <class FilterType, class DestImage>
inline
void applyFourierFilterFamily(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results)
{
    int w= srcLowerRight.x - srcUpperLeft.x;

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (&(*(srcUpperLeft + Diff2D(w, 0))) == &(*(srcUpperLeft + Diff2D(0, 1))))
        applyFourierFilterFamilyImpl(srcUpperLeft, srcLowerRight, sa,
                                     filters, results);
    else
    {
        int h = srcLowerRight.y - srcUpperLeft.y;
        FFTWComplexImage workImage(w, h);
        copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
                  destImage(workImage));

        FFTWComplexImage const & cworkImage = workImage;
        applyFourierFilterFamilyImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                                     filters, results);
    }
}

template <class FilterType, class DestImage>
void applyFourierFilterFamilyImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results)
{
    // FIXME: sa is not used
    // (maybe check if StandardAccessor, else copy?)    

    // make sure the filter images have the right dimensions
    vigra_precondition((srcLowerRight - srcUpperLeft) == filters.imageSize(),
                       "applyFourierFilterFamily called with src image size != filters.imageSize()!");

    // make sure the result image array has the right dimensions
    results.resize(filters.size());
    results.resizeImages(filters.imageSize());

    // FFT from srcImage to freqImage
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

    FFTWComplexImage freqImage(w, h);
    FFTWComplexImage result(w, h);

    fftw_plan forwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)&(*srcUpperLeft),
                               (fftw_complex *)freqImage.begin(),
                               FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(forwardPlan);
    fftw_destroy_plan(forwardPlan);

    fftw_plan backwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)result.begin(),
                               (fftw_complex *)result.begin(),
                               FFTW_BACKWARD, FFTW_ESTIMATE );
    typedef typename
        NumericTraits<typename DestImage::Accessor::value_type>::isScalar
        isScalarResult;

    // convolve with filters in freq. domain
    for (unsigned int i= 0;  i < filters.size(); i++)
    {
        combineTwoImages(srcImageRange(freqImage), srcImage(filters[i]),
                         destImage(result), std::multiplies<FFTWComplex<> >());

        // FFT back into spatial domain (inplace in destImage)
        fftw_execute(backwardPlan);

        // normalization (after FFTs), maybe stripping imaginary part
        applyFourierFilterImplNormalization(result,
                                            results[i].upperLeft(), results[i].accessor(),
                                            isScalarResult());
    }
    fftw_destroy_plan(backwardPlan);
}

/********************************************************/
/*                                                      */
/*                fourierTransformReal                  */
/*                                                      */
/********************************************************/

/** \brief Real Fourier transforms for even and odd boundary conditions
           (aka. cosine and sine transforms).


    If the image is real and has even symmetry, its Fourier transform
    is also real and has even symmetry. The Fourier transform of a real image with odd
    symmetry is imaginary and has odd symmetry. In either case, only about a quarter
    of the pixels need to be stored because the rest can be calculated from the symmetry
    properties. This is especially useful, if the original image is implicitly assumed
    to have reflective or anti-reflective boundary conditions. Then the "negative"
    pixel locations are defined as

    \code
    even (reflective boundary conditions):      f[-x] = f[x]     (x = 1,...,N-1)
    odd (anti-reflective boundary conditions):  f[-1] = 0
                                                f[-x] = -f[x-2]  (x = 2,...,N-1)
    \endcode

    end similar at the other boundary (see the FFTW documentation for details).
    This has the advantage that more efficient Fourier transforms that use only
    real numbers can be implemented. These are also known as cosine and sine transforms
    respectively.

    If you use the odd transform it is important to note that in the Fourier domain,
    the DC component is always zero and is therefore dropped from the data structure.
    This means that index 0 in an odd symmetric Fourier domain image refers to
    the <i>first</i> harmonic. This is especially important if an image is first
    cosine transformed (even symmetry), then in the Fourier domain multiplied
    with an odd symmetric filter (e.g. a first derivative) and finally transformed
    back to the spatial domain with a sine transform (odd symmetric). For this to work
    properly the image must be shifted left or up by one pixel (depending on whether
    the x- or y-axis is odd symmetric) before the inverse transform can be applied.
    (see example below).

    The real Fourier transform functions are named <tt>fourierTransformReal??</tt>
    where the questions marks stand for either <tt>E</tt> or <tt>O</tt> indicating
    whether the x- and y-axis is to be transformed using even or odd symmetry.
    The same functions can be used for both the forward and inverse transforms,
    only the normalization changes. For signal processing, the following
    normalization factors are most appropriate:

    \code
                          forward             inverse
    ------------------------------------------------------------
    X even, Y even           1.0         4.0 * (w-1) * (h-1)
    X even, Y odd           -1.0        -4.0 * (w-1) * (h+1)
    X odd,  Y even          -1.0        -4.0 * (w+1) * (h-1)
    X odd,  Y odd            1.0         4.0 * (w+1) * (h+1)
    \endcode

    where <tt>w</tt> and <tt>h</tt> denote the image width and height.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcTraverser, class SrcAccessor,
                  class DestTraverser, class DestAccessor>
        void
        fourierTransformRealEE(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm);

        fourierTransformRealEO, fourierTransformRealOE, fourierTransformRealOO likewise
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcTraverser, class SrcAccessor,
                  class DestTraverser, class DestAccessor>
        void
        fourierTransformRealEE(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm);

        fourierTransformRealEO, fourierTransformRealOE, fourierTransformRealOO likewise
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> \<vigra/fftw3.hxx\><br>
        Namespace: vigra

    \code
    vigra::FImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // forward cosine transform == reflective boundary conditions
    fourierTransformRealEE(srcImageRange(spatial), destImage(fourier), (fftw_real)1.0);

    // multiply with a first derivative of Gaussian in x-direction
    for(int y = 0; y < height; ++y)
    {
        for(int x = 1; x < width; ++x)
        {
            double dx = x * M_PI / (width - 1);
            double dy = y * M_PI / (height - 1);
            fourier(x-1, y) = fourier(x, y) * dx * std::exp(-(dx*dx + dy*dy) * scale*scale / 2.0);
        }
        fourier(width-1, y) = 0.0;
    }

    // inverse transform -- odd symmetry in x-direction, even in y,
    //                      due to symmetry of the filter
    fourierTransformRealOE(srcImageRange(fourier), destImage(spatial),
                           (fftw_real)-4.0 * (width+1) * (height-1));
    \endcode
*/
doxygen_overloaded_function(template <...> void fourierTransformReal)

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEE(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealEE(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEE(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_REDFT00, FFTW_REDFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEE(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_REDFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_REDFT00);
}

/********************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOE(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealOE(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOE(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_RODFT00, FFTW_REDFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOE(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_REDFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_REDFT00);
}

/********************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEO(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealEO(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEO(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_REDFT00, FFTW_RODFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEO(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_RODFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_RODFT00);
}

/********************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOO(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealOO(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOO(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_RODFT00, FFTW_RODFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOO(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_RODFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_RODFT00);
}

/*******************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
void
fourierTransformRealWorkImageImpl(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                                  DestTraverser dul, DestAccessor dest,
                                  fftw_real norm, fftw_r2r_kind kindx, fftw_r2r_kind kindy)
{
    FFTWRealImage workImage(slr - sul);
    copyImage(srcIterRange(sul, slr, src), destImage(workImage));
    FFTWRealImage const & cworkImage = workImage;
    fourierTransformRealImpl(cworkImage.upperLeft(), cworkImage.lowerRight(),
                             dul, dest, norm, kindx, kindy);
}


template <class DestTraverser, class DestAccessor>
void
fourierTransformRealImpl(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         DestTraverser dul, DestAccessor dest,
         fftw_real norm, fftw_r2r_kind kindx, fftw_r2r_kind kindy)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    BasicImage<fftw_real> res(w, h);

    fftw_plan plan = fftw_plan_r2r_2d(h, w,
                         (fftw_real *)&(*sul), (fftw_real *)res.begin(),
                         kindy, kindx, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if(norm != 1.0)
        transformImage(srcImageRange(res), destIter(dul, dest),
                       std::bind1st(std::multiplies<fftw_real>(), 1.0 / norm));
    else
        copyImage(srcImageRange(res), destIter(dul, dest));
}


//@}
#endif

} // namespace vigra

#endif // VIGRA_MULTI_FFT_HXX
