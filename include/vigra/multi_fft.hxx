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

class FFTWPlan
{
    enum PlanType { NONE, FLOAT, DOUBLE, LONGDOUBLE };

public:
    typedef ArrayVector<int> Shape;

    void * plan;
    PlanType plan_type;
    Shape shape, instrides, outstrides;
    
    template <unsigned int N, class Real, class C1, class C2>
    FFTWPlan(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
             MultiArrayView<N, FFTWComplex<Real>, C2> out,
             int SIGN, unsigned int planner_flags,
             bool execute_immediately = false);

	template <class Real>
	void destroy();
    
	template <>
	void destroy<double>()
	{
		if(plan == 0)
			return;
		vigra_precondition(plan_type == DOUBLE,
			"FFTWPlan::destroy<double>(): plan has wrong type.");
		fftw_destroy_plan((fftw_plan)plan);
		plan = 0;
		plan_type = NONE;
	}
    
	template <>
	void destroy<float>()
	{
		if(plan == 0)
			return;
		vigra_precondition(plan_type == FLOAT,
			"FFTWPlan::destroy<float>(): plan has wrong type.");
		fftwf_destroy_plan((fftwf_plan)plan);
		plan = 0;
		plan_type = NONE;
	}
    
	template <>
	void destroy<long double>()
	{
		if(plan == 0)
			return;
		vigra_precondition(plan_type == LONGDOUBLE,
			"FFTWPlan::destroy<long double>(): plan has wrong type.");
		fftwl_destroy_plan((fftwl_plan)plan);
		plan = 0;
		plan_type = NONE;
	}
    
    ~FFTWPlan()
    {
        switch(plan_type)
        {
            case DOUBLE:     fftw_destroy_plan((fftw_plan)plan);   break;
//            case FLOAT:      fftwf_destroy_plan((fftwf_plan)plan); break;
//            case LONGDOUBLE: fftwl_destroy_plan((fftwl_plan)plan); break;
        }
    }

    template <unsigned int N, class Real, class C1, class C2>
    void execute(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                 MultiArrayView<N, FFTWComplex<Real>, C2> out) const;

  private:
  
    FFTWPlan(FFTWPlan const &);
    FFTWPlan & operator=(FFTWPlan const &);
    
    void createPlan(unsigned int N, int* shape, 
                    fftw_complex * in,  int* instrides,  int instep,
                    fftw_complex * out, int* outstrides, int outstep,
                    int sign, unsigned int planner_flags)
    {
        plan = fftw_plan_many_dft(N, shape, 1,
                                  in, instrides, instep, 0,
                                  out, outstrides, outstep, 0,
                                  sign, planner_flags);
        vigra_postcondition(plan != 0,
            "FFTWPlan: Unable to create plan.");
        plan_type = DOUBLE;
    }

    void createPlan(unsigned int N, int* shape, 
                    fftwf_complex * in,  int* instrides,  int instep,
                    fftwf_complex * out, int* outstrides, int outstep,
                    int sign, unsigned int planner_flags)
    {
        plan = fftwf_plan_many_dft(N, shape, 1,
                                   in, instrides, instep, 0,
                                   out, outstrides, outstep, 0,
                                   sign, planner_flags);
        vigra_postcondition(plan != 0,
            "FFTWPlan: Unable to create plan.");
        plan_type = FLOAT;
    }

    void createPlan(unsigned int N, int* shape, 
                    fftwl_complex * in,  int* instrides,  int instep,
                    fftwl_complex * out, int* outstrides, int outstep,
                    int sign, unsigned int planner_flags)
    {
        plan = fftwl_plan_many_dft(N, shape, 1,
                                   in, instrides, instep, 0,
                                   out, outstrides, outstep, 0,
                                   sign, planner_flags);
        vigra_postcondition(plan != 0,
            "FFTWPlan: Unable to create plan.");
        plan_type = LONGDOUBLE;
    }
    
    void executePlan(fftw_complex * in,  fftw_complex * out) const
    {
        vigra_precondition(plan_type == DOUBLE,
            "FFTWPlan::execute(): type mismatch between plan and data.");
        fftw_execute_dft((fftw_plan)plan, in, out);
    }
    
    void executePlan(fftwf_complex * in,  fftwf_complex * out) const
    {
        vigra_precondition(plan_type == FLOAT,
            "FFTWPlan::execute(): type mismatch between plan and data.");
        fftwf_execute_dft((fftwf_plan)plan, in, out);
    }
    
    void executePlan(fftwl_complex * in,  fftwl_complex * out) const
    {
        vigra_precondition(plan_type == LONGDOUBLE,
            "FFTWPlan::execute(): type mismatch between plan and data.");
        fftwl_execute_dft((fftwl_plan)plan, in, out);
    }
    
    void executePlanImmediately(fftw_complex * in,  fftw_complex * out) const
    {
        fftw_execute((fftw_plan)plan);
    }
    
    void executePlanImmediately(fftwf_complex * in,  fftwf_complex * out) const
    {
        fftwf_execute((fftwf_plan)plan);
    }
    
    void executePlanImmediately(fftwl_complex * in,  fftwl_complex * out) const
    {
        fftwl_execute((fftwl_plan)plan);
    }
    
  public:
    static int bestPaddingSize(int s)
    {
        // numbers below 30000 whose prime decomposition only contains powers of 2, 3, 5, 7
        static int goodSizes[482] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 
            16, 18, 20, 21, 24, 25, 27, 28, 30, 32, 35, 36, 40, 42, 45, 48, 49, 50, 
            54, 56, 60, 63, 64, 70, 72, 75, 80, 81, 84, 90, 96, 98, 100, 105, 108, 
            112, 120, 125, 126, 128, 135, 140, 144, 147, 150, 160, 162, 168, 175, 
            180, 189, 192, 196, 200, 210, 216, 224, 225, 240, 243, 245, 250, 252, 
            256, 270, 280, 288, 294, 300, 315, 320, 324, 336, 343, 350, 360, 375, 
            378, 384, 392, 400, 405, 420, 432, 441, 448, 450, 480, 486, 490, 500, 
            504, 512, 525, 540, 560, 567, 576, 588, 600, 625, 630, 640, 648, 672, 
            675, 686, 700, 720, 729, 735, 750, 756, 768, 784, 800, 810, 840, 864, 
            875, 882, 896, 900, 945, 960, 972, 980, 1000, 1008, 1024, 1029, 1050, 
            1080, 1120, 1125, 1134, 1152, 1176, 1200, 1215, 1225, 1250, 1260, 1280, 
            1296, 1323, 1344, 1350, 1372, 1400, 1440, 1458, 1470, 1500, 1512, 1536, 
            1568, 1575, 1600, 1620, 1680, 1701, 1715, 1728, 1750, 1764, 1792, 1800, 
            1875, 1890, 1920, 1944, 1960, 2000, 2016, 2025, 2048, 2058, 2100, 2160, 
            2187, 2205, 2240, 2250, 2268, 2304, 2352, 2400, 2401, 2430, 2450, 2500, 
            2520, 2560, 2592, 2625, 2646, 2688, 2700, 2744, 2800, 2835, 2880, 2916, 
            2940, 3000, 3024, 3072, 3087, 3125, 3136, 3150, 3200, 3240, 3360, 3375, 
            3402, 3430, 3456, 3500, 3528, 3584, 3600, 3645, 3675, 3750, 3780, 3840, 
            3888, 3920, 3969, 4000, 4032, 4050, 4096, 4116, 4200, 4320, 4374, 4375, 
            4410, 4480, 4500, 4536, 4608, 4704, 4725, 4800, 4802, 4860, 4900, 5000, 
            5040, 5103, 5120, 5145, 5184, 5250, 5292, 5376, 5400, 5488, 5600, 5625, 
            5670, 5760, 5832, 5880, 6000, 6048, 6075, 6125, 6144, 6174, 6250, 6272, 
            6300, 6400, 6480, 6561, 6615, 6720, 6750, 6804, 6860, 6912, 7000, 7056, 
            7168, 7200, 7203, 7290, 7350, 7500, 7560, 7680, 7776, 7840, 7875, 7938, 
            8000, 8064, 8100, 8192, 8232, 8400, 8505, 8575, 8640, 8748, 8750, 8820, 
            8960, 9000, 9072, 9216, 9261, 9375, 9408, 9450, 9600, 9604, 9720, 9800, 
            10000, 10080, 10125, 10206, 10240, 10290, 10368, 10500, 10584, 10752, 
            10800, 10935, 10976, 11025, 11200, 11250, 11340, 11520, 11664, 11760, 
            11907, 12000, 12005, 12096, 12150, 12250, 12288, 12348, 12500, 12544, 
            12600, 12800, 12960, 13122, 13125, 13230, 13440, 13500, 13608, 13720, 
            13824, 14000, 14112, 14175, 14336, 14400, 14406, 14580, 14700, 15000, 
            15120, 15309, 15360, 15435, 15552, 15625, 15680, 15750, 15876, 16000, 
            16128, 16200, 16384, 16464, 16800, 16807, 16875, 17010, 17150, 17280, 
            17496, 17500, 17640, 17920, 18000, 18144, 18225, 18375, 18432, 18522, 
            18750, 18816, 18900, 19200, 19208, 19440, 19600, 19683, 19845, 20000, 
            20160, 20250, 20412, 20480, 20580, 20736, 21000, 21168, 21504, 21600, 
            21609, 21870, 21875, 21952, 22050, 22400, 22500, 22680, 23040, 23328, 
            23520, 23625, 23814, 24000, 24010, 24192, 24300, 24500, 24576, 24696, 
            25000, 25088, 25200, 25515, 25600, 25725, 25920, 26244, 26250, 26460, 
            26880, 27000, 27216, 27440, 27648, 27783, 28000, 28125, 28224, 28350, 
            28672, 28800, 28812, 29160, 29400 }; 


        if(s <= 0 || s > 29400)
            return s;
        return *std::upper_bound(goodSizes, goodSizes+482, s, std::less_equal<int>());
    }
        
};

template <unsigned int N, class Real, class C1, class C2>
FFTWPlan::FFTWPlan(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                   MultiArrayView<N, FFTWComplex<Real>, C2> out,
                   int SIGN, unsigned int planner_flags,
                   bool execute_immediately)
: plan(0),
  plan_type(NONE),
  shape(N), instrides(N), outstrides(N)
{
    typedef typename FFTWComplex<Real>::complex_type Complex;
    
    vigra_precondition(in.shape() == out.shape(),
        "FFTWPlan: input and output must have the same shape.");
    vigra_precondition(in.strideOrdering() == out.strideOrdering(),
        "FFTWPlan: input and output must have the same stride ordering.");

    MultiArrayView<N, FFTWComplex<Real>, StridedArrayTag> ins  = in.permuteStridesDescending();
    MultiArrayView<N, FFTWComplex<Real>, StridedArrayTag> outs = out.permuteStridesDescending();
    
	std::copy(ins.shape().begin(), ins.shape().end(), shape.begin());
	std::copy(ins.stride().begin(), ins.stride().end(), instrides.begin());
	std::copy(outs.stride().begin(), outs.stride().end(), outstrides.begin());
    
    Shape itotal(ins.shape().begin(), ins.shape().end()), 
		  ototal(ins.shape().begin(), ins.shape().end());
    for(int j=1; j<N; ++j)
    {
        itotal[j] = ins.stride(j-1) / ins.stride(j);
        ototal[j] = outs.stride(j-1) / outs.stride(j);
    }
    
    createPlan(N, shape.begin(), 
               (Complex*)ins.data(), itotal.begin(), ins.stride(N-1),
               (Complex*)outs.data(), ototal.begin(), outs.stride(N-1),
               SIGN, planner_flags);
    if(execute_immediately)
        executePlanImmediately((Complex*)ins.data(), (Complex*)outs.data());
}

template <unsigned int N, class Real, class C1, class C2>
void FFTWPlan::execute(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                       MultiArrayView<N, FFTWComplex<Real>, C2> out) const
{
    typedef typename FFTWComplex<Real>::complex_type Complex;
    
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

	executePlan((Complex*)ins.data(), (Complex*)outs.data());
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
    FFTWPlan plan(in, out, FFTW_FORWARD, options.planner_flags, /*execute_immediately*/true);
	plan.destroy<Real>();
}

template <unsigned int N, class Real, class C1, class C2>
inline void 
fourierTransformInverse(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                        MultiArrayView<N, FFTWComplex<Real>, C2> out,
                        FourierOptions const & options = FourierOptions())
{
    FFTWPlan plan(in, out, FFTW_BACKWARD, options.planner_flags, /*execute_immediately*/true);
	plan.destroy<Real>();
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
    FFTWPlan plan(out, out, FFTW_FORWARD, options.planner_flags, /*execute_immediately*/true);
	plan.destroy<Real>();
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
