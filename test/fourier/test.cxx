/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
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

#include "unittest.hxx"
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/functorexpression.hxx>
#include <vigra/fftw3.hxx>
#include <vigra/impex.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/gaborfilter.hxx>
#include <vigra/multi_fft.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/convolution.hxx>
#include "test.hxx"

using namespace std;
using namespace vigra;
using namespace vigra::functor;

struct Compare1
{
    double s;

    typedef FFTWComplexImage::value_type value_type;
    typedef FFTWComplexImage::value_type result_type;

    Compare1(double is)
    : s(is)
    {}

    result_type operator()(value_type v1, value_type v2) const
    { return v1 - v2/s; }
};

struct FFTWComplexTest
{
    FFTWComplex<> clx0, clx1, clx2, clx3;

    template <class VECTOR>
    void printVector(VECTOR const & v)
    {
        std::cerr << "(";
        for(int i=0; i<v.size(); ++i)
            std::cerr << v[i] << ", ";
        std::cerr << ")\n";
    }

    FFTWComplexTest()
    : clx0(), clx1(2.0), clx2(0.0, -2.0), clx3(2.0, -2.0)
    {}

    void testConstruction()
    {
        shouldEqual(clx0.re() , 0.0);
        shouldEqual(clx0.im() , 0.0);
        shouldEqual(clx1.re() , 2.0);
        shouldEqual(clx1.im() , 0.0);
        shouldEqual(clx2.re() , 0.0);
        shouldEqual(clx2.im() , -2.0);
        shouldEqual(clx3.re() , 2.0);
        shouldEqual(clx3.im() , -2.0);
        shouldEqual(clx0[0] , 0.0);
        shouldEqual(clx0[1] , 0.0);
        shouldEqual(clx1[0] , 2.0);
        shouldEqual(clx1[1] , 0.0);
        shouldEqual(clx2[0] , 0.0);
        shouldEqual(clx2[1] , -2.0);
        shouldEqual(clx3[0] , 2.0);
        shouldEqual(clx3[1] , -2.0);
    }

    void testComparison()
    {
        should(clx0 == FFTWComplex<>(0.0, 0.0));
        should(clx1 == FFTWComplex<>(2.0, 0.0));
        should(clx2 == FFTWComplex<>(0.0, -2.0));
        should(clx3 == FFTWComplex<>(2.0, -2.0));
        should(clx0 != clx1);
        should(clx0 != clx2);
        should(clx1 != clx2);
        should(clx1 != clx3);
    }

    void testArithmetic()
    {
        shouldEqual(abs(clx0) , 0.0);
        shouldEqual(abs(clx1) , 2.0);
        shouldEqual(abs(clx2) , 2.0);
        shouldEqual(abs(clx3) , sqrt(8.0));
        should(conj(clx3) == FFTWComplex<>(2.0, 2.0));

        shouldEqual(clx1.phase() , 0.0);
        shouldEqual(sin(clx2.phase()) , -1.0);
        shouldEqualTolerance(sin(clx3.phase()), -sqrt(2.0)/2.0, 1.0e-7);

        should(FFTWComplex<>(2.0, -2.0) == clx1 + clx2);

        should(clx0 == clx1 - clx1);
        should(clx0 == clx2 - clx2);
        should(FFTWComplex<>(2.0, 2.0) == clx1 - clx2);

        should(2.0*clx1 == FFTWComplex<>(4.0, 0.0));
        should(2.0*clx2 == FFTWComplex<>(0.0, -4.0));
        should(clx1*clx2 == FFTWComplex<>(0.0, -4.0));
        should(clx3*conj(clx3) == clx3.squaredMagnitude());

        should(clx1/2.0 == FFTWComplex<>(1.0, 0.0));
        should(clx2/2.0 == FFTWComplex<>(0.0, -1.0));
        should(clx2/clx1 == FFTWComplex<>(0.0, -1.0));
        should(clx3*conj(clx3) == clx3.squaredMagnitude());

    }

    void testAccessor()
    {
        FFTWComplexImage img(2,2);
        img = clx2;
        img(0,0) = clx3;
        FFTWComplexImage::Iterator i = img.upperLeft();

        FFTWComplexImage::Accessor get = img.accessor();
        FFTWRealAccessor<> real;
        FFTWImaginaryAccessor<> imag;
        FFTWMagnitudeAccessor<> mag;
        FFTWPhaseAccessor<> phase;
        FFTWWriteRealAccessor<> writeReal;

        should(get(i) == clx3);
        should(get(i, Diff2D(1,1)) == clx2);
        should(real(i) == 2.0);
        should(real(i, Diff2D(1,1)) == 0.0);
        should(imag(i) == -2.0);
        should(imag(i, Diff2D(1,1)) == -2.0);
        shouldEqual(mag(i) , sqrt(8.0));
        shouldEqual(mag(i, Diff2D(1,1)) , 2.0);
        shouldEqualTolerance(sin(phase(i)), -sqrt(2.0)/2.0, 1.0e-7);
        shouldEqual(sin(phase(i, Diff2D(1,1))), -1.0);

        writeReal.set(2.0, i);
        writeReal.set(2.0, i, Diff2D(1,1));
        should(get(i) == clx1);
        should(get(i, Diff2D(1,1)) == clx1);
    }

    void testForwardBackwardTrans()
    {
        const int w=256, h=256;

        FFTWComplexImage in(w, h);
        for (int y=0; y<in.height(); y++)
            for (int x=0; x<in.width(); x++)
            {
                in(x,y)= rand()/(double)RAND_MAX;
            }

        FFTWComplexImage out(w, h);
        
        fourierTransform(srcImageRange(in), destImage(out));
        fourierTransformInverse(srcImageRange(out), destImage(out));

        FindAverage<FFTWComplex<>::value_type> average;
        inspectImage(srcImageRange(out, FFTWImaginaryAccessor<>()), average);

        shouldEqualTolerance(average(), 0.0, 1e-14);

        combineTwoImages(srcImageRange(in), srcImage(out), destImage(out),
                         Compare1((double)w*h));

        average = FindAverage<FFTWMagnitudeAccessor<>::value_type>();
        inspectImage(srcImageRange(out, FFTWMagnitudeAccessor<>()), average);

        shouldEqualTolerance(average(), 0.0, 1e-14);

        for (int y=0; y<in.height(); y++)
            for (int x=0; x<in.width(); x++)
            {
                in(x,y)[1]= rand()/(double)RAND_MAX;
            }

        fourierTransform(srcImageRange(in), destImage(out));
        fourierTransformInverse(srcImageRange(out), destImage(out));

        combineTwoImages(srcImageRange(in), srcImage(out), destImage(out),
                         Compare1((double)w*h));

        FindAverage<FFTWComplex<> > caverage;
        inspectImage(srcImageRange(out), caverage);

        shouldEqualTolerance(caverage().magnitude(), 0.0, 1e-14);
    }

    void testRearrangeQuadrants()
    {
        double t4[] = { 0, 1, 2, 2,
                        1, 1, 2, 2,
                        3, 3, 4, 4,
                        3, 3, 4, 4};
        double t4res[] = { 4, 4, 3, 3,
                           4, 4, 3, 3,
                           2, 2, 0, 1,
                           2, 2, 1, 1};

        FFTWComplexImage in4(4,4), out4(4,4);
        copy(t4, t4+16, in4.begin());
        moveDCToCenter(srcImageRange(in4), destImage(out4));
        moveDCToUpperLeft(srcImageRange(out4), destImage(in4));
        for(int i=0; i<16; ++i)
        {
            should(out4.begin()[i] == t4res[i]);
            should(in4.begin()[i] == t4[i]);
        }

        double t5[] = { 0, 1, 1, 2, 2,
                        1, 1, 1, 2, 2,
                        1, 1, 1, 2, 2,
                        3, 3, 3, 4, 4,
                        3, 3, 3, 4, 4};
        double t5res[] = { 4, 4, 3, 3, 3,
                           4, 4, 3, 3, 3,
                           2, 2, 0, 1, 1,
                           2, 2, 1, 1, 1,
                           2, 2, 1, 1, 1};

        FFTWComplexImage in5(5,5), out5(5,5);
        copy(t5, t5+25, in5.begin());
        moveDCToCenter(srcImageRange(in5), destImage(out5));
        moveDCToUpperLeft(srcImageRange(out5), destImage(in5));
        for(int i=0; i<25; ++i)
        {
            should(out5.begin()[i] == t5res[i]);
            should(in5.begin()[i] == t5[i]);
        }
    }
};

struct MultiFFTTest
{
	typedef double R;
	typedef FFTWComplex<R> C;
	typedef MultiArray<2, R, FFTWAllocator<R> > DArray2;
	typedef MultiArray<2, C, FFTWAllocator<C> > CArray2;
	typedef MultiArrayShape<2>::type Shape2;
	typedef MultiArray<3, R, FFTWAllocator<R> > DArray3;
	typedef MultiArray<3, C, FFTWAllocator<C> > CArray3;
	typedef MultiArrayShape<3>::type Shape3;

	void testFFTShift()
	{
		Shape3 s(5,5,5);
		MultiArray<3, unsigned int> a(s), ref(s), iref(s);

		for(int z=0; z<5; ++z)
		{
			for(int y=0; y<5; ++y)
			{
				for(int x=0; x<5; ++x)
				{
					unsigned int v = ((z > 2) ? 4 : 0) + ((y > 2) ? 2 : 0) + ((x > 2) ? 1 : 0);
					a(x,y,z) = iref(x,y,z) = v;
					v = ((z < 2) ? 4 : 0) + ((y < 2) ? 2 : 0) + ((x < 2) ? 1 : 0);
					ref(x,y,z) = v;
				}
			}
		}

		moveDCToCenter(a);
		should(a == ref);
		moveDCToUpperLeft(a);
		should(a == iref);

		MultiArrayView<3, unsigned int> b = a.subarray(Shape3(1,1,1), s);
		moveDCToCenter(b);
		should(b == ref.subarray(Shape3(0,0,0), Shape3(4,4,4)));
		moveDCToUpperLeft(b);
		should(b == iref.subarray(Shape3(1,1,1), s));
		moveDCToCenter(b);
		moveDCToCenter(b);
		should(b == iref.subarray(Shape3(1,1,1), s));

		MultiArrayShape<1>::type s1_e(10), s1_5(5), s1_6(6);
		MultiArray<1, double> e1(s1_e), k1_5(s1_5), k1_6(s1_6);

		for(int k=0; k<k1_5.size(); ++k)
			k1_5(k) = k+1;
		detail::fftEmbedKernel(k1_5, e1);

		double k1_5ref[] = {3, 4, 5, 0, 0, 0, 0, 0, 1, 2};
		shouldEqualSequence(e1.data(), e1.data()+e1.size(), k1_5ref);

		for(int k=0; k<k1_6.size(); ++k)
			k1_6(k) = k+1;
		detail::fftEmbedKernel(k1_6, e1);

		double k1_6ref[] = {4, 5, 6, 0, 0, 0, 0, 1, 2, 3};
		shouldEqualSequence(e1.data(), e1.data()+e1.size(), k1_6ref);

		detail::fftEmbedArray(k1_5, e1);
		double a1_5ref[] = {3, 2, 1, 2, 3, 4, 5, 4, 3, 2 };
		shouldEqualSequence(e1.data(), e1.data()+e1.size(), a1_5ref);

		detail::fftEmbedArray(k1_6, e1);
		double a1_6ref[] = {3, 2, 1, 2, 3, 4, 5, 6, 5, 4 };
		shouldEqualSequence(e1.data(), e1.data()+e1.size(), a1_6ref);

		MultiArrayShape<2>::type s2_e(8, 8), s2_4(4, 4), s2_5(5, 5);
		MultiArray<2, double> e2(s2_e), k2_4(s2_4), k2_5(s2_5);

		for(int k=0; k<k2_4.size(); ++k)
			k2_4(k) = k+1;
		detail::fftEmbedKernel(k2_4, e2);

		double k2_4ref[] = {11,12, 0, 0, 0, 0, 9,10,
			                15,16, 0, 0, 0, 0,13,14,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 3, 4, 0, 0, 0, 0, 1, 2,
		                     7, 8, 0, 0, 0, 0, 5, 6};
		shouldEqualSequence(e2.data(), e2.data()+e2.size(), k2_4ref);

		for(int k=0; k<k2_5.size(); ++k)
			k2_5(k) = k+1;
		detail::fftEmbedKernel(k2_5, e2);

		double k2_5ref[] = {13,14,15, 0, 0, 0,11,12,
			                18,19,20, 0, 0, 0,16,17,
			                23,24,25, 0, 0, 0,21,22,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 0, 0, 0, 0, 0, 0, 0, 0,
			                 3, 4, 5, 0, 0, 0, 1, 2,
		                     8, 9,10, 0, 0, 0, 6, 7};
		shouldEqualSequence(e2.data(), e2.data()+e2.size(), k2_5ref);

		detail::fftEmbedArray(k2_4, e2);
		double a2_4ref[] = {11,10, 9,10,11,12,11,10,
			                 7, 6, 5, 6, 7, 8, 7, 6,
			                 3, 2, 1, 2, 3, 4, 3, 2,
			                 7, 6, 5, 6, 7, 8, 7, 6,
			                11,10, 9,10,11,12,11,10,
			                15,14,13,14,15,16,15,14,
			                11,10, 9,10,11,12,11,10,
		                     7, 6, 5, 6, 7, 8, 7, 6};
		shouldEqualSequence(e2.data(), e2.data()+e2.size(), a2_4ref);

		detail::fftEmbedArray(k2_5, e2);
		double a2_5ref[] = { 7, 6, 7, 8, 9,10, 9, 8,
			                 2, 1, 2, 3, 4, 5, 4, 3,
			                 7, 6, 7, 8, 9,10, 9, 8,
			                12,11,12,13,14,15,14,13,
			                17,16,17,18,19,20,19,18,
			                22,21,22,23,24,25,24,23,
			                17,16,17,18,19,20,19,18,
		                    12,11,12,13,14,15,14,13};
		shouldEqualSequence(e2.data(), e2.data()+e2.size(), a2_5ref);
	}

	void testFFT2D()
	{
		Shape2 s(256, 256);

		FFTWComplexImage in(s[0], s[1]);
        for (int y=0; y<in.height(); y++)
            for (int x=0; x<in.width(); x++)
            {
                in(x,y)= rand()/(double)RAND_MAX;
            }

        FFTWComplexImage out(in.size());
		CArray2 aout(s);

        fourierTransform(srcImageRange(in), destImage(out));
		fourierTransform(MultiArrayView<2, C>(s, const_cast<C*>(in.data())), 
			             aout);

		shouldEqualSequence(aout.data(), aout.data()+aout.size(), out.data());

		DArray2 rin(s);
		copyImage(srcImageRange(in, FFTWRealAccessor<>()), destImage(rin));

		fourierTransform(rin, aout);
		shouldEqualSequence(aout.data(), aout.data()+aout.size(), out.data());
	}

	void testFFT3D()
	{
		Shape3 s(32, 24, 16);
		CArray3 r(s), ir(s), irr(s);

		fourierTransform(MultiArrayView<3, double>(s, f3data), r);

		DArray3 re(s);
		copyMultiArray(srcMultiArrayRange(r, FFTWRealAccessor<>()), destMultiArray(re));
		shouldEqualSequenceTolerance(re.data(), re.data()+re.size(), f3ref, 1e-10);

		FindMinMax<double> minmax;
		inspectMultiArray(srcMultiArrayRange(r, FFTWImaginaryAccessor<>()), minmax);
		shouldEqualTolerance(minmax.min, 0.0, 1e-10);
		shouldEqualTolerance(minmax.max, 0.0, 1e-10);

		fourierTransformInverse(r, ir);

		copyMultiArray(srcMultiArrayRange(ir, FFTWRealAccessor<>()), destMultiArray(re));
		shouldEqualSequenceTolerance(re.data(), re.data()+re.size(), f3data, 1e-10);

		inspectMultiArray(srcMultiArrayRange(ir, FFTWImaginaryAccessor<>()), minmax);
		shouldEqualTolerance(minmax.min, 0.0, 1e-10);
		shouldEqualTolerance(minmax.max, 0.0, 1e-10);
	}

	void testPadding()
	{
		shouldEqual(0, detail::fftwPaddingSize(0));
		shouldEqual(1, detail::fftwPaddingSize(1));
		shouldEqual(3, detail::fftwPaddingSize(3));
		shouldEqual(256, detail::fftwPaddingSize(255));
		shouldEqual(256, detail::fftwPaddingSize(256));
		shouldEqual(260, detail::fftwPaddingSize(257));
		shouldEqual(0, detail::fftwEvenPaddingSize(0));
		shouldEqual(2, detail::fftwEvenPaddingSize(1));
		shouldEqual(4, detail::fftwEvenPaddingSize(3));
		shouldEqual(256, detail::fftwEvenPaddingSize(255));
		shouldEqual(256, detail::fftwEvenPaddingSize(256));
		shouldEqual(260, detail::fftwEvenPaddingSize(257));

		Shape3 s(113, 256, 257);
		shouldEqual(Shape3(117, 256, 260), fftwBestPaddedShape(s));
		shouldEqual(Shape3(120, 256, 260), fftwBestPaddedShapeR2C(s));
	}

	void testConvolveFFT()
	{
		typedef MultiArrayView<2, double> MV;
		ImageImportInfo info("ghouse.gif");
		Shape2 s(info.width(), info.height());
		DArray2 in(s), out(s), ref(s), out2(s), ref2(s);
		importImage(info, destImage(in));

		double scale = 2.0;
		gaussianSmoothing(srcImageRange(in), destImage(ref), scale);
		gaussianSmoothing(srcImageRange(in), destImage(ref2), 2.0*scale);

		Kernel2D<double> gauss;
		gauss.initGaussian(scale);
		MV kernel(Shape2(gauss.width(), gauss.height()), &gauss[gauss.upperLeft()]);
		convolveFFT(in, kernel, out);

		shouldEqualSequenceTolerance(out.data(), out.data()+out.size(),
			                         ref.data(), 1e-14);
		
		Kernel2D<double> gauss2;
		gauss2.initGaussian(2.0*scale);
		MV kernel2(Shape2(gauss2.width(), gauss2.height()), &gauss2[gauss2.upperLeft()]);

		MV kernels[] = { kernel, kernel2 };
		MV outs[] = { out, out2 };
		convolveFFTMany(in, kernels, kernels+2, outs);

		shouldEqualSequenceTolerance(out.data(), out.data()+out.size(),
			                         ref.data(), 1e-14);
		shouldEqualSequenceTolerance(out2.data(), out2.data()+out2.size(),
			                         ref2.data(), 1e-14);
	}

	void testConvolveFourierKernel()
	{
		ImageImportInfo info("ghouse.gif");
		Shape2 s(info.width(), info.height());
		DArray2 in(s), out(s), ref(s), tmp(s);
		importImage(info, destImage(in));

		double scale = 2.0;
		gaussianSmoothing(srcImageRange(in), destImage(ref), scale);

		Shape2 paddedShape = fftwBestPaddedShapeR2C(s + Shape2(16)),
			   kernelShape(paddedShape),
			   complexShape(paddedShape);
		kernelShape[0] = kernelShape[0] / 2 + 1;
		complexShape[0] += 1;
		Shape2 center = div(complexShape, Shape2::value_type(2));

		CArray2 kernel(complexShape);
		for(int y=0; y<complexShape[1]; ++y)
		{
			for(int x=0; x<complexShape[0]; ++x)
			{
				double xx = 2.0 * M_PI * (x - center[0]) / complexShape[0];
				double yy = 2.0 * M_PI * (y - center[1]) / complexShape[1];
				double r2 = sq(xx) + sq(yy);
				kernel(x,y) = std::exp(-0.5 * sq(scale) * r2);
			}
		}
		moveDCToUpperLeft(kernel);

		FFTWConvolvePlan<2, double> plan;
		plan.initFourierKernel(in, kernel.subarray(Shape2(), kernelShape), out);
		plan.executeFourierKernel(in, kernel.subarray(Shape2(), kernelShape), out);

		shouldEqualSequenceTolerance(out.data(), out.data()+out.size(),
			                         ref.data(), 1e-2);

		Kernel1D<double> gauss, grad;
		gauss.initGaussian(scale);
		grad.initGaussianDerivative(scale, 1);

		separableConvolveX(srcImageRange(in), destImage(tmp), kernel1d(grad));
		separableConvolveY(srcImageRange(tmp), destImage(ref), kernel1d(gauss));

		for(int y=0; y<complexShape[1]; ++y)
		{
			for(int x=0; x<complexShape[0]; ++x)
			{
				double xx = 2.0 * M_PI * (x - center[0]) / complexShape[0];
				double yy = 2.0 * M_PI * (y - center[1]) / complexShape[1];
				double r2 = sq(xx) + sq(yy);
				kernel(x,y) = C(0, xx*std::exp(-0.5 * sq(scale) * r2));
			}
		}
		moveDCToUpperLeft(kernel);

		plan.executeFourierKernel(in, kernel.subarray(Shape2(), kernelShape), out);

		ref -= out;

		FindMinMax<double> minmax;
		FindAverage<double> average;
		inspectImage(srcImageRange(ref), minmax);
		inspectImage(srcImageRange(ref), average);
		
		should(std::max(minmax.max, -minmax.min) < 0.2);
		should(average.average() < 0.001);
	}
};

struct FFTWTestSuite
: public test_suite
{

    FFTWTestSuite()
    : test_suite("FFTWTest")
    {

        add( testCase(&FFTWComplexTest::testConstruction));
        add( testCase(&FFTWComplexTest::testComparison));
        add( testCase(&FFTWComplexTest::testArithmetic));
        add( testCase(&FFTWComplexTest::testAccessor));
        add( testCase(&FFTWComplexTest::testForwardBackwardTrans));
        add( testCase(&FFTWComplexTest::testRearrangeQuadrants));

        add( testCase(&MultiFFTTest::testFFTShift));
        add( testCase(&MultiFFTTest::testFFT2D));
        add( testCase(&MultiFFTTest::testFFT3D));
        add( testCase(&MultiFFTTest::testPadding));
        add( testCase(&MultiFFTTest::testConvolveFFT));
        add( testCase(&MultiFFTTest::testConvolveFourierKernel));
    }
};

int initializing;

struct CompareFunctor
{
	double sumDifference_;

	CompareFunctor(): sumDifference_(0) {}

	void operator()(const fftw_real &a, const fftw_real &b)
		{ sumDifference_+= abs(a-b); }

    double operator()()
		{ return sumDifference_; }
};

struct GaborTests
{
	ImageImportInfo info;
	int w, h;
	BasicImage<fftw_real> image;

	GaborTests()
		: info("ghouse.gif"),
		  w(info.width()), h(info.height()),
		  image(w, h)
	{
		importImage(info, destImage(image));
	}

	template<class Iterator, class Accessor>
	void checkImage(triple<Iterator, Iterator, Accessor> src, const char *filename)
	{
		if (initializing)
		{
			exportImage(src, ImageExportInfo(filename));
			cout << "wrote " << filename << endl;
		}
		else
		{
			ImageImportInfo info(filename);
			FImage image(info.width(), info.height());
			importImage(info, destImage(image));

			CompareFunctor cmp;

			inspectTwoImages(src, srcImage(image), cmp);
			cout << "difference to " << filename << ": " << cmp() << endl;
			shouldEqualTolerance(cmp(), 0.0, 1e-4);
		}
	}

	void testImages()
	{
		BasicImage<fftw_real> filter(w, h);
		int directionCount= 8;
		int dir= 1, scale= 1;
		double angle = dir * M_PI / directionCount;
		double centerFrequency = 3.0/8.0 / VIGRA_CSTD::pow(2.0,scale);
		createGaborFilter(destImageRange(filter),
						  angle, centerFrequency,
						  angularGaborSigma(directionCount, centerFrequency),
						  radialGaborSigma(centerFrequency));

		checkImage(srcImageRange(filter), "filter.xv");

		cout << "Applying filter...\n";
		FFTWComplexImage result(w, h);
		applyFourierFilter(srcImageRange(image), srcImage(filter), destImage(result));
		checkImage(srcImageRange(result, FFTWMagnitudeAccessor<>()), "gaborresult.xv");

		BasicImage<fftw_real> realPart(w, h);
		applyFourierFilter(srcImageRange(image), srcImage(filter), destImage(realPart));

		CompareFunctor cmp;
		inspectTwoImages(srcImageRange(result, FFTWRealAccessor<>()),
						 srcImage(realPart), cmp);
		cout << "difference between real parts: " << cmp() << endl;
		shouldEqualTolerance(cmp(), 0.0, 1e-4);

		cout << "testing vector results..\n";
		FVector2Image vectorResult(w,h);
		applyFourierFilter(srcImageRange(image), srcImage(filter),
						   destImage(vectorResult));
		CompareFunctor iCmp;
		inspectTwoImages(srcImageRange(result, FFTWImaginaryAccessor<>()),
						 srcImage(vectorResult,
								  VectorComponentAccessor<FVector2Image::PixelType>(1)),
						 iCmp);
		cout << "difference between imaginary parts: " << iCmp() << endl;
		shouldEqualTolerance(cmp(), 0.0, 1e-4);

		cout << "applying on ROI...\n";
		BasicImage<fftw_real> bigImage(w+20, h+20);
		BasicImage<fftw_real>::Iterator bigUL= bigImage.upperLeft() + Diff2D(5, 5);
		copyImage(srcImageRange(image), destIter(bigUL));
		applyFourierFilter(srcIterRange(bigUL, bigUL + Diff2D(w, h)),
						   srcImage(filter), destImage(result));
		checkImage(srcImageRange(result, FFTWMagnitudeAccessor<>()), "gaborresult.xv");
#if 0
		cout << "Creating plans with measurement...\n";
		fftwnd_plan forwardPlan=
			fftw2d_create_plan(h, w, FFTW_FORWARD, FFTW_MEASURE );
		fftwnd_plan backwardPlan=
			fftw2d_create_plan(h, w, FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE);

		cout << "Applying again...\n";
		applyFourierFilter(srcImageRange(image), srcImage(filter), destImage(result),
						   forwardPlan, backwardPlan);

		checkImage(srcImageRange(result, FFTWMagnitudeAccessor()), "gaborresult.xv");
#endif
	}

	void testFamily()
	{
		cout << "testing 8x2 GaborFilterFamily...\n";
		GaborFilterFamily<FImage> filters(w, h, 8, 2);
		ImageArray<FFTWComplexImage> results((unsigned int)filters.size(), filters.imageSize());
		applyFourierFilterFamily(srcImageRange(image), filters, results);
		checkImage(srcImageRange(results[filters.filterIndex(1,1)], FFTWMagnitudeAccessor<>()),
				   "gaborresult.xv");

		ImageArray<FImage> realParts((unsigned int)filters.size(), filters.imageSize());
		applyFourierFilterFamily(srcImageRange(image), filters, realParts);

		CompareFunctor cmp;
		inspectTwoImages(srcImageRange(results[3], FFTWRealAccessor<>()),
						 srcImage(realParts[3]), cmp);
		cout << "difference between real parts: " << cmp() << endl;
		shouldEqualTolerance(cmp(), 0.0, 1e-4);
	}
};

struct GaborTestSuite
: public test_suite
{
    GaborTestSuite()
    : test_suite("FFTWWithGaborTest")
    {
        add(testCase(&GaborTests::testImages));
        add(testCase(&GaborTests::testFamily));
    }
};

int main(int argc, char **argv)
{
    initializing= argc>1 && std::string(argv[1]) == "init"; // global variable,
    // initializing means: don't compare with image files but create them

    FFTWTestSuite fftwTest;

    int failed = fftwTest.run(testsToBeExecuted(argc, argv));

    std::cout << fftwTest.report() << std::endl;

    GaborTestSuite gaborTest;

    failed = failed || gaborTest.run(testsToBeExecuted(argc, argv));

    std::cout << gaborTest.report() << std::endl;

    return (failed != 0);
}
