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

// -*- c++ -*-
// $Id$

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <functional>
#include <cmath>

#include "vigra/unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/convolution.hxx" 
#include "vigra/navigator.hxx"
#include "vigra/random.hxx"

#include "vigra/impex.hxx"
#include "vigra/imageinfo.hxx"
#include "vigra/basicimage.hxx"
#include "vigra/diff2d.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/multi_resize.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/bordertreatment.hxx"
#include "vigra/recursive_multi_convolution.hxx"

using namespace vigra;
using namespace vigra::functor;

struct MultiArraySeparableConvolutionTest
{

    typedef float PixelType;

    typedef TinyVector<PixelType,3> VectorPixelType;
    typedef MultiArray<3, VectorPixelType> Image3x3;

    typedef MultiArray<3,PixelType> Image3D;
    typedef Image3D::difference_type Size3;

    typedef BasicImage<PixelType> Image2D;
    typedef Image2D::difference_type Size2;


    // - - - - - - - - - - - - - - - - - - - - - - - -

    template<class Image>
    void makeRandom( Image &image )
    {
        typedef typename Image::value_type T;

        int size = image.size();
        for( int k = 0; k < size; ++k ) 
        {
            typedef typename NumericTraits<typename Image::value_type>::isIntegral isIntegral;
            if(isIntegral::value)
                image[k] = (T)randomMT19937().uniformInt(256);
            else
                image[k] = (T)randomMT19937().uniform();
        }
    }

    void makeWedge( Image3D &image )
    {
        const Size3 shape = image.shape();
        const int width = shape[0];
        const int height = shape[1];
        const int depth = shape[2];
        for( int z = 0; z < depth; ++z ) 
        {
            for( int y = 0; y < height; ++y ) 
            {
                for( int x = 0; x < width; ++x ) 
                {
                    const Image3D::value_type val = Image3D::value_type(x + y + z);
                    image( x, y, z ) = val;
                }
            }
        }
    }


    void makeBox( Image3D &image )
    {
        const int b = 8;
        const Size3 shape = image.shape();
        const int width = shape[0];
        const int height = shape[1];
        const int depth = shape[2];
        for( int z = 0; z < depth; ++z ) 
        {
            for( int y = 0; y < height; ++y ) 
            {
                for( int x = 0; x < width; ++x ) 
                {
                
                    Image3D::value_type val = 80;

                    if( (x > b) && x < (width-b) &&
                            (y > b) && y < (height-b) &&
                            (z > b)    && z < (depth-b) ) 
                    {
                        val = 220;
                    }
                    image( x, y, z ) = val;
                }
            }
        }
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - -

    void test_1DValidity( const Image3D &src, float ksize )
    {
        Image3D d1( src.shape() );
        Image3D dn( src.shape() );
        Image3D dest3( src.shape() );

        for( int d = 0; d < 3; ++d ) 
        {
            std::vector<vigra::Kernel1D<float> > kernels( 3 );
            kernels[d].initGaussianDerivative( ksize, 1 );

            separableConvolveMultiArray( srcMultiArrayRange(src),
                                         destMultiArray(dn),
                                         kernels.begin() );

            convolveMultiArrayOneDimension( srcMultiArrayRange(src),
                                            destMultiArray(d1),
                                            d,
                                            kernels[d] );
 
            shouldEqualSequence(    dn.begin(), dn.end(), d1.begin() );
        }

        for( int d = 0; d < 3; ++d ) 
        {
            std::vector<vigra::Kernel1D<float> > kernels( 3 );
            kernels[d].initGaussianDerivative( ksize, 1 );

            separableConvolveMultiArray(src, dn, kernels.begin());

            convolveMultiArrayOneDimension(src, d1, d, kernels[d] );
 
            shouldEqualSequence(    dn.begin(), dn.end(), d1.begin() );
        }
    }

    // - - - - - - - - - - - - - - - - - - - - - - - -

    void test_1DValidityB( const Image3D &src, float ksize )
    {
        Image3D dst1( src.shape() );
        Image3D dst2( src.shape() );

        const int depth = src.shape()[2];

        vigra::Kernel1D<float> kernx;
        vigra::Kernel1D<float> kerny;
        kernx.initGaussian( ksize );
        kerny.initGaussianDerivative( ksize, 1 );

        int z;

        // x convolution
        convolveMultiArrayOneDimension( srcMultiArrayRange(src),
                                        destMultiArray(dst1),
                                        0, kernx );

        for( z = 0; z < depth; ++z ) 
        {
            BasicImageView<Image3D::value_type> sslice =
                makeBasicImageView( src.bindOuter(z) );
            BasicImageView<Image3D::value_type> dslice =
                makeBasicImageView( dst2.bindOuter(z) );
            
            vigra::separableConvolveX( srcImageRange(sslice), destImage(dslice),
                                                                 vigra::kernel1d(kernx) );
        }

        shouldEqualSequence(    dst1.begin(), dst1.end(), dst2.begin() );


        // y convolution
        convolveMultiArrayOneDimension( srcMultiArrayRange(src),
                                        destMultiArray(dst1),
                                        1, kerny );

        for( z = 0; z < depth; ++z ) 
        {
            BasicImageView<Image3D::value_type> sslice =
                makeBasicImageView( src.bindOuter(z) );
            BasicImageView<Image3D::value_type> dslice =
                makeBasicImageView( dst2.bindOuter(z) );
            
            vigra::separableConvolveY( srcImageRange(sslice), destImage(dslice),
                                                                 vigra::kernel1d(kerny) );
        }

        shouldEqualSequence(    dst1.begin(), dst1.end(), dst2.begin() );
    }

    void test_2DValidity( const Image3D &src, float ksize )
    {
        Image3D d2( src.shape() );
        Image3D dn( src.shape() );

        int depth = src.shape()[2];

        std::vector<vigra::Kernel1D<float> > kernels( 3 );
        kernels[0].initGaussian( ksize );
        kernels[1].initGaussianDerivative( ksize, 1 );

        for( int z = 0; z < depth; ++z )
        {
            BasicImageView<Image3D::value_type> sslice =
                makeBasicImageView( src.bindOuter(z) );
            BasicImageView<Image3D::value_type> dslice =
                makeBasicImageView( d2.bindOuter(z) );
            
            vigra::convolveImage( srcImageRange(sslice), destImage(dslice),
                                    kernels[0], kernels[1] );
        }

        separableConvolveMultiArray( srcMultiArrayRange(src),
                                     destMultiArray(dn),
                                     kernels.begin() );

        shouldEqualSequence( dn.begin(), dn.end(), d2.begin() );
    }


    // - - - - - - - - - - - - - - - - - - - - - - - -

    void test_inplacenessN( const Image3D &src, float ksize )
    {
        Image3D da( src.shape() );
        Image3D db( src.shape() );

        std::vector<vigra::Kernel1D<float> > kernels( 3 );
        kernels[0].initGaussian( ksize );
        kernels[1].initGaussianDerivative( ksize, 1 );

        vigra::separableConvolveMultiArray( srcMultiArrayRange(src),
                                            destMultiArray(da),
                                            kernels.begin() );

        copyMultiArray(srcMultiArrayRange(src), destMultiArray(db));

        vigra::separableConvolveMultiArray( srcMultiArrayRange(db),
                                            destMultiArray(db),
                                            kernels.begin() );

        shouldEqualSequenceTolerance( da.begin(), da.end(),
                                        db.begin(),
                                        1e-5 );
    }


    // - - - - - - - - - - - - - - - - - - - - - - - -

    void testSmoothing()
    {
        makeRandom(srcImage);

        ArrayVector<Kernel1D<double> > kernels(3);
        kernels[0].initGaussian(1.0);
        kernels[1].initAveraging(1);
        kernels[2].initGaussian(2.0);

        //kernels[1].setBorderTreatment(BORDER_TREATMENT_REFLECT);

        Image3D res(shape);
        separableConvolveMultiArray(srcMultiArrayRange(srcImage), destMultiArray(res), 
                                    kernels.begin());

        typedef Shape3 S;
        int w = shape[0], h = shape[1], d = shape[2];
        S start[] = { S(0, 0, 25), S(30, 0, 0), S(0, 30, 0), 
                      S(2, h-1, 20), S(15, 1, d-2), S(2, 14, 1), 
                      S(0,0,0), S(0, 0, 34), S(0, 12, 0) };
        S stop[]  = { S(w, h, 26), S(31, h, d), S(w, 31, d), 
                      S(w-2, h, 30), S(28, h-1, d), S(w-2, 44, d-1), 
                      S(w,h,d), S(w, 1, 39), S(w-1, 26, 1) };
        for(int k=0; k<9; ++k)
        {
            Image3D subarray(stop[k]-start[k]);
            separableConvolveMultiArray(srcImage, subarray, 
                                        kernels.begin(), start[k], stop[k]);

            shouldEqualSequenceTolerance(subarray.begin(), subarray.end(), 
                                         res.subarray(start[k], stop[k]).begin(), 1e-6);
        }
    }

    void test_inplaceness1( const Image3D &src, float ksize, bool useDerivative )
    {
        Image3D da( src.shape() );
        Image3D db( src.shape() );

        Kernel1D<float> kernel;
        if( ! useDerivative )
            kernel.initGaussian( ksize );
        else
            kernel.initGaussianDerivative( ksize, 1 );


        for( int i = 0; i < 3; ++i ) 
        {
            const int d = 2-i;

            convolveMultiArrayOneDimension(src, da, d, kernel );

            copyMultiArray(src, db);

            convolveMultiArrayOneDimension(db, db, d, kernel );

            shouldEqualSequence( da.begin(), da.end(), db.begin() );
        }
    }


    // - - - - - - - - - - - - - - - - - - - - - - - -

    void test_gradient1( const Image3D &base, bool useGaussian )
    {
        const double sigma = kernelSize/2;
        const int b = useGaussian ? int( 0.5 + 3*sigma ) : 1;
        Image3D src( base.shape() );
        Image3x3 grad( src.shape() );
        makeWedge( src );

        if( ! useGaussian )
            symmetricGradientMultiArray(src, grad);
        else
            gaussianGradientMultiArray(src, grad, sigma, ConvolutionOptions<3>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));

        Image3x3::value_type v;
        v[0] = 1; v[1] = 1; v[2] = 1;
        const float v2 = dot(v,v);

        const Size3 shape = src.shape();
        const int width = shape[0];
        const int height = shape[1];
        const int depth = shape[2];
        for( int z = b; z < depth-b; ++z ) 
        {
            for( int y = b; y < height-b; ++y ) 
            {
                for( int x = b; x < width-b; ++x ) 
                {
                    shouldEqualTolerance( dot(grad(x,y,z), v), v2, 1e-5 );
                }
            }
        }

    }


    void test_gradient_magnitude()
    {
        using namespace functor;
        
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape), mgrad(shape), rmgrad(shape);
        MultiArray<2, TinyVector<double, 2> > grad(shape);
        
        makeRandom(src);
        
        gaussianGradientMagnitude(src, rmgrad, 1.0);
        gaussianGradientMultiArray(src, grad, 1.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        transformMultiArray(grad, mgrad, norm(Arg1()));

        shouldEqualSequence(mgrad.data(), mgrad.data()+size, rmgrad.data());

        rmgrad.init(0);
        gaussianGradientMagnitude(src, rmgrad, 1.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        shouldEqualSequence(mgrad.data(), mgrad.data()+size, rmgrad.data());

        MultiArray<2, TinyVector<double, 3> > rgb(shape);
        MultiArrayView<3, Multiband<double> > expanded(rgb.expandElements(2));
       
        makeRandom(expanded);
        
        mgrad.init(0);
        gaussianGradientMagnitude(srcImageRange(rgb), destImage(mgrad), 1.0);
        rmgrad.init(0);
        gaussianGradientMagnitude(rgb, rmgrad, 1.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        shouldEqualSequenceTolerance(mgrad.data(), mgrad.data()+size, rmgrad.data(), 1e-14);

        rmgrad.init(0);
        gaussianGradientMagnitude<2>(expanded, rmgrad, 1.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        shouldEqualSequenceTolerance(mgrad.data(), mgrad.data()+size, rmgrad.data(), 1e-14);

        MultiArray<3, Multiband<double> > spectral(Shape3(shape[0], shape[1], 10));
        MultiArrayView<3, Multiband<double> > spectral_expanded(spectral);
        gaussianGradientMagnitude<2>(spectral_expanded, rmgrad, 1.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
    }

    void test_laplacian()
    {
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape), laplacian(shape);
        MultiArray<2, double> rlaplacian(shape[0], shape[1]);
        
        makeRandom(src);
        
        laplacianOfGaussian(srcImageRange(src), destImage(rlaplacian), 2.0);
        laplacianOfGaussianMultiArray(srcMultiArrayRange(src), destMultiArray(laplacian), 2.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));

        shouldEqualSequenceTolerance(laplacian.data(), laplacian.data()+size, rlaplacian.data(), 1e-12);

        laplacian = 0;
        laplacianOfGaussianMultiArray(src, laplacian, 2.0,ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        shouldEqualSequenceTolerance(laplacian.data(), laplacian.data()+size, rlaplacian.data(), 1e-12);
    }

    void test_divergence()
    {
        // test divergence of gradient - theoretically, it equals the Laplacian, but in practice this requires
        //                               large kernel windows, a high tolerance, and doesn't hold near the border
        MultiArray<2, double> src;
        importImage("oi_single.gif", src);

        MultiArray<2, double> laplacian(src.shape()), divergence(src.shape());
        MultiArray<2, TinyVector<double, 2> > grad(src.shape());

        laplacianOfGaussianMultiArray(src, laplacian, 2.0, ConvolutionOptions<2>().filterWindowSize(5).setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));

        gaussianGradientMultiArray(src, grad, sqrt(2.0), ConvolutionOptions<2>().filterWindowSize(5).setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        gaussianDivergenceMultiArray(grad, divergence, sqrt(2.0), ConvolutionOptions<2>().filterWindowSize(5).setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));

        divergence -= laplacian;
        MultiArrayView <2, double> center(divergence.subarray(Shape2(10,10), Shape2(-10,-10)));

        using namespace multi_math;
        should(all(abs(center) < 0.001));
    }

    void test_hessian()
    {
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape);
        MultiArray<2, TinyVector<double, 3> > hessian(shape), hessian1(shape);
        BasicImage<TinyVector<double, 3> > rhessian(shape[0], shape[1]);
        
        makeRandom(src);
        
        typedef VectorComponentAccessor<TinyVector<double, 3> > BandAccessor;
        hessianMatrixOfGaussian(srcImageRange(src), 
                                destImage(rhessian, BandAccessor(0)), 
                                destImage(rhessian, BandAccessor(1)), 
                                destImage(rhessian, BandAccessor(2)), 
                                2.0);
        hessianOfGaussianMultiArray(srcMultiArrayRange(src), destMultiArray(hessian), 2.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));

        TinyVector<double, 3> epsilon(1e-12, 1e-12, 1e-12);
        shouldEqualSequenceTolerance(hessian.data(), hessian.data()+size, rhessian.data(), epsilon);

        hessianOfGaussianMultiArray(src, hessian1, 2.0 );
        shouldEqualSequenceTolerance(hessian1.data(), hessian1.data()+size, rhessian.data(), epsilon);
    }

    void test_structureTensor()
    {
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape);
        MultiArray<2, TinyVector<double, 3> > st(shape), st1(shape);
        BasicImage<TinyVector<double, 3> > rst(shape[0], shape[1]);
        
        makeRandom(src);
        
        typedef VectorComponentAccessor<TinyVector<double, 3> > BandAccessor;
        structureTensor(srcImageRange(src), 
                        destImage(rst, BandAccessor(0)), 
                        destImage(rst, BandAccessor(1)), 
                        destImage(rst, BandAccessor(2)), 
                        1.5, 3.0);
        structureTensorMultiArray(srcMultiArrayRange(src), destMultiArray(st), 1.5, 3.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));

        TinyVector<double, 3> epsilon(1e-12, 1e-12, 1e-12);
        shouldEqualSequenceTolerance(st.data(), st.data()+size, rst.data(), epsilon);

        structureTensorMultiArray(src, st1, 1.5, 3.0, ConvolutionOptions<2>().setKernelApproximation(MULTI_CONVOLUTION_KERNEL_FIR));
        shouldEqualSequenceTolerance(st1.data(), st1.data()+size, rst.data(), epsilon);
    }

    //--------------------------------------------

    const Size3 shape;
    Image3D srcImage;
    const float kernelSize;

    MultiArraySeparableConvolutionTest()
        : shape( 60, 70, 50 ),
            srcImage( shape ),
            kernelSize( 1.8f )
    {
        makeBox( srcImage );
    }


    void test_Valid1() 
    {
        test_1DValidity( srcImage, kernelSize );
    }

    void test_Valid3() 
    {
     test_1DValidityB( srcImage, kernelSize );
    }

    void test_Valid2() 
    {
        test_2DValidity( srcImage, kernelSize );
    }

    void test_InplaceN() 
    {
        test_inplacenessN( srcImage, kernelSize );
    }

    void test_Inplace1() 
    {
        test_inplaceness1( srcImage, kernelSize, false );
        test_inplaceness1( srcImage, kernelSize, true );
    }

    void test_gradient1() 
    {
        test_gradient1( srcImage, false );
        test_gradient1( srcImage, true );
    }
};                //-- struct MultiArraySeparableConvolutionTest

//--------------------------------------------------------

struct MultiArraySeparableConvolutionTestSuite
: public vigra::test_suite
{
    MultiArraySeparableConvolutionTestSuite()
        : vigra::test_suite("MultiArraySeparableConvolutionTestSuite")
        {
                add( testCase( &MultiArraySeparableConvolutionTest::test_Valid1 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_Valid2 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_Valid3 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_InplaceN ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_Inplace1 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::testSmoothing ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_gradient1 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_laplacian ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_divergence ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_hessian ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_structureTensor ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_gradient_magnitude ) );
    }
}; // struct MultiArraySeparableConvolutionTestSuite

//--------------------------------------------------------

typedef double pixel_type;

typedef vigra::MultiArrayShape<2>::type                         shape_2d;
typedef vigra::MultiArray<2, pixel_type>                        array_2d;
typedef vigra::MultiArray<2, vigra::TinyVector<pixel_type, 2> > grad_2d;
typedef vigra::MultiArray<2, vigra::TinyVector<pixel_type, 3> > symm_2d;
typedef vigra::BasicImageView<pixel_type>                       image_2d;
typedef vigra::ConvolutionOptions<2>                            options_2d;

struct delta_dist
{
    double scale;
    double offset;
    delta_dist(double s, double of = 0) : scale(std::abs(s)), offset(of) {}
    template<class V>
    double operator()(V a, V b) const
    {
        double delta = std::abs(double(a) - double(b));
        if (offset == 0)
            return delta;
        else
            return offset - delta;
    }
};

struct rel_dist
{
    double scale;
    double offset;
    rel_dist(double s, double of = 0) : scale(std::abs(s)), offset(of) {}
    template<class V>
    double operator()(V a, V b) const
    {
        double delta = std::abs((double(a) - double(b)) / double(a));
        if (offset == 0)
            return delta;
        else
            return offset - delta;
    }
};

struct logger
{
    const bool do_log;
    std::ostringstream os;

    logger(bool d = true) : do_log(d) {}

    operator bool()
    {
        return do_log;
    }

    void next()
    {
        if(os.str().size())
            os << ", ";
    }

    void operator()(const std::string & msg)
    {
        next();
        os << msg;
    }
    template<class X>
    void operator()(const std::string & name, const X & value)
    {
        next();
        os << name << " = " << value;
    }

    std::string str() {
        return os.str();
    }
};

template<class IM>
double min_max_delta(const IM & image)
{
    vigra::FindMinMax<double> minmax;
    vigra::inspectMultiArray(image, minmax);
    return minmax.max - minmax.min;
}

template<class IM>
double max(const IM & image)
{
    vigra::FindMinMax<double> minmax;
    vigra::inspectMultiArray(image, minmax);
    return minmax.max;
}

template<class IM>
double min(const IM & image)
{
    vigra::FindMinMax<double> minmax;
    vigra::inspectMultiArray(image, minmax);
    return minmax.min;
}

template<class IM>
void write_out(const IM & image, const char *const name, bool do_it = true)
{
    if (do_it)
    {
        vigra::exportImage(srcImageRange(image),
                        vigra::ImageExportInfo(name));
    }
}

template<class F>
double compare(F func, const char *const func_name,
        const array_2d & image_a, const array_2d & image_b,
        double sigma_pix = 0, const char* name = 0, bool do_it = true) {
    array_2d delta_image(image_a);
    vigra::combineTwoMultiArrays(srcMultiArrayRange(image_a),
                                 srcMultiArray(image_b),
                                 destMultiArray(delta_image),
                                 func);

    // crop off border artifacts
    image_2d image_view_delta = vigra::makeBasicImageView(delta_image);
    // accounting for x step size by using "pixel" counts ...
    int crop = 1 + int(3 * sigma_pix + 0.5);
    int wcrop = image_view_delta.width() - 2 * crop;
    shouldMsg(wcrop > 1, ("no width left after std. dev. crop at"
            " borders, width = " + asString(wcrop)).c_str());

    array_2d delta_cropped(shape_2d(wcrop, image_view_delta.height()));
    copyImage(srcIterRange(
        image_view_delta.upperLeft() + vigra::Diff2D(crop, 0),
        image_view_delta.lowerRight() - vigra::Diff2D(crop, 0)),
        destImage(delta_cropped));

    if (!name)
        name = func_name;
    write_out(delta_cropped, name, do_it);

    return max(delta_cropped);
}

void resize_0(const array_2d & a, array_2d & b)
{
    vigra::resizeImageNoInterpolation(srcImageRange(a), destImageRange(b));
}
template<unsigned order>
void resize_n(const array_2d & a, array_2d & b)
{
    typedef typename vigra::BSpline<order, double> my_spline;

    if(order <=3) // arbitrarily use the different APIs for small and large orders
        vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(a),
                                                   destMultiArrayRange(b),
                                                   my_spline());
    else
        vigra::resizeMultiArraySplineInterpolation(a, b, my_spline());
}

struct t_func
{
    static const double zeros[2]; // = { 0, 0 };
    virtual double operator()(array_2d & image, const options_2d & opt) const = 0;
    virtual double operator()(array_2d & image,
                              const double sigma_d[2],
                              const double step_size[2],
                              const MultiConvolutionKernel kernel_approx) const
    {
        // test various ways to create option setter functions.
        options_2d opt;
        opt.resolutionStdDev(sigma_d).stepSize(step_size);
        std::vector<double> vstep;
        vstep.push_back(step_size[0]);
        vstep.push_back(step_size[1]);
        opt.resolutionStdDev(sigma_d).stepSize((const double *const)step_size);
        opt.resolutionStdDev(sigma_d).stepSize(vstep);
        opt.setKernelApproximation(kernel_approx);

        return operator()(image, opt);
    }
    virtual double operator()(array_2d & test_image, double im_scale, MultiConvolutionKernel kernel_approx) const
    {
        double w[2] = { 1.0 / im_scale, 1 };
        return operator()(test_image, zeros, w, kernel_approx);
    }
    virtual double operator()(array_2d & image, MultiConvolutionKernel kernel_approx) const
    {
        return operator()(image, 1, kernel_approx);
    }
    virtual std::string name() const = 0;
    virtual void log(logger & log_write) const
    {
        log_write("test operator", name());
    }
    virtual ~t_func() {}
};
const double t_func::zeros[2] = { 0, 0 };

struct gsmaa_f : public t_func
{
    double sigma;
    gsmaa_f(double s) : sigma(s) {}
    std::string name() const {
        return "gaussianSmoothMultiArray";
    }
    double operator()(array_2d & test_image,  const options_2d & opt) const
    {
        vigra::gaussianSmoothMultiArray(vigra::srcMultiArrayRange(test_image),
                                        vigra::destMultiArray(test_image),
                                        sigma, opt);
        return sigma;
    }
};

struct ggma_f : public t_func
{
    double sigma;
    int axis;
    ggma_f(double s, int a) : sigma(s), axis(a) {}
    std::string name() const
    {
        return "gaussianGradientMultiArray, axis " + asString(axis);
    }
    double operator()(array_2d & test_image,  const options_2d & opt) const {

        grad_2d grad_data(test_image.shape());
        vigra::gaussianGradientMultiArray(vigra::srcMultiArrayRange(test_image),
                                          vigra::destMultiArray(grad_data),
                                          sigma, opt);
        // copy gradient #axis to image:
        typedef grad_2d::value_type value_type;
        typedef vigra::VectorComponentAccessor<value_type> accessor;
        vigra::copyMultiArray(vigra::srcMultiArrayRange(grad_data, accessor(axis)),
                              vigra::destMultiArray(test_image));
        return sigma;
    }
};

struct sgma_f : public t_func {
    int axis;
    sgma_f(int a) : axis(a) {}
    std::string name() const {
        return "symmetricGradientMultiArray, axis " + asString(axis);
    }
    double operator()(array_2d & test_image,  const options_2d & opt) const
    {
        grad_2d grad_data(test_image.shape());
        vigra::symmetricGradientMultiArray(
                                          vigra::srcMultiArrayRange(test_image),
                                          vigra::destMultiArray(grad_data),
                                          opt);

        // copy gradient #axis to image:
        typedef grad_2d::value_type value_type;
        typedef vigra::VectorComponentAccessor<value_type> accessor;
        vigra::copyMultiArray(vigra::srcMultiArrayRange(grad_data, accessor(axis)),
                              vigra::destMultiArray(test_image));
        return 0;
    }
};

struct lgma_f : public t_func
{
    double sigma;
    lgma_f(double s) : sigma(s) {}
    std::string name() const {
        return "laplacianOfGaussianMultiArray";
    }
    double operator()(array_2d & test_image,  const options_2d & opt) const
    {
        array_2d dest_image(test_image);
        vigra::laplacianOfGaussianMultiArray(
                                          vigra::srcMultiArrayRange(test_image),
                                          vigra::destMultiArray(dest_image),
                                          sigma, opt);
        // copy result to image:
        vigra::copyMultiArray(vigra::srcMultiArrayRange(dest_image),
                              vigra::destMultiArray(test_image));
        return sigma;
    }
};

struct hgma_f : public t_func
{
    double sigma;
    int entry;
    hgma_f(double s, int e) : sigma(s), entry(e) {}
    std::string name() const
    {
        return "hessianOfGaussianMultiArray, entry " + asString(entry);
    }
    double operator()(array_2d & test_image,  const options_2d & opt) const
    {
        symm_2d hess_data(test_image.shape());
        vigra::hessianOfGaussianMultiArray(
                                          vigra::srcMultiArrayRange(test_image),
                                          vigra::destMultiArray(hess_data),
                                          sigma, opt);
        // copy hessian entry to image:
        typedef symm_2d::value_type value_type;
        typedef vigra::VectorComponentAccessor<value_type> accessor;
        vigra::copyMultiArray(vigra::srcMultiArrayRange(hess_data, accessor(entry)),
                              vigra::destMultiArray(test_image));
        return sigma;
    }
};

struct stma_f : public t_func
{
    double inner;
    double outer;
    int entry;
    stma_f(double in, double out, int e) : inner(in), outer(out), entry(e) {}
    std::string name() const
    {
        return "structureTensorMultiArray, entry " + asString(entry);
    }
    double operator()(array_2d & test_image,  const options_2d & opt) const
    {

        symm_2d st_data(test_image.shape());
        vigra::structureTensorMultiArray(vigra::srcMultiArrayRange(test_image),
                                         vigra::destMultiArray(st_data),
                                         inner, outer, opt);
        // copy st entry to image:
        typedef symm_2d::value_type value_type;
        typedef vigra::VectorComponentAccessor<value_type> accessor;
        vigra::copyMultiArray(vigra::srcMultiArrayRange(st_data, accessor(entry)),
                              vigra::destMultiArray(test_image));
        return std::sqrt(inner * inner + outer * outer);
    }
};

const t_func * new_test_alloc(double sigma, double outer, int test_nr)
{
    switch (test_nr)
    {
    case 0:
            return new gsmaa_f(sigma);
    case 1:
            return new ggma_f(sigma, 0);
    case 2:
            return new ggma_f(sigma, 1);
    case 11:
            return new sgma_f(0);
    case 12:
            return new sgma_f(1);
    case 20:
            return new lgma_f(sigma);
    case 21:
            return new hgma_f(sigma, 0);
    case 22:
            return new hgma_f(sigma, 1);
    case 23:
            return new hgma_f(sigma, 2);
    case 31:
            return new stma_f(sigma, outer, 0);
    case 32:
            return new stma_f(sigma, outer, 1);
    case 33:
            return new stma_f(sigma, outer, 2);
        }
    return 0;
}

const t_func * new_test(double sigma, double outer, int test_nr,
                                                              logger & log_name)
{
    const t_func *const test = new_test_alloc(sigma, outer, test_nr);
    if (test)
        test->log(log_name);
    else
        log_name("[no test operator found in new_test_alloc()]");
    return test;
}

struct cmp_double
{
    double tol;
    cmp_double(double t) : tol(t) {}
    void check(double val) const
    {
        shouldMsg(!tol || (std::abs(val) <= tol),
                  ("relative difference above tolerance: "
                  "accepted = " + asString(tol) +
                  ", actual = " + asString(val)).c_str());
    }
};

struct cmp_data
{
    bool write_im;
    cmp_double global_tol;
    cmp_double local_tol;
    cmp_data(bool w, double glob, double loc)
        : write_im(w), global_tol(glob), local_tol(loc) {}
    operator bool() const
    {
        return write_im;
    }
};

void test_compare(const array_2d & x1_image_b, const array_2d & res_image,
                  double sigma_pix, const cmp_data & cmp,
                  logger & log_write)
{
    double mm_delta = min_max_delta(x1_image_b);

    double dist = compare(delta_dist(1), "delta_dist",
                          res_image, x1_image_b, sigma_pix,
                          "delta_dist_im.png", cmp);
    double rel = compare(rel_dist(1), "rel_dist", res_image, x1_image_b, sigma_pix,
                         "rel_dist_im.png", cmp);

    log_write("globally relative error", dist / mm_delta);
    log_write("locally relative error", rel);
    cmp.global_tol.check(dist / mm_delta);
    cmp.local_tol.check(rel);
}

unsigned resized(double scale, int size)
{
    return static_cast<unsigned>(std::floor(1 + scale * (size - 1)));
}

shape_2d resized_shape(const vigra::ImageImportInfo & size_info, double scale_x,
                       double scale_y)
{
    return shape_2d(resized(scale_x, size_info.width()),
                    resized(scale_y, size_info.height()));
}

void test_upscaled(void (*resize)(const array_2d &, array_2d &),
                   const vigra::ImageImportInfo & size_info,
                   const array_2d & test_image,
                   const t_func & test_f,
                   double im_scale, const cmp_data & cmp,
                   logger & log_write, logger & log_name,
                   MultiConvolutionKernel kernel_approx)
{

    log_name("upscaled test");
    if (log_name)
        return;
    // image inflated by im_scale in x direction:
    array_2d x_scaled_image(resized_shape(size_info, im_scale, 1));
    (*resize)(test_image, x_scaled_image);

    write_out(x_scaled_image, "test_any_scaled.png", cmp);
    double sigma = test_f(x_scaled_image, im_scale, kernel_approx);

    write_out(x_scaled_image, "test_any_scaled_new_f.png", cmp);

    array_2d x1_image_b(test_image);
    test_f(x1_image_b, kernel_approx);
    write_out(x1_image_b, "test_any_old_f.png", cmp);

    // compare with resampled image
    array_2d res_image(test_image);
    (*resize)(x_scaled_image, res_image);
    write_out(res_image, "test_any_res_f.png", cmp);

    test_compare(x1_image_b, res_image, sigma, cmp, log_write);
}

void test_downscaled(void (*resize)(const array_2d &, array_2d &),
                     const vigra::ImageImportInfo & size_info,
                     const array_2d & m2_image_old,
                     const t_func & test_f,
                     double im_scale, const cmp_data & cmp,
                     logger & log_write, logger & log_name,
                     MultiConvolutionKernel kernel_approx)
{
    const double y_scale = 3;
    const double x_scale = im_scale * y_scale;
    log_name("downscaled test (factor " + asString(y_scale) + ")");
    if (log_name)
        return;

    int w = (int(size_info.width() - 1) / int(x_scale)) * int(x_scale) + 1;
    int h = (int(size_info.height() - 1) / int(y_scale)) * int(y_scale) + 1;

    array_2d test_image(shape_2d(w, h));
    (*resize)(m2_image_old, test_image);

    // pre-sample image with gaussian, std. dev. == scale:
    const double std_dev_factor = 1;
    array_2d pre_scale_image(test_image);
    double sigmas[2] = { 0, 0 };
    sigmas[0] = std_dev_factor * x_scale;
    sigmas[1] = std_dev_factor * y_scale;
    vigra::gaussianSmoothMultiArray(test_image, pre_scale_image,
                                    options_2d().stdDev(sigmas));
    // downscale:
    array_2d downscaled_image(resized_shape(size_info,
                        1 / x_scale, 1 / y_scale));
    (*resize)(pre_scale_image, downscaled_image);
    write_out(downscaled_image, "test_downscaled.png", cmp);

    const double step_size[2] = { x_scale, y_scale };

    double sigma = test_f(downscaled_image, sigmas, step_size, kernel_approx);

    write_out(downscaled_image, "test_downscaled_new_f.png", cmp);

    array_2d x1_image_b(test_image);
    test_f(x1_image_b, kernel_approx);

    write_out(x1_image_b, "test_any_old_f.png", cmp);

    // compare with resampled image
    array_2d res_image_b(downscaled_image);
    (*resize)(x1_image_b, res_image_b);

    write_out(res_image_b, "test_any_res_b_f.png", cmp);

    test_compare(res_image_b, downscaled_image, sigma / x_scale, cmp, log_write);
}

void test_any(void (*resize)(const array_2d &, array_2d &),
              const vigra::ImageImportInfo & size_info,
              const array_2d & test_image,
              const t_func & test_f,
              double im_scale,
              const cmp_data & cmp,
              logger & log_write,
              logger & log_name,
              int test_type,
              MultiConvolutionKernel kernel_approx)
{
    if (test_type == 0)
    {
        test_upscaled(resize, size_info, test_image, test_f, im_scale,
                      cmp, log_write, log_name, kernel_approx);
    }
    else if (test_type == 1)
    {
        test_downscaled(resize, size_info, test_image, test_f, im_scale,
                        cmp, log_write, log_name, kernel_approx);
    }
}

typedef const double test_data[9];

struct args
{
    const int argc;
    test_data & argv;
    int pos;
    args(int a_c, test_data & a_v) : argc(a_c), argv(a_v), pos(-1) {}
    template<class X>
    X operator()(X default_value)
    {
        ++pos;
        return static_cast<X>(
                 (argc > pos) ? argv[pos] : static_cast<double>(default_value));
    }
};

std::string perform_test(int argc, test_data & argv,
                         const vigra::ImageImportInfo & import_info,
                         const array_2d & test_image,
                         MultiConvolutionKernel kernel_approx,
                         bool run_test = true)
{
    logger log_name(!run_test);
    logger log_write;
    args cmd_line(argc, argv);
    const double sigma        = cmd_line(2.0);
    const double im_scale     = cmd_line(3.0);
    const int intp_type       = cmd_line(0);
    const bool write_im       = cmd_line(1) == 1;
    const int test_nr         = cmd_line(0);
    const double outer        = cmd_line(sigma);
    const int test_type       = cmd_line(0);
    const double global_tol   = cmd_line(0.0);
    const double local_tol    = cmd_line(0.0);

    const cmp_data cmp(write_im, global_tol, local_tol);

    const t_func *const test = new_test(sigma, outer, test_nr, log_name);
    if (! test)
    {
        shouldMsg(!run_test, ("unknown test number " + asString(test_nr)
                                            + " for new_test_alloc()").c_str());
        return "(unknown test number)";
    }
    log_name("sigma", sigma);
    log_name("image scaling factor", im_scale);
    log_name("outer sigma", outer);
    log_name("global_tol", global_tol);
    log_name("local_tol", local_tol);
    
    if (intp_type == 0)
    {
        log_name("resizing without interpolation");
        test_any(&resize_0, import_info, test_image, *test, im_scale,
                    cmp, log_write, log_name, test_type, kernel_approx);
    }
    else if (intp_type == 3)
    {
        log_name("resizing with spline3");
        test_any(&resize_n<3>, import_info, test_image, *test, im_scale,
                    cmp, log_write, log_name, test_type, kernel_approx);
    }
    else if (intp_type == 5)
    {
        log_name("resizing with spline5");
        test_any(&resize_n<5>, import_info, test_image, *test, im_scale,
                    cmp, log_write, log_name, test_type, kernel_approx);
    }
    if (log_name)
        return log_name.str();
    else
        return log_name.str() + ":\n" + log_write.str();
}

struct scaled_test : public vigra::detail::test_functor<scaled_test>
{
    static const int argc;
    test_data & argv;
    const vigra::ImageImportInfo & import_info;
    const array_2d & test_image;
    MultiConvolutionKernel kernel_approx_;

    scaled_test(test_data & a,  const vigra::ImageImportInfo & ii,
                                                             const array_2d & t, MultiConvolutionKernel kernel_approx)
        : argv(a), import_info(ii), test_image(t), kernel_approx_(kernel_approx) {}
    void operator()()
    {
        // std::cout << perform_test(argc, argv, import_info, test_image)
        // << "\n";
        perform_test(argc, argv, import_info, test_image, kernel_approx_);
    }
    std::string str()
    {
        return "MultiArraySeparableConvolutionScaledTestSuite ["
               + perform_test(argc, argv, import_info, test_image, kernel_approx_, false) + "]";
    }
};
const int scaled_test::argc = sizeof(test_data) / sizeof(double);

test_data tests[] =
{
/*    sigma        write_im     test_type           */
/*         im_scale    test_nr      global_tol      */
/*             intp_type   outer        local_tol   */
    { 2,   3,  5,  0,  0,  0,   0,  0,  0.007 },
    { 1,   3,  0,  0,  1,  0,   0,  0.023,  0 },
    { 1,   3,  0,  0,  2,  0,   0,  0.011,  0 },
    { 4.5, 3,  0,  0,  1,  0,   0,  0.004,  0 },
    { 4.5, 3,  0,  0,  2,  0,   0,  0.002,  0 },
    { 5, 0.99, 5,  0, 11,  0,   0,  0.065,  0 },
    { 5, 0.99, 5,  0, 12,  0,   0,  0.135,  0 },
    { 3,   3,  0,  0, 20,  0,   0,  0.014,  0 },
    { 3,   3,  0,  0, 21,  0,   0,  0.019,  0 },
    { 3,   3,  0,  0, 22,  0,   0,  0.005,  0 },
    { 3,   3,  0,  0, 23,  0,   0,  0.001,  0 },
    { 4,   3,  0,  0, 31,  2,   0,  0.012,  0 },
    { 4,   3,  0,  0, 32,  2,   0,  0.004,  0 },
    { 4,   3,  0,  0, 33,  2,   0,  0.001,  0 },
    { 15,  3,  0,  0, 0,   0,   1,  0,  0.012 },
    { 15,  3,  5,  0, 0,   0,   1,  0,  0.012 },
    { 15,  2,  5,  0, 1,   0,   1,  0.005,  0 },
    { 15,  2,  5,  0, 2,   0,   1,  0.003,  0 },
    { 15,  2,  0,  0, 20,  0,   1,  0.028,  0 },
    { 15,  3,  0,  0, 21,  0,   1,  0.045,  0 },
    { 15,  3,  0,  0, 22,  0,   1,  0.023,  0 },
    { 15,  3,  0,  0, 23,  0,   1,  0.024,  0 },
    { 15,  3,  0,  0, 31,  1.1, 1,  0.025,  0 },
    { 15,  3,  0,  0, 32,  1.1, 1,  0.035,  0 },
    { 15,  3,  0,  0, 33,  1.1, 1,  0.006,  0 },
    { 0,   0,  0,  0,  0,  0,   0,  0,      0 }
};


//--------------------------------------------------------


struct MultiArraySeparableConvolutionScaledTestSuite : public vigra::test_suite
{
    vigra::ImageImportInfo import_info;
    array_2d               test_image;

    MultiArraySeparableConvolutionScaledTestSuite(MultiConvolutionKernel kernel_approx, const char *suit_name = "MultiArraySeparableConvolutionScaledTestSuite")
        : vigra::test_suite(suit_name),
          import_info("oi_single.gif"),
          test_image(shape_2d(import_info.width(), import_info.height()))
        {
            vigra::importImage(import_info, destImage(test_image));

            for (test_data* p = tests; (*p)[0]; ++p)
            {
                scaled_test* test
                                 = new scaled_test(*p, import_info, test_image, kernel_approx);
                add(vigra::create_test_case(*test, test->str().c_str()));
            }
        }
}; // struct MultiArraySeparableConvolutionScaledTestSuite

//--------------------------------------------------------


template <typename ConvolutionKernel>
struct MultiArraySeparableRecursiveConvolutionTest
{

    typedef float PixelType;

    typedef TinyVector<PixelType,3> VectorPixelType;
    typedef MultiArray<3, VectorPixelType> Image3x3;

    typedef MultiArray<3,PixelType> Image3D;
    typedef Image3D::difference_type Size3;

    typedef BasicImage<PixelType> Image2D;
    typedef Image2D::difference_type Size2;


    // - - - - - - - - - - - - - - - - - - - - - - - -

    template<class Image>
    void makeRandom( Image &image )
    {
        typedef typename Image::value_type T;

        int size = image.size();
        for( int k = 0; k < size; ++k ) 
        {
            typedef typename NumericTraits<typename Image::value_type>::isIntegral isIntegral;
            if(isIntegral::value)
                image[k] = (T)randomMT19937().uniformInt(256);
            else
                image[k] = (T)randomMT19937().uniform();
        }
    }

    void testBorder(BorderTreatmentMode border_treatment)
    {
        ConvolutionKernel kernel;
        kernel.initGaussian(10.0);
        MultiArray<1, double> src(100), src_copy(100 + kernel.right() - kernel.left()), res(100), res_copy(100 + kernel.right() - kernel.left());

        makeRandom(src);

        vigra::detail::copyLineWithBorderTreatment(src.traverser_begin(), src.traverser_end(), StandardValueAccessor<double>(),
            src_copy.traverser_begin(), StandardValueAccessor<double>(), 0, 100, kernel.left(), kernel.right(), border_treatment);

        kernel.setBorderTreatment(border_treatment);
        separableConvolveMultiArray(srcMultiArrayRange(src), destMultiArray(res), kernel);

        kernel.setBorderTreatment(BORDER_TREATMENT_ZEROPAD);
        separableConvolveMultiArray(srcMultiArrayRange(src_copy), destMultiArray(res_copy), kernel);

        shouldEqualSequence(res.begin(), res.end(), res_copy.begin() - kernel.left());
    }

    void testBorder(void)
    {
        testBorder(BORDER_TREATMENT_ZEROPAD);
        testBorder(BORDER_TREATMENT_REPEAT);
        testBorder(BORDER_TREATMENT_REFLECT);
        testBorder(BORDER_TREATMENT_WRAP);
    }

    void testScaling(void)
    {
        MultiArray<1, double> src(100), res(100), res_scaled(100);
        makeRandom(src);

        ConvolutionKernel kernel;

        kernel.initGaussian(5.0);
        separableConvolveMultiArray(srcMultiArrayRange(src), destMultiArray(res), kernel);

        kernel.scale(0.5);
        separableConvolveMultiArray(srcMultiArrayRange(src), destMultiArray(res_scaled), kernel);

        for (unsigned int i = 0; i < 100; ++i)
            res_scaled(i) *= 2;

        shouldEqualSequence(res.begin(), res.end(), res_scaled.begin());
    }

    void testGaussian(unsigned int order, double sigma) {
        ConvolutionKernel kernel_iir;
        Kernel1D<double> kernel_fir;
        MultiArray<1, double> src(100), res_iir(100), res_fir(100);

        for (unsigned int i = 0; i < 100; ++i) {
            src[i] = res_iir[i] = res_fir[i] = 0.0;
        }

        src[50] = 1.0;

        kernel_iir.initGaussianDerivative(sigma, order);
        kernel_fir.initGaussianDerivative(sigma, order, 1.0, 9);

        separableConvolveMultiArray(srcMultiArrayRange(src), destMultiArray(res_iir), kernel_iir);
        separableConvolveMultiArray(srcMultiArrayRange(src), destMultiArray(res_fir), kernel_fir);

        double diff;
        double tol = 1e-12;

        if (kernel_iir.order == 2)
            tol = 3e-3;
        else if(kernel_iir.order == 3)
            tol = 1e-3;
        else if(kernel_iir.order == 4)
            tol = 1e-4;

        for (unsigned int i = 0; i < 100; ++i) {
            diff = abs(res_fir[i] - res_iir[i]);
            if (diff > tol) {
                std::ostringstream msg;
                msg << "Assertion failed: Sequence items differ at index " << i << " abs(" << res_iir[i] << " - " << res_fir[i] << ") = " << diff << " > " << tol;
                vigra_fail(msg.str());
            }
        }
    }

    void testSmoothing(void) {
        testGaussian(0, 5.0);
        testGaussian(0, 10.0);
    }

    void test1stDeriv(void) {
        testGaussian(1, 5.0);
        testGaussian(1, 10.0);
    }

    void test2ndDeriv(void) {
        testGaussian(2, 5.0);
        testGaussian(2, 10.0);       
    }

    MultiArraySeparableRecursiveConvolutionTest()
    {
    }
};


template <typename ConvolutionKernel>
struct MultiArraySeparableRecursiveConvolutionTestSuite
: public vigra::test_suite
{
    MultiArraySeparableRecursiveConvolutionTestSuite(const char *name = "MultiArraySeparableRecursiveConvolutionTestSuite")
        : vigra::test_suite(name)
        {
                add( testCase( &MultiArraySeparableRecursiveConvolutionTest<ConvolutionKernel>::testBorder ) );
                add( testCase( &MultiArraySeparableRecursiveConvolutionTest<ConvolutionKernel>::testScaling ) );
                add( testCase( &MultiArraySeparableRecursiveConvolutionTest<ConvolutionKernel>::testSmoothing ) );
                add( testCase( &MultiArraySeparableRecursiveConvolutionTest<ConvolutionKernel>::test1stDeriv ) );
                add( testCase( &MultiArraySeparableRecursiveConvolutionTest<ConvolutionKernel>::test2ndDeriv ) );
    }
}; // struct MultiArraySeparableRecursiveConvolutionTestSuite

//--------------------------------------------------------


struct MultiArraySeparableRecursiveConvolutionImageTest
: public vigra::test_suite
{
    typedef MultiArray<2, double> Image;
    typedef vigra::MultiArrayShape<2>::type shape_2d;
    typedef RecursiveConvolutionKernel<double, detail::deriche_2_tag> DericheKernel2nd;
    typedef RecursiveConvolutionKernel<double, detail::deriche_3_tag> DericheKernel3rd;
    typedef RecursiveConvolutionKernel<double, detail::deriche_4_tag> DericheKernel4th;

    template <class RecursiveConvolutionKernel>
    static void convolve(Image &srcimg, Image &destimg, unsigned order, double sigma)
    {
        RecursiveConvolutionKernel kernel;
        kernel.initGaussianDerivative(sigma, order);
        separableConvolveMultiArray(srcMultiArrayRange(srcimg), destMultiArray(destimg), kernel);
    }

    template<class RecursiveConvolutionKernel>
    static void test_equal(const char * src, const char *dst, unsigned order, double sigma)
    {
        ImageImportInfo src_info(src);

        Image srcimg(shape_2d(src_info.width(), src_info.height()));
        Image dstimg(shape_2d(src_info.width(), src_info.height()));
        Image resimg(shape_2d(src_info.width(), src_info.height()));

        importImage(src_info, destImage(srcimg));
        importImage(ImageImportInfo(dst), destImage(dstimg));
        resimg.init(42.42);

        convolve<RecursiveConvolutionKernel>(srcimg, resimg, order, sigma);
        //exportImage(srcImageRange(resimg), vigra::ImageExportInfo(dst));

        shouldEqualSequence(dstimg.data(), dstimg.data()+dstimg.size(), resimg.data());
    }

    static void test_deriche(unsigned deriche_order, unsigned deriv_order, double sigma)
    {
        const char *srcname = "lenna.gif";
        std::ostringstream dstname;

        dstname << "deriche_" << deriche_order << "_" << deriv_order << "_" << sigma << ".xv";

        switch (deriche_order) {
            case 2: test_equal<DericheKernel2nd>(srcname, dstname.str().c_str(), deriv_order, sigma); break;
            case 3: test_equal<DericheKernel3rd>(srcname, dstname.str().c_str(), deriv_order, sigma); break;
            case 4: test_equal<DericheKernel4th>(srcname, dstname.str().c_str(), deriv_order, sigma); break;
            default: vigra_fail("MultiArraySeparableRecursiveConvolutionImageTest::test_deriche: Invalid deriv_order.");
        } 
    }

    static void test_deriche_nth(unsigned deriche_order)
    {
        for (unsigned deriv_order = 0; deriv_order < 3; ++deriv_order) {
            test_deriche(deriche_order, deriv_order, 5.0);
            test_deriche(deriche_order, deriv_order, 10.0);
            test_deriche(deriche_order, deriv_order, 14.0);
        }
    }

    static void test_deriche_2nd(void)
    {
        test_deriche_nth(2);
    }

    static void test_deriche_3rd(void)
    {
        test_deriche_nth(3);
    }

    static void test_deriche_4th(void)
    {
        test_deriche_nth(4);
    }

    MultiArraySeparableRecursiveConvolutionImageTest(const char *name = "MultiArraySeparableRecursiveConvolutionImageTest")
        : vigra::test_suite(name)
        {
                add( testCase( &test_deriche_2nd) );
                add( testCase( &test_deriche_3rd) );
                add( testCase( &test_deriche_4th) );
    }
};



int main(int argc, char ** argv)
{
    int failed = 0;

    // run the multi-array separable convolution test suite
    MultiArraySeparableConvolutionTestSuite test_mascts;
    failed += test_mascts.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_mascts.report() << std::endl;

    // run the multi-array separable scaled-convolution test suite for FIR filters
    MultiArraySeparableConvolutionScaledTestSuite test_mascsts_fir(MULTI_CONVOLUTION_KERNEL_FIR);
    failed += test_mascsts_fir.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_mascsts_fir.report() << std::endl;

    // run the multi-array separable scaled-convolution test suite for Deriche IIR filters
    MultiArraySeparableConvolutionScaledTestSuite test_mascsts_iir_deriche(MULTI_CONVOLUTION_KERNEL_IIR_DERICHE, "MultiArraySeparableConvolutionScaledTestSuiteRecursiveDeriche");
    failed += test_mascsts_iir_deriche.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_mascsts_iir_deriche.report() << std::endl;

    // run the multi-array separable recursive scaled-convolution test suite for Deriche filters
    MultiArraySeparableRecursiveConvolutionTestSuite<RecursiveConvolutionKernel<double, detail::deriche_2_tag>> test_masrcts_deriche2nd("MultiArraySeparableRecursiveConvolutionTestSuiteDeriche2ndOrder");
    failed += test_masrcts_deriche2nd.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_masrcts_deriche2nd.report() << std::endl;

    // run the multi-array separable recursive scaled-convolution test suite for Deriche filters
    MultiArraySeparableRecursiveConvolutionTestSuite<RecursiveConvolutionKernel<double, detail::deriche_3_tag>> test_masrcts_deriche3rd("MultiArraySeparableRecursiveConvolutionTestSuiteDeriche3rdOrder");
    failed += test_masrcts_deriche3rd.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_masrcts_deriche3rd.report() << std::endl;

    // run the multi-array separable recursive scaled-convolution test suite for Deriche filters
    MultiArraySeparableRecursiveConvolutionTestSuite<RecursiveConvolutionKernel<double, detail::deriche_4_tag>> test_masrcts_deriche4th("MultiArraySeparableRecursiveConvolutionTestSuiteDeriche4thOrder");
    failed += test_masrcts_deriche4th.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_masrcts_deriche4th.report() << std::endl;

    // run the multi-array separable recursive image test suite
    MultiArraySeparableRecursiveConvolutionImageTest test_recursive_image;
    failed += test_recursive_image.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test_recursive_image.report() << std::endl;


    return (failed != 0);
}
