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

#include "unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/convolution.hxx" 
#include "vigra/navigator.hxx"
#include "vigra/random.hxx"

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
        int size = image.size();
        for( int k = 0; k < size; ++k ) 
        {
            typedef typename NumericTraits<typename Image::value_type>::isIntegral isIntegral;
            if(isIntegral::value)
                image[k] = randomMT19937().uniformInt(256);
            else
                image[k] = randomMT19937().uniform();
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

            vigra::convolveMultiArrayOneDimension( srcMultiArrayRange(src),
                                                    destMultiArray(da),
                                                    d,
                                                    kernel );

            copyMultiArray(srcMultiArrayRange(src), destMultiArray(db));

            vigra::convolveMultiArrayOneDimension( srcMultiArrayRange(db),
                                                    destMultiArray(db),
                                                    d,
                                                    kernel );

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
            symmetricGradientMultiArray( srcMultiArrayRange(src), destMultiArray(grad) );
        else
            gaussianGradientMultiArray( srcMultiArrayRange(src), destMultiArray(grad), sigma );

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

        MultiArray<2, double > src(shape), mgrad(shape);
        MultiArray<2, TinyVector<double, 2> > grad(shape);
        BasicImage<double> rmgrad(shape[0], shape[1]);
        
        makeRandom(src);
        
        gaussianGradientMagnitude(srcImageRange(src), destImage(rmgrad), 1.0);
        gaussianGradientMultiArray(srcMultiArrayRange(src), destMultiArray(grad), 1.0 );
        transformMultiArray(srcMultiArrayRange(grad), destMultiArray(mgrad), norm(Arg1()));

        shouldEqualSequence(mgrad.data(), mgrad.data()+size, rmgrad.data());
    }

    void test_laplacian()
    {
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape), laplacian(shape);
        BasicImage<double> rlaplacian(shape[0], shape[1]);
        
        makeRandom(src);
        
        laplacianOfGaussian(srcImageRange(src), destImage(rlaplacian), 2.0);
        laplacianOfGaussianMultiArray(srcMultiArrayRange(src), destMultiArray(laplacian), 2.0 );

        shouldEqualSequenceTolerance(laplacian.data(), laplacian.data()+size, rlaplacian.data(), 1e-12);
    }

    void test_hessian()
    {
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape);
        MultiArray<2, TinyVector<double, 3> > hessian(shape);
        BasicImage<TinyVector<double, 3> > rhessian(shape[0], shape[1]);
        
        makeRandom(src);
        
        typedef VectorComponentAccessor<TinyVector<double, 3> > BandAccessor;
        hessianMatrixOfGaussian(srcImageRange(src), 
                                destImage(rhessian, BandAccessor(0)), 
                                destImage(rhessian, BandAccessor(1)), 
                                destImage(rhessian, BandAccessor(2)), 
                                2.0);
        hessianOfGaussianMultiArray(srcMultiArrayRange(src), destMultiArray(hessian), 2.0 );

        TinyVector<double, 3> epsilon(1e-12, 1e-12, 1e-12);
        shouldEqualSequenceTolerance(hessian.data(), hessian.data()+size, rhessian.data(), epsilon);
    }

    void test_structureTensor()
    {
        MultiArrayShape<2>::type shape(30,40);
        int size = shape[0]*shape[1];

        MultiArray<2, double > src(shape);
        MultiArray<2, TinyVector<double, 3> > st(shape);
        BasicImage<TinyVector<double, 3> > rst(shape[0], shape[1]);
        
        makeRandom(src);
        
        typedef VectorComponentAccessor<TinyVector<double, 3> > BandAccessor;
        structureTensor(srcImageRange(src), 
                        destImage(rst, BandAccessor(0)), 
                        destImage(rst, BandAccessor(1)), 
                        destImage(rst, BandAccessor(2)), 
                        1.5, 3.0);
        structureTensorMultiArray(srcMultiArrayRange(src), destMultiArray(st), 1.5, 3.0 );

        TinyVector<double, 3> epsilon(1e-12, 1e-12, 1e-12);
        shouldEqualSequenceTolerance(st.data(), st.data()+size, rst.data(), epsilon);
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
                add( testCase( &MultiArraySeparableConvolutionTest::test_gradient1 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_laplacian ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_hessian ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_structureTensor ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_gradient_magnitude ) );
        }
}; // struct MultiArraySeparableConvolutionTestSuite

//--------------------------------------------------------

int main(int argc, char ** argv)
{
    // run the multi-array separable convolution test suite
    MultiArraySeparableConvolutionTestSuite test1;
    int failed = test1.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test1.report() << std::endl;
    return (failed != 0);
}


