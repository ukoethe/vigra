// -*- c++ -*-
// $Id$

#include "unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/convolution.hxx" 
#include "vigra/navigator.hxx"
#include "vigra/functorexpression.hxx"

using namespace vigra;
using namespace vigra::functor;

struct MultiArrayPointoperatorsTest
{

    typedef float PixelType;
    typedef MultiArray<3,PixelType> Image3D;
    typedef MultiArrayView<3,PixelType> View3D;
    typedef Image3D::difference_type Size3;

    Image3D img;

    MultiArrayPointoperatorsTest()
    : img(Size3(5,4,3))
    {
        unsigned int i;
        PixelType c = 0.1;
        for(i=0; i<img.elementCount(); ++i, ++c)
            img.data()[i] = c;
    }

    void testInit()
    {
        Image3D res(img.shape());
        const Image3D::value_type ini = 1.1;
        should(res.shape() == Size3(5,4,3));

        initMultiArray(destMultiArrayRange(res), ini);
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), ini);
    }

    void testCopy()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        copyMultiArray(srcMultiArrayRange(img), destMultiArray(res));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z));
    }

    void testCopyOuterExpansion()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        
        copyMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,0,0));
    }

    void testCopyInnerExpansion()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        
        copyMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(0,0,z));
    }

    void testTransform()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        transformMultiArray(srcMultiArrayRange(img), destMultiArray(res),
                            Arg1() + Arg1());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z));
    }

    void testTransformOuterExpand()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        
        transformMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res),
                            Arg1() + Arg1());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,0,0));
    }

    void testTransformInnerExpand()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        
        transformMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res),
                            Arg1() + Arg1());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(0,0,z));
    }

    void testTransformOuterReduce()
    {
        Image3D res(Size3(5,1,1));
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        transformMultiArray(srcMultiArrayRange(img), destMultiArrayRange(res),
                            reduceFunctor(Arg1() + Arg2(), 0.0));
        
        int x,y,z;
        for(x=0; x<img.shape(0); ++x)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(z=0; z<img.shape(2); ++z)
                    sum += img(x,y,z);
            shouldEqual(res(x,0,0), sum);
        }
    }

    void testTransformInnerReduce()
    {
        Image3D res(Size3(1,1,3));
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        transformMultiArray(srcMultiArrayRange(img), destMultiArrayRange(res),
                            reduceFunctor(Arg1() + Arg2(), 0.0));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    sum += img(x,y,z);
            shouldEqual(res(0,0,z), sum);
        }
    }

    void testCombine2()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArray(img), 
                              destMultiArray(res),
                              Arg1() + Arg2());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z));
    }

    void testCombine2OuterExpand()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(x,0,0));

        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(x,0,0));

        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 3.0*img(x,0,0));
    }

    void testCombine2InnerExpand()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(0,0,z));

        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(0,0,z));

        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 3.0*img(0,0,z));
    }

    void testCombine2OuterReduce()
    {
        Image3D res(Size3(5,1,1));
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
                              reduceFunctor(Arg1() + Arg2() + Arg3(), 0.0));
        
        int x,y,z;
        for(x=0; x<img.shape(0); ++x)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(z=0; z<img.shape(2); ++z)
                    sum += img(x,y,z);
            shouldEqual(res(x,0,0), 2.0*sum);
        }
    }

    void testCombine2InnerReduce()
    {
        Image3D res(Size3(1,1,3));
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
                              reduceFunctor(Arg1() + Arg2() + Arg3(), 0.0));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    sum += img(x,y,z);
            shouldEqual(res(0,0,z), 2.0*sum);
        }
    }

    void testCombine3()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineThreeMultiArrays(srcMultiArrayRange(img), 
                                srcMultiArray(img), srcMultiArray(img), 
                                destMultiArray(res),
                                Arg1() + Arg2() + Arg3());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 3.0*img(x,y,z));
    }
};

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

  void makeWedge( Image3D &image )
  {
    const Size3 size = image.shape();
    const int width = size[0];
    const int height = size[1];
    const int depth = size[2];
    for( int z = 0; z < depth; ++z ) 
    {
      for( int y = 0; y < height; ++y ) 
      {
        for( int x = 0; x < width; ++x ) 
        {
          const Image3D::value_type val = x + y + z;
          image( x, y, z ) = val;
        }
      }
    }
  }


  void makeBox( Image3D &image )
  {
    const int b = 8;
    const Size3 size = image.shape();
    const int width = size[0];
    const int height = size[1];
    const int depth = size[2];
    for( int z = 0; z < depth; ++z ) 
    {
      for( int y = 0; y < height; ++y ) 
      {
        for( int x = 0; x < width; ++x ) 
        {
        
          Image3D::value_type val = 80;

          if( (x > b) && x < (width-b) &&
              (y > b) && y < (height-b) &&
              (z > b)  && z < (depth-b) ) 
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
    Image3D d1( src.size() );
    Image3D dn( src.size() );
    Image3D dest3( src.size() );

    int depth = src.size()[2];

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
 
      shouldEqualSequence(  dn.begin(), dn.end(),
                            d1.begin() );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - -

  void test_1DValidityB( const Image3D &src, float ksize )
  {
    Image3D dst1( src.size() );
    Image3D dst2( src.size() );

    const int depth = src.size()[2];

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

    shouldEqualSequence(  dst1.begin(), dst1.end(),
                          dst2.begin() );


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

    shouldEqualSequence(  dst1.begin(), dst1.end(),
                          dst2.begin() );
  }

  void test_2DValidity( const Image3D &src, float ksize )
  {
    Image3D d2( src.size() );
    Image3D dn( src.size() );

    int depth = src.size()[2];

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

    shouldEqualSequence( dn.begin(), dn.end(),
                         d2.begin() );
  }


  // - - - - - - - - - - - - - - - - - - - - - - - -

  void test_inplacenessN( const Image3D &src, float ksize )
  {
    Image3D da( src.size() );
    Image3D db( src.size() );

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
    Image3D da( src.size() );
    Image3D db( src.size() );

    Kernel1D<float> kernel;
    if( ! useDerivative )
      kernel.initGaussian( ksize );
    else
      kernel.initGaussianDerivative( ksize, 1 );


    for( int i = 0; i < 3; ++i ) {
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

      shouldEqualSequence( da.begin(), da.end(),
                           db.begin() );
    }
  }


  // - - - - - - - - - - - - - - - - - - - - - - - -

  void test_gradient1( const Image3D &base, bool useGaussian )
  {
    const double sigma = kernelSize/2;
    const int b = useGaussian ? int( 0.5 + 3*sigma ) : 1;
    Image3D src( base.size() );
    Image3x3 grad( src.size() );
    makeWedge( src );

    if( ! useGaussian )
      symmetricGradientMultiArray( srcMultiArrayRange(src),
                                   destMultiArray(grad) );
    else
      gaussianGradientMultiArray( srcMultiArrayRange(src),
                                  destMultiArray(grad), sigma );

    Image3x3::value_type v;
    v[0] = 1; v[1] = 1; v[2] = 1;
    const float v2 = dot(v,v);

    const Size3 size = src.shape();
    const int width = size[0];
    const int height = size[1];
    const int depth = size[2];
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


  void test_gradient_magnitude( const Image3D &src )
  {
    // just a test for mere compileability
    Image3D dst( src.size() );
    Image3x3 grad( src.size() );
    symmetricGradientMultiArray( srcMultiArrayRange(src),
                                 destMultiArray(grad) );

    transformMultiArray( srcMultiArrayRange(grad),
                          destMultiArray(dst),
                          VectorNormFunctor<Image3x3::value_type>() );
  }

  //--------------------------------------------

  const Size3 size;
  Image3D srcImage;
  const float kernelSize;

  MultiArraySeparableConvolutionTest()
    : size( 60, 70, 50 ),
      srcImage( size ),
      kernelSize( 1.8 )
  {
    makeBox( srcImage );
  }


  void test_Valid1() {
    test_1DValidity( srcImage, kernelSize );
  }

  void test_Valid3() {
   test_1DValidityB( srcImage, kernelSize );
  }

  void test_Valid2() {
    test_2DValidity( srcImage, kernelSize );
  }

  void test_InplaceN() {
    test_inplacenessN( srcImage, kernelSize );
  }

  void test_Inplace1() {
    test_inplaceness1( srcImage, kernelSize, false );
    test_inplaceness1( srcImage, kernelSize, true );
  }

  void test_gradient1() {
    test_gradient1( srcImage, false );
    test_gradient1( srcImage, true );
  }

  void test_gradient_magnitude() {
    test_gradient_magnitude( srcImage );
  }

};        //-- struct MultiArraySeparableConvolutionTest

//--------------------------------------------------------

struct MultiArrayPointOperatorsTestSuite
: public vigra::test_suite
{
  MultiArrayPointOperatorsTestSuite()
    : vigra::test_suite("MultiArrayPointOperatorsTestSuite")
    {
        // add( testCase( &MultiArrayTest::test_default_ctor ) );
        add( testCase( &MultiArrayPointoperatorsTest::testInit ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCopy ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCopyOuterExpansion ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCopyInnerExpansion ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransform ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformOuterExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformInnerExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformOuterReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformInnerReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2 ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2OuterExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2InnerExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2OuterReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2InnerReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine3 ) );
    }
}; // struct MultiArrayPointOperatorsTestSuite

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
        add( testCase( &MultiArraySeparableConvolutionTest::test_gradient_magnitude ) );
    }
}; // struct MultiArraySeparableConvolutionTestSuite

//--------------------------------------------------------

int main()
{
  // run the multi-array point operator test suite
  MultiArrayPointOperatorsTestSuite test0;
  int failed = test0.run();
  std::cout << test0.report() << std::endl;
  // run the multi-array separable convolution test suite
  MultiArraySeparableConvolutionTestSuite test1;
  failed += test1.run();
  std::cout << test1.report() << std::endl;
  return (failed != 0);
}


