//-- -*- c++ -*-
// $Id$

#include "unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/convolution.hxx" 
#include "vigra/navigator.hxx"

using namespace vigra;


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

  void makeBox( Image3D &image )
  {
    const int b = 8;
    const Size3 size = image.shape();
    const int width = size[0];
    const int height = size[1];
    const int depth = size[2];
    for( int z = 0; z < depth; ++z ) {
      for( int y = 0; y < height; ++y ) {
	for( int x = 0; x < width; ++x ) {
	
	  Image3D::value_type val = 80;

	  if( (x > b) && x < (width-b) &&
	      (y > b) && y < (height-b) &&
	      (z > b)  && z < (depth-b) ) {
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
    std::cerr << "Testing 1D validity... ";

    Image3D d1( src.size() );
    Image3D dn( src.size() );
    Image3D dest3( src.size() );

    int depth = src.size()[2];

    for( int d = 0; d < 3; ++d ) {
      std::cerr << "  Dim: " << d;
      std::vector<vigra::Kernel1D<float> > kernels( 3 );
#if 0
      kernels[d].initGaussian( ksize );
#else
      kernels[d].initGaussianDerivative( ksize, 1 );
#endif

      separableConvolveMultiarray( srcMultiArrayRange(src),
				   destMultiArray(dn),
				   kernels.begin() );

      separableConvolveMultiarray( srcMultiArrayRange(src),
				   destMultiArray(d1),
				   d,
				   kernels[d] );

      shouldEqualSequence(  dn.begin(), dn.end(),
			    d1.begin() );
    }
    std::cerr << std::endl;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - -

  void test_1DValidityB( const Image3D &src, float ksize )
  {
    std::cerr << "Testing 1D validity II... ";

    Image3D dst1( src.size() );
    Image3D dst2( src.size() );

    const int depth = src.size()[2];

    vigra::Kernel1D<float> kernx;
    vigra::Kernel1D<float> kerny;
    kernx.initGaussian( ksize );
    kerny.initGaussianDerivative( ksize, 1 );

    int z;

    //-- x
    separableConvolveMultiarray( srcMultiArrayRange(src),
				 destMultiArray(dst1),
				 0, kernx );

    for( z = 0; z < depth; ++z ) {
      BasicImageView<Image3D::value_type> sslice =
	makeBasicImageView( src.bindOuter(z) );
      BasicImageView<Image3D::value_type> dslice =
	makeBasicImageView( dst2.bindOuter(z) );
      
      vigra::separableConvolveX( srcImageRange(sslice), destImage(dslice),
				 vigra::kernel1d(kernx) );
    }

    shouldEqualSequence(  dst1.begin(), dst1.end(),
			  dst2.begin() );


    //-- y
    separableConvolveMultiarray( srcMultiArrayRange(src),
				 destMultiArray(dst1),
				 1, kerny );

    for( z = 0; z < depth; ++z ) {
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

  // - - - - - - - - - - - - - - - - - - - - - - - -

  void test_2DValidity( const Image3D &src, float ksize )
  {
    std::cerr << "Testing 2D validity... ";

    Image3D d2( src.size() );
    Image3D dn( src.size() );

    int depth = src.size()[2];

    std::vector<vigra::Kernel1D<float> > kernels( 3 );
    kernels[0].initGaussian( ksize );
    kernels[1].initGaussianDerivative( ksize, 1 );

    for( int z = 0; z < depth; ++z ) {
      BasicImageView<Image3D::value_type> sslice =
	makeBasicImageView( src.bindOuter(z) );
      BasicImageView<Image3D::value_type> dslice =
	makeBasicImageView( d2.bindOuter(z) );
      
      vigra::convolveImage( srcImageRange(sslice), destImage(dslice),
			    kernels[0], kernels[1] );
    }

    separableConvolveMultiarray( srcMultiArrayRange(src),
				 destMultiArray(dn),
				 kernels.begin() );

    shouldEqualSequence( dn.begin(), dn.end(),
			 d2.begin() );
    
    std::cerr << std::endl;
  }

  Size3 size;
  Image3D srcImage;
  const float kernelSize;

  MultiArraySeparableConvolutionTest()
    : size( 60, 60, 40 ),
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

};	//-- struct MultiArraySeparableConvolutionTest

//----------------------------------------------------------

struct MultiArraySeparableConvolutionTestSuite
: public vigra::test_suite
{
  MultiArraySeparableConvolutionTestSuite()
    : vigra::test_suite("MultiArraySeparableConvolutionTestSuite")
    {
      //        add( testCase( &MultiArrayTest::test_default_ctor ) );
        add( testCase( &MultiArraySeparableConvolutionTest::test_Valid1 ) );
        add( testCase( &MultiArraySeparableConvolutionTest::test_Valid2 ) );
        add( testCase( &MultiArraySeparableConvolutionTest::test_Valid3 ) );
    }
};	//-- struct MultiArraySeparableConvolutionTestSuite

//----------------------------------------------------------

int main()
{
  // run the multi-array separable convolutiontestsuite
  MultiArraySeparableConvolutionTestSuite test1;
  int failed = test1.run();
  std::cout << test1.report() << std::endl;
  
  return (failed != 0);
}


