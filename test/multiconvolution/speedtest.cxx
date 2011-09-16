// -*- c++ -*-
// $Id$

#include "unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_pointoperators.hxx"
//#include "vigra/multi_convolution.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/convolution.hxx" 
#include "vigra/navigator.hxx"
#include "vigra/functorexpression.hxx"

#include <ctime>

#define W 100
#define H 30
#define D 170


using namespace vigra;
using namespace vigra::functor;

namespace Impls
{

  /////////////////////////////////////////////////////////////////////////////////////////

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void 
convolveCDR( SrcIterator si, SrcShape const & shape, SrcAccessor src, 
         DestIterator di, DestAccessor dest, KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp;

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    { // only operate on first dimension here
        SNavigator snav( si, shape, 0 );
        DNavigator dnav( di, shape, 0 );        

        for( ; snav.hasMore(); snav++, dnav++ ) 
        {
      convolveLine( srcIterRange( snav.begin(), snav.end(), src ),
            destIter( dnav.begin(), dest ),
            kernel1d( *kit ) );
        }
        ++kit;
    }

    // operate on further dimensions
    for( int d = 1; d < N; ++d, ++kit ) 
    {
        DNavigator dnav( di, shape, d );

        tmp.resize( shape[d] );

        for( ; dnav.hasMore(); dnav++ ) 
        {
             convolveLine( srcIterRange( dnav.begin(), dnav.end(), dest ),
                           destIter( tmp.begin(), StandardValueAccessor<TmpType>() ),
                           kernel1d( *kit ) );

             // copy temp result to target object
             copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
                       dnav.begin(), dest);
        }
    }
}



template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void convolveCDR( 
        triple<SrcIterator, SrcShape, SrcAccessor> const & source, 
        pair<DestIterator, DestAccessor> const & dest, KernelIterator kit )
{
  convolveCDR( source.first, source.second, source.third, 
           dest.first, dest.second, kit );
}




  /////////////////////////////////////////////////////////////////////////////////////////

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void 
convolveCopyDest( SrcIterator si, SrcShape const & shape, SrcAccessor src, 
          DestIterator di, DestAccessor dest, KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    { // only operate on first dimension here
        SNavigator snav( si, shape, 0 );
        DNavigator dnav( di, shape, 0 );        

        for( ; snav.hasMore(); snav++, dnav++ ) 
        {
      convolveLine( srcIterRange( snav.begin(), snav.end(), src ),
            destIter( tmp.begin(), StandardValueAccessor<TmpType>() ),
            kernel1d( *kit ) );

      // copy temp result to target object
      copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
            dnav.begin(), dest);
        }
        ++kit;
    }

    // operate on further dimensions
    for( int d = 1; d < N; ++d, ++kit ) 
    {
        DNavigator dnav( di, shape, d );

        tmp.resize( shape[d] );

        for( ; dnav.hasMore(); dnav++ ) 
        {
             convolveLine( srcIterRange( dnav.begin(), dnav.end(), dest ),
                           destIter( tmp.begin(), StandardValueAccessor<TmpType>() ),
                           kernel1d( *kit ) );

             // copy temp result to target object
             copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
                       dnav.begin(), dest);
        }
    }
}



template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void convolveCopyDest( 
        triple<SrcIterator, SrcShape, SrcAccessor> const & source, 
        pair<DestIterator, DestAccessor> const & dest, KernelIterator kit )
{
  convolveCopyDest( source.first, source.second, source.third, 
            dest.first, dest.second, kit );
}





  /////////////////////////////////////////////////////////////////////////////////////////

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void 
convolveCopySrc( SrcIterator si, SrcShape const & shape, SrcAccessor src, 
         DestIterator di, DestAccessor dest, KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    { // only operate on first dimension here
        SNavigator snav( si, shape, 0 );
        DNavigator dnav( di, shape, 0 );        

        for( ; snav.hasMore(); snav++, dnav++ ) 
        {
      // copy source to temp 
      copyLine( snav.begin(), snav.end(), src,
            tmp.begin(), StandardValueAccessor<TmpType>() );

      convolveLine( srcIterRange(tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>()),
            destIter(dnav.begin(), dest),
            kernel1d( *kit ) );
        }
        ++kit;
    }

    // operate on further dimensions
    for( int d = 1; d < N; ++d, ++kit ) 
    {
        DNavigator dnav( di, shape, d );

        tmp.resize( shape[d] );

        for( ; dnav.hasMore(); dnav++ ) 
        {
      // copy source to temp 
      copyLine( dnav.begin(), dnav.end(), dest,
            tmp.begin(), StandardValueAccessor<TmpType>() );

      convolveLine( srcIterRange(tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>()),
            destIter(dnav.begin(), dest),
            kernel1d( *kit ) );
        }
    }
}



template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void convolveCopySrc( 
        triple<SrcIterator, SrcShape, SrcAccessor> const & source, 
        pair<DestIterator, DestAccessor> const & dest, KernelIterator kit )
{
  convolveCopySrc( source.first, source.second, source.third, 
            dest.first, dest.second, kit );
}





  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

} //-- namespace Impls


//---------------------------------------------------------------

struct MultiArraySepConvSpeedTest
{
  typedef float PixelType;
  typedef MultiArray<3,PixelType> Image3D;
  typedef Image3D::difference_type Size3;

  const Size3 size;
  Image3D img, dst1, dst2, dst3;
  ArrayVector<vigra::Kernel1D<float> > kernels;

  MultiArraySepConvSpeedTest()
    : size(W,H,D),
      img( size ),
      dst1( size ),
      dst2( size ),
      dst3( size ),
      kernels( 3 )
  {
    std::cout <<"   MultiArraySepConvSpeedTest called" << std::endl;
    const double sigma = 2.3;
    makeBox( img );
    kernels[0].initGaussian( sigma );
    kernels[1].initGaussian( sigma );
    kernels[2].initGaussian( sigma );
  }


  void testCpySrc( const Image3D &src, Image3D &dst )
  {
    Impls::convolveCopySrc( srcMultiArrayRange(src),
                destMultiArray(dst),
                kernels.begin() );
  }


  void testCpyDest( const Image3D &src, Image3D &dst )
  {
    Impls::convolveCopyDest( srcMultiArrayRange(src),
                 destMultiArray(dst),
                 kernels.begin() );
  }


  void testCDR( const Image3D &src, Image3D &dst )
  {
    Impls::convolveCDR( srcMultiArrayRange(src),
            destMultiArray(dst),
            kernels.begin() );
  }


  void testCorrectness()
  {
    testCDR(img,dst1);
    testCpyDest(img,dst2);
    testCpySrc(img,dst3);

    //    shouldEqualSequence( dst1.
  }


  /*
  void speedTest( const char *name, void (MultiArraySepConvSpeedTest::*f)( const Image3D &, Image3D &),
          const Image3D &src, Image3D &dst )
  {
    int t = clock();
    this->*f( src, dst );
    t = clock() - t;

    std::cout << "Timed function: " << name << std::endl << "   = " << t  << std::endl;
  }
  */


#define Speedy(f,name) int t = clock(); f; t = clock() - t;         \
    std::cout << "Timed function: " << name << std::endl << "   = " << t  << std::endl;

  void test1()
  {
    Speedy( testCDR(img,dst1), "testCDR" );
  }

  void test2()
  {
    Speedy( testCpyDest(img,dst2), "testCpyDest" );
  }

  void test3()
  {
    Speedy( testCpySrc(img,dst3), "testCpySrc" );
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

};


struct MultiArraySepConvSpeedTestSuite
: public vigra::test_suite
{
    MultiArraySepConvSpeedTestSuite()
    : vigra::test_suite("MultiArraySeparableConvolutionTestSuite")
    {
        add( testCase( &MultiArraySepConvSpeedTest::test3 ) );
        add( testCase( &MultiArraySepConvSpeedTest::test1 ) );
        add( testCase( &MultiArraySepConvSpeedTest::test2 ) );
        add( testCase( &MultiArraySepConvSpeedTest::testCorrectness ) );
    }
};




int main()
{
  // run the multi-array separable convolutiontestsuite
  MultiArraySepConvSpeedTestSuite test1;
  int failed = test1.run();
  std::cout << test1.report() << std::endl;
  return (failed != 0);
}
