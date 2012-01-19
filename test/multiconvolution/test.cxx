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

#include "unittest.hxx"
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
            separableConvolveMultiArray(srcMultiArrayRange(srcImage), destMultiArray(subarray), 
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
                add( testCase( &MultiArraySeparableConvolutionTest::testSmoothing ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_gradient1 ) );
                add( testCase( &MultiArraySeparableConvolutionTest::test_laplacian ) );
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
    vigra::inspectMultiArray(srcMultiArrayRange(image), minmax);
    return minmax.max - minmax.min;
}

template<class IM>
double max(const IM & image)
{
    vigra::FindMinMax<double> minmax;
    vigra::inspectMultiArray(srcMultiArrayRange(image), minmax);
    return minmax.max;
}

template<class IM>
double min(const IM & image)
{
    vigra::FindMinMax<double> minmax;
    vigra::inspectMultiArray(srcMultiArrayRange(image), minmax);
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
    vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(a),
                                               destMultiArrayRange(b),
                                               my_spline());
}

struct t_func
{
    static const double zeros[2]; // = { 0, 0 };
    virtual double operator()(array_2d & image, const options_2d & opt) const = 0;
    virtual double operator()(array_2d & image,
                              const double sigma_d[2],
                              const double step_size[2]) const
    {
        // test various ways to create option setter functions.
        options_2d opt;
        opt.resolutionStdDev(sigma_d).stepSize(step_size);
        std::vector<double> vstep;
        vstep.push_back(step_size[0]);
        vstep.push_back(step_size[1]);
        opt.resolutionStdDev(sigma_d).stepSize((const double *const)step_size);
        opt.resolutionStdDev(sigma_d).stepSize(vstep);

        return operator()(image, opt);
    }
    virtual double operator()(array_2d & test_image, double im_scale) const
    {
        double w[2] = { 1.0 / im_scale, 1 };
        return operator()(test_image, zeros, w);
    }
    virtual double operator()(array_2d & image) const
    {
        return operator()(image, 1);
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
                   logger & log_write, logger & log_name)
{

    log_name("upscaled test");
    if (log_name)
        return;
    // image inflated by im_scale in x direction:
    array_2d x_scaled_image(resized_shape(size_info, im_scale, 1));
    (*resize)(test_image, x_scaled_image);

    write_out(x_scaled_image, "test_any_scaled.png", cmp);
    double sigma = test_f(x_scaled_image, im_scale);

    write_out(x_scaled_image, "test_any_scaled_new_f.png", cmp);

    array_2d x1_image_b(test_image);
    test_f(x1_image_b);
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
                     logger & log_write, logger & log_name)
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
    vigra::gaussianSmoothMultiArray(vigra::srcMultiArrayRange(test_image),
                                    vigra::destMultiArray(pre_scale_image),
                                    options_2d().stdDev(sigmas));
    // downscale:
    array_2d downscaled_image(resized_shape(size_info,
                        1 / x_scale, 1 / y_scale));
    (*resize)(pre_scale_image, downscaled_image);
    write_out(downscaled_image, "test_downscaled.png", cmp);

    const double step_size[2] = { x_scale, y_scale };

    double sigma = test_f(downscaled_image, sigmas, step_size);

    write_out(downscaled_image, "test_downscaled_new_f.png", cmp);

    array_2d x1_image_b(test_image);
    test_f(x1_image_b);

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
              int test_type)
{
    if (test_type == 0)
    {
        test_upscaled(resize, size_info, test_image, test_f, im_scale,
                      cmp, log_write, log_name);
    }
    else if (test_type == 1)
    {
        test_downscaled(resize, size_info, test_image, test_f, im_scale,
                        cmp, log_write, log_name);
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
                    cmp, log_write, log_name, test_type);
    }
    else if (intp_type == 3)
    {
        log_name("resizing with spline3");
        test_any(&resize_n<3>, import_info, test_image, *test, im_scale,
                    cmp, log_write, log_name, test_type);
    }
    else if (intp_type == 5)
    {
        log_name("resizing with spline5");
        test_any(&resize_n<5>, import_info, test_image, *test, im_scale,
                    cmp, log_write, log_name, test_type);
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

    scaled_test(test_data & a,  const vigra::ImageImportInfo & ii,
                                                             const array_2d & t)
        : argv(a), import_info(ii), test_image(t) {}
    void operator()()
    {
        // std::cout << perform_test(argc, argv, import_info, test_image)
        // << "\n";
        perform_test(argc, argv, import_info, test_image);
    }
    std::string str()
    {
        return "MultiArraySeparableConvolutionScaledTestSuite ["
               + perform_test(argc, argv, import_info, test_image, false) + "]";
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

    MultiArraySeparableConvolutionScaledTestSuite()
        : vigra::test_suite("MultiArraySeparableConvolutionScaledTestSuite"),
          import_info("oi_single.gif"),
          test_image(shape_2d(import_info.width(), import_info.height()))
        {
            vigra::importImage(import_info, destImage(test_image));

            for (test_data* p = tests; (*p)[0]; ++p)
            {
                scaled_test* test
                                 = new scaled_test(*p, import_info, test_image);
                add(vigra::create_test_case(*test, test->str().c_str()));
            }
        }
}; // struct MultiArraySeparableConvolutionScaledTestSuite

//--------------------------------------------------------

int main(int argc, char ** argv)
{
    // run the multi-array separable convolution test suite
    MultiArraySeparableConvolutionTestSuite test1;
    int failed = test1.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test1.report() << std::endl;

    // run the multi-array separable scaled-convolution test suite
    MultiArraySeparableConvolutionScaledTestSuite test2;
    failed += test2.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test2.report() << std::endl;

    return (failed != 0);
}
