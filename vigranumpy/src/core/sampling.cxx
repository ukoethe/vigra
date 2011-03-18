/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpysampling_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/affinegeometry.hxx>
#include <vigra/basicgeometry.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/splines.hxx>
#include <vigra/mathutil.hxx>
#include <vigra/multi_resize.hxx>
#include <vigra/splineimageview.hxx>
#include <vigra/resampling_convolution.hxx>

namespace python = boost::python;

namespace vigra
{
    
template < class PixelType >
NumpyAnyArray pythonResampleImage(NumpyArray<3, Multiband<PixelType> > image, double factor, NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");
    int height,width;
    if(factor < 1.0)
    {
        width = (int) std::ceil(factor * image.shape(0));
        height= (int) std::ceil(factor * image.shape(1));
    }
    else
    {
        width = (int) std::ceil(factor * image.shape(0));
        height= (int) std::ceil(factor * image.shape(1));
    }
    // res.reshapeIfEmpty(MultiArrayShape<3>::type(width,height,image.shape(2)),
                          // "resampleImage(): Output images has wrong dimensions");

    TaggedShape newShape(image.taggedShape());
    newShape[0] = width;
    newShape[1] = height;
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty(newShape,"resampleImage(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resampleImage(srcImageRange(bimage), destImage(bres), factor);
    }

    return res;
}

enum RotationDirection
{
    ROTATE_CW,
    ROTATE_CCW,
    UPSIDE_DOWN
};

template < class PixelType>
NumpyAnyArray pythonFixedRotateImage(NumpyArray<3, Multiband<PixelType> > image, RotationDirection dir, NumpyArray<3,Multiband<PixelType> > res)
{
    int degree=0;
    switch(dir)
    {
    case ROTATE_CW:
        degree=270;
        break;
    case ROTATE_CCW:
        degree=90;
        break;
    case UPSIDE_DOWN:
        degree=180;
        break;
    }
    
    TaggedShape newShape(image.taggedShape());
    if(degree % 180 == 0)
    {
        // res.reshapeIfEmpty(image.shape(),"rotateImageSimple(): Output images has wrong dimensions");
        res.reshapeIfEmpty(newShape,"rotateImageSimple(): Output images has wrong dimensions");
    }
    else
    {
        //FIXME: adjust axistags ? 
        // res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(1),image.shape(0),image.shape(2)),"rotateImage(): Output image has wrong dimensions");
        newShape[0] = image.shape(1);
        newShape[1] = image.shape(0);
        res.reshapeIfEmpty(newShape,"rotateImage(): Output image has wrong dimensions");
    }
    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        rotateImage(srcImageRange(bimage),destImage(bres),degree);
    }
    return res;
}

template < class PixelType>
NumpyAnyArray pythonFreeRotateImageDegree(NumpyArray<3, Multiband<PixelType> > image, double degree,RotationDirection dir, int splineOrder,NumpyArray<3,Multiband<PixelType> > res)
{
    return pythonFreeRotateImageRadiant(image,degree*M_PI/180.0,dir,splineOrder,res);
}

template < class PixelType>
NumpyAnyArray pythonFreeRotateImageRadiant(NumpyArray<3, Multiband<PixelType> > image, double radiant, RotationDirection dir, int splineOrder,NumpyArray<3,Multiband<PixelType> > res)
{
    //reshape, if empty. Otherwise accept res dimensions
    // if(!res.hasData())
        // res.reshapeIfEmpty(image.shape(),"rotateImageRadiant(): Output images has wrong dimensions");
    if(!res.hasData())
        res.reshapeIfEmpty(image.taggedShape(),"rotateImageRadiant(): Output images has wrong dimensions");

    vigra_precondition(res.shape(2)==image.shape(2),"rotateImageRadiant(): number of channels of image and result have to be equal");
    //res.init(NumericTraits< PixelType>::zero());
    if(dir==ROTATE_CW)
        radiant=-radiant;

    //Define the transformation
    linalg::TemporaryMatrix< double >  transform=translationMatrix2D(TinyVector<double,2>(res.shape(0)/2.0,res.shape(1)/2.0))*
        rotationMatrix2DRadians(radiant,TinyVector<double,2>(0.0,0.0))*
        translationMatrix2D(TinyVector<double,2>(-image.shape(0)/2.0,-image.shape(1)/2.0));

    
    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        switch (splineOrder)
        {
        case 0:
            {
                SplineImageView< 0, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 1:
            {
                SplineImageView< 1, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 2:
            {
                SplineImageView< 2, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
            }
            break;
        case 3:
            { 
                SplineImageView< 3, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 4:
            {
                SplineImageView< 4, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 5:
            {
                SplineImageView< 5, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        default:
            PyErr_SetString(PyExc_ValueError, "Spline order not supported.");
            python::throw_error_already_set();
        }
    }
    return res;
}

template < class PixelType>
NumpyAnyArray pythonResizeImageNoInterpolation(NumpyArray<3, Multiband<PixelType> > image, 
                                               python::object destSize,
                                               NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    // res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        // "Output image has wrong dimensions");
    TaggedShape newShape(image.taggedShape());
    newShape[0] = size[0];
    newShape[1] = size[1];
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty( newShape, "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resizeImageNoInterpolation(srcImageRange(bimage),destImageRange(bres));
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonResizeImageLinearInterpolation(NumpyArray<3, Multiband<PixelType> > image,
                                                   python::object destSize,
                                                   NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");

    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");

    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    // res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        // "Output image has wrong dimensions");
    TaggedShape newShape(image.taggedShape());
    newShape[0] = size[0];
    newShape[1] = size[1];
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty( newShape, "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resizeImageLinearInterpolation(srcImageRange(bimage), destImageRange(bres));
    }
    return res;
}

template < class PixelType, int dim >
NumpyAnyArray pythonResizeImageSplineInterpolation(NumpyArray<dim, Multiband<PixelType> > image,
                                                   python::object destSize,
                                                   int splineOrder=3, NumpyArray<dim, Multiband<PixelType> > res=python::object())
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");

    // TinyVector<UInt32,dim> out_shape;
    // if(!res.hasData())
    // {
        // typedef typename MultiArrayShape<dim-1>::type shape;
        // shape size;
        // size=python::extract<shape>(destSize)();
        // for(int ii=0;ii<dim-1;++ii)
            // out_shape[ii]=size[ii];
        // out_shape[dim-1]=image.shape(dim-1);
    // }
    // else
    // {
        // for(int ii=0;ii<dim;++ii)
            // out_shape[ii]=res.shape(ii);
    // }
    // res.reshapeIfEmpty(out_shape, "Output image has wrong dimensions");
    
    typedef typename MultiArrayShape<dim-1>::type OutShape;
    OutShape out_shape;
    if(!res.hasData())
    {
        out_shape=python::extract<OutShape>(destSize)();
    }
    else
    {
        for(int k=0; k<dim-1; ++k)
            out_shape[k] = res.shape(k);
    }
    
    TaggedShape newShape(image.taggedShape());
    newShape = out_shape;
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty( newShape, "Output image has wrong dimensions");

    for(int k=0;k<image.shape(dim-1);++k)
    {
        
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        switch (splineOrder)
        {
            case 1:
            {
                BSpline< 1, double > spline;
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                                               destMultiArrayRange(bres), spline);
                break;
            }
            case 2:
            {
                BSpline< 2, double > spline;
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                    destMultiArrayRange(bres), spline);
                break;
            }
            case 3:
            {
                BSpline< 3, double > spline;
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                    destMultiArrayRange(bres), spline);
                break;
            }
            case 4:
            {
                BSpline< 4, double > spline;
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                    destMultiArrayRange(bres), spline);
                break;
            }
            case 5:
            {
                BSpline< 5, double > spline;
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                    destMultiArrayRange(bres), spline);
                break;
            }
            default:
            {
                PyErr_SetString(PyExc_ValueError, "Spline order not supported.");
                python::throw_error_already_set();
            }
        }
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonResizeImageCatmullRomInterpolation(NumpyArray<3, Multiband<PixelType> > image,
                                                       python::object destSize,
                                                       NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 3) && (image.shape(1) > 3),
        "The input image must have a size of at least 4x4.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    vigra_precondition((size[0] > 1) && (size[1] > 1),
        "The destination image must have a size of at least 2x2.");
    // res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        // "Output image has wrong dimensions");
    TaggedShape newShape(image.taggedShape());
    newShape[0] = size[0];
    newShape[1] = size[1];
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty( newShape, "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);

        resizeImageCatmullRomInterpolation(srcImageRange(bimage),destImageRange(bres));
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonResizeImageCoscotInterpolation(NumpyArray<3, Multiband<PixelType> > image,
                                                   python::object destSize,
                                                   NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 3) && (image.shape(1) > 3),
        "The input image must have a size of at least 4x4.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    vigra_precondition((size[0] > 1) && (size[1] > 1),
        "The destination image must have a size of at least 2x2.");
    // res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        // "Output image has wrong dimensions");
    TaggedShape newShape(image.taggedShape());
    newShape[0] = size[0];
    newShape[1] = size[1];
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty( newShape, "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resizeImageCoscotInterpolation(srcImageRange(bimage),destImageRange(bres));
    }
    return res;
}

template <class PixelType>
NumpyAnyArray resamplingGaussian2D(NumpyArray<3, Multiband<PixelType> > image, 
    double sigmax, unsigned int derivativeOrderX, double samplingRatioX, double offsetX,
    double sigmay, unsigned int derivativeOrderY, double samplingRatioY, double offsetY, 
    NumpyArray<3, Multiband<PixelType> > res = python::object())
{
    vigra_precondition(samplingRatioX > 0 ,
       "resamplingGaussian(): samplingRatioX must be > 0.");
    vigra_precondition(samplingRatioY > 0 ,
       "resamplingGaussian(): samplingRatioY must be > 0.");
    Rational<int> xratio(samplingRatioX), yratio(samplingRatioY),
                  xoffset(offsetX), yoffset(offsetY);
    Gaussian< double > smoothx(sigmax, derivativeOrderX);
    Gaussian< double > smoothy(sigmay, derivativeOrderY);

	// res.reshapeIfEmpty(MultiArrayShape<3>::type(rational_cast< int >(image.shape(0)*xratio), 
                                                // rational_cast< int >(image.shape(1)*yratio), 
	                                            // image.shape(2)), 
                       // "resamplingGaussian2D(): Output array has wrong shape.");
    TaggedShape newShape(image.taggedShape());
    newShape[0] = rational_cast< int >(image.shape(0)*xratio);
    newShape[1] = rational_cast< int >(image.shape(1)*yratio);
    //FIXME: adjust axistag resolution
    res.reshapeIfEmpty( newShape, "resamplingGaussian2D(): Output array has wrong shape.");

	for(int k=0; k<image.shape(2); ++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    resamplingConvolveImage(srcImageRange(bimage), destImageRange(bres),
	            smoothx, xratio, xoffset, smoothy, yratio, yoffset);
	}
    return res;
}


/********************************************************************/
/*                                                                  */
/*                         SplineImageView                          */
/*                                                                  */
/********************************************************************/

template <class SplineView>
NumpyArray<2, Singleband<typename SplineView::value_type> >
SplineView_coefficientImage(SplineView const & self)
{
    NumpyArray<2, Singleband<typename SplineView::value_type> > res(self.shape());
    copyImage(srcImageRange(self.image()), destImage(res));
    return res;
}

template <class SplineView>
NumpyArray<2, Singleband<typename SplineView::value_type> >
SplineView_interpolatedImage(SplineView const & self, double xfactor, double yfactor, unsigned int xorder, unsigned int yorder)
{
    vigra_precondition(xfactor > 0.0 && yfactor > 0.0,
        "SplineImageView.interpolatedImage(xfactor, yfactor): factors must be positive.");
    int wn = int((self.width() - 1.0) * xfactor + 1.5);
    int hn = int((self.height() - 1.0) * yfactor + 1.5);
    NumpyArray<2, Singleband<typename SplineView::value_type> > res(MultiArrayShape<2>::type(wn, hn));
    for(int yn = 0; yn < hn; ++yn)
    {
        double yo = yn / yfactor;
        for(int xn = 0; xn < wn; ++xn)
        {
            double xo = xn / xfactor;
            res(xn, yn) = self(xo, yo, xorder, yorder);
        }
    }
    return res;
}

#define VIGRA_SPLINE_IMAGE(what, dx, dy) \
template <class SplineView> \
NumpyArray<2, Singleband<typename SplineView::value_type> > \
SplineView_##what##Image(SplineView const & self, double xfactor, double yfactor) \
{ \
    return SplineView_interpolatedImage(self, xfactor, yfactor, dx, dy); \
}

VIGRA_SPLINE_IMAGE(dx,  1, 0)
VIGRA_SPLINE_IMAGE(dxx, 2, 0)
VIGRA_SPLINE_IMAGE(dx3, 3, 0)
VIGRA_SPLINE_IMAGE(dy,  0, 1)
VIGRA_SPLINE_IMAGE(dyy, 0, 2)
VIGRA_SPLINE_IMAGE(dy3, 0, 3)
VIGRA_SPLINE_IMAGE(dxy, 1, 1)
VIGRA_SPLINE_IMAGE(dxxy, 2, 1)
VIGRA_SPLINE_IMAGE(dxyy, 1, 2)

#undef VIGRA_SPLINE_IMAGE

#define VIGRA_SPLINE_GRADIMAGE(what) \
template <class SplineView> \
NumpyArray<2, Singleband<typename SplineView::value_type> > \
SplineView_##what##Image(SplineView const & self, double xfactor, double yfactor) \
{ \
    vigra_precondition(xfactor > 0.0 && yfactor > 0.0, \
        "SplineImageView." #what "Image(xfactor, yfactor): factors must be positive."); \
    int wn = int((self.width() - 1.0) * xfactor + 1.5); \
    int hn = int((self.height() - 1.0) * yfactor + 1.5); \
    NumpyArray<2, Singleband<typename SplineView::value_type> > res(MultiArrayShape<2>::type(wn, hn)); \
    for(int yn = 0; yn < hn; ++yn) \
    { \
        double yo = yn / yfactor; \
        for(int xn = 0; xn < wn; ++xn) \
        { \
            double xo = xn / xfactor; \
            res(xn, yn) = self.what(xo, yo); \
        } \
    } \
    return res; \
}

VIGRA_SPLINE_GRADIMAGE(g2)
VIGRA_SPLINE_GRADIMAGE(g2x)
VIGRA_SPLINE_GRADIMAGE(g2y)

#undef VIGRA_SPLINE_GRADIMAGE

//FIXME: SplineImageView::coefficientArray() should be changed so that it can 
//       accept NumpyArray directly
template <class SplineView>
PyObject *
SplineView_facetCoefficients(SplineView const & self, double x, double y)
{
    BasicImage<double> coeff;
    self.coefficientArray(x, y, coeff);
    
    NumpyArray<2, double> res(MultiArrayShape<2>::type(coeff.width(), coeff.height()));
    copyImage(srcImageRange(coeff), destImage(res));
    
    python_ptr module(PyImport_ImportModule("numpy"), python_ptr::keep_count);
    pythonToCppException(module);
    python_ptr matrix(PyObject_GetAttrString(module, "matrix"), python_ptr::keep_count);
    pythonToCppException(matrix);

    return PyArray_View(res.pyArray(), 0, (PyTypeObject *)matrix.get());
}

template <class SplineView, class T>
SplineView *
pySplineView(NumpyArray<2, Singleband<T> > const & img)
{
    return new SplineView(srcImageRange(img), 0);
}

template <class SplineView, class T>
SplineView *
pySplineView1(NumpyArray<2, Singleband<T> > const & img, bool skipPrefilter)
{
    return new SplineView(srcImageRange(img), skipPrefilter);
}

template <class SplineView>
python::class_<SplineView> &
defSplineView(char const * name)
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);
    
    typedef typename SplineView::value_type Value;
    typedef typename SplineView::difference_type Shape;
    
    Value (SplineView::*callfct)(double, double) const = &SplineView::operator();
    Value (SplineView::*callfct2)(double, double, unsigned int, unsigned int) const = &SplineView::operator();

    static python::class_<SplineView> theclass(name, python::no_init);
    theclass
        .def("__init__", python::make_constructor(registerConverters(&pySplineView<SplineView, UInt8>)),
             "Construct a SplineImageView for the given image::\n\n"
             "    SplineImageView(image, skipPrefilter = False)\n\n"
             "Currently, 'image' can have dtype numpy.uint8, numpy.int32, and numpy.float32. "
             "If 'skipPrefilter' is True, image values are directly used as spline "
             "coefficients, so that the view performs approximation rather than interploation.\n\n")
        .def("__init__", python::make_constructor(registerConverters(&pySplineView<SplineView, npy_int32>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView<SplineView, float>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView1<SplineView, UInt8>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView1<SplineView, npy_int32>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView1<SplineView, float>)))
        .def("size", &SplineView::shape)
        .def("shape", &SplineView::shape, "The shape of the underlying image.\n\n")
        .def("width", &SplineView::width, "The width of the underlying image.\n\n")
        .def("height", &SplineView::height, "The height of the underlying image.\n\n")
        .def("isInside", &SplineView::isInside, 
             "Check if a coordinate is inside the underlying image.\n\n"
             "SplineImageView.isInside(x, y) -> bool\n\n")
        .def("isValid", &SplineView::isValid, 
             "Check if a coordinate is within the valid range of the SplineImageView.\n\n"
             "SplineImageView.isValid(x, y) -> bool\n\n"
             "Thanks to reflective boundary conditions, the valid range is three times "
             "as big as the size of the underlying image.\n\n")
        .def("__getitem__", (Value (SplineView::*)(Shape const &) const)&SplineView::operator(),
             "Return the value of the spline at a real-valued coordinate.\n\n"
             "Usage:\n\n"
             "    s = SplineImageView3(image)\n"
             "    value = s[10.1, 11.3]\n\n")
        .def("__call__", callfct,
             "Return the value of the spline or one of its derivatives at a real-valued coordinate.\n\n"
             "Usage:\n\n"
             "    s = SplineImageView3(image)\n"
             "    value = s(10.1, 11.3)\n"
             "    xorder = 1   # derivative order along x axis\n"
             "    yorder = 0   # derivative order along y axis\n"
             "    derivative = s(10.1, 11.3, xorder, yorder)\n\n")
        .def("__call__", callfct2)
        .def("dx", (Value (SplineView::*)(double, double) const)&SplineView::dx, args("x", "y"),
             "Return first derivative in x direction at a real-valued coordinate.\n\n"
             "SplineImageView.dx(x, y) -> value\n\n")
        .def("dy", (Value (SplineView::*)(double, double) const)&SplineView::dy, args("x", "y"),
             "Return first derivative in y direction at a real-valued coordinate.\n\n"
             "SplineImageView.dy(x, y) -> value\n\n")
        .def("dxx", (Value (SplineView::*)(double, double) const)&SplineView::dxx, args("x", "y"),
             "Return second derivative in x direction at a real-valued coordinate.\n\n"
             "SplineImageView.dxx(x, y) -> value\n\n")
        .def("dxy", (Value (SplineView::*)(double, double) const)&SplineView::dxy, args("x", "y"),
             "Return mixed second derivative at a real-valued coordinate.\n\n"
             "SplineImageView.dxy(x, y) -> value\n\n")
        .def("dyy", (Value (SplineView::*)(double, double) const)&SplineView::dyy, args("x", "y"),
             "Return second derivative in y direction at a real-valued coordinate.\n\n"
             "SplineImageView.dyy(x, y) -> value\n\n")
        .def("dx3", (Value (SplineView::*)(double, double) const)&SplineView::dx3, args("x", "y"),
             "Return third derivative in x direction at a real-valued coordinate.\n\n"
             "SplineImageView.dx3(x, y) -> value\n\n")
        .def("dxxy", (Value (SplineView::*)(double, double) const)&SplineView::dxxy, args("x", "y"),
             "Return mixed third derivative at a real-valued coordinate.\n\n"
             "SplineImageView.dxxy(x, y) -> value\n\n")
        .def("dxyy", (Value (SplineView::*)(double, double) const)&SplineView::dxyy, args("x", "y"),
             "Return mixed third derivative at a real-valued coordinate.\n\n"
             "SplineImageView.dxyy(x, y) -> value\n\n")
        .def("dy3", (Value (SplineView::*)(double, double) const)&SplineView::dy3, args("x", "y"),
             "Return third derivative in y direction at a real-valued coordinate.\n\n"
             "SplineImageView.dy3(x, y) -> value\n\n")
        .def("g2", (Value (SplineView::*)(double, double) const)&SplineView::g2, args("x", "y"),
             "Return gradient squared magnitude at a real-valued coordinate.\n\n"
             "SplineImageView.g2(x, y) -> value\n\n")
        .def("g2x", (Value (SplineView::*)(double, double) const)&SplineView::g2x, args("x", "y"),
             "Return first derivative in x direction of the gradient squared magnitude at a real-valued coordinate.\n\n"
             "SplineImageView.g2x(x, y) -> value\n\n")
        .def("g2y", (Value (SplineView::*)(double, double) const)&SplineView::g2y, args("x", "y"),
             "Return first derivative in y direction of the gradient squared magnitude at a real-valued coordinate.\n\n"
             "SplineImageView.g2y(x, y) -> value\n\n")
        .def("dxImage", &SplineView_dxImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dx`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dxImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dyImage", &SplineView_dyImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dy`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dyImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dxxImage", &SplineView_dxxImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dxx`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dxxImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dxyImage", &SplineView_dxyImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dxy`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dxyImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dyyImage", &SplineView_dyyImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dyy`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dyyImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dx3Image", &SplineView_dx3Image<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dx3`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dx3Image(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dxxyImage", &SplineView_dxxyImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dxxy`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dxxyImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dxyyImage", &SplineView_dxyyImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dxyy`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dxyyImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("dy3Image", &SplineView_dy3Image<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`dy3`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.dy3Image(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("g2Image", &SplineView_g2Image<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`g2`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.g2Image(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("g2xImage", &SplineView_g2xImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`g2x`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.g2xImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("g2yImage", &SplineView_g2yImage<SplineView>, (arg("xfactor") = 2.0, arg("yfactor") = 2.0),
             "Like :meth:`g2y`, but returns an entire image with the given sampling factors. For example,\n\n"
             "SplineImageView.g2yImage(2.0, 2.0) -> image\n\n"
             "creates an derivative image with two-fold oversampling in both directions.\n\n")
        .def("coefficientImage", &SplineView_coefficientImage<SplineView>)
        .def("interpolatedImage", &SplineView_interpolatedImage<SplineView>,
             (arg("xfactor") = 2.0, arg("yfactor") = 2.0, arg("xorder")=0, arg("yorder")=0),
             "Return an interpolated image or derivative image with the given sampling factors "
             "and derivative orders. For example, we get a two-fold oversampled image "
             "with the x-derivatives in each pixel by:\n\n"
             "SplineImageView.interpolatedImage(2.0, 2.0, 1, 0) -> image\n\n")
        .def("facetCoefficients", &SplineView_facetCoefficients<SplineView>,
             "SplineImageView.facetCoefficients(x, y) -> matrix\n\n"
             "Return the facet coefficient matrix so that spline values can be computed\n"
             "explicitly. The matrix has size (order+1)x(order+1), where order is \n"
             "the order of the spline. The matrix must be multiplied from left and right\n"
             "with the powers of the local facet x- and y-coordinates respectively\n"
             "(note that local facet coordinates are in the range [0,1] for odd order\n"
             "splines and [-0.5, 0.5] for even order splines).\n\n"
             "Usage for odd spline order:\n\n"
             "    s = SplineImageView3(image)\n"
             "    c = s.coefficients(10.1, 10.7)\n"
             "    x = matrix([1, 0.1, 0.1**2, 0.1**3])\n"
             "    y = matrix([1, 0.7, 0.7**2, 0.7**3])\n"
             "    assert abs(x * c * y.T - s[10.1, 10.7]) < smallNumber\n"
             "\n"
             "Usage for even spline order:\n\n"
             "    s = SplineImageView2(image)\n"
             "    c = s.coefficients(10.1, 10.7)\n"
             "    x = matrix([1, 0.1, 0.1**2])\n"
             "    y = matrix([1, -0.3, (-0.3)**2])\n"
             "    assert abs(x * c * y.T - s[10.1, 10.7]) < smallNumber\n\n")
        ;

    return theclass;
}


void defineSampling()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);

    enum_<RotationDirection>("RotationDirection")
        .value("CLOCKWISE",ROTATE_CW)
        .value("COUNTER_CLOCKWISE",ROTATE_CCW)
        .value("UPSIDE_DOWN",UPSIDE_DOWN);
    def("rotateImageRadiant",
        registerConverters(&pythonFreeRotateImageRadiant<float>),
        (arg("image"), arg("radiant"), arg("direction")=ROTATE_CW, arg("splineOrder")=0, arg("out")=object()),
        "Rotate an image by an arbitrary angle around its center using splines for interpolation.\n"
        "\n"
        "The angle may be given in radiant (parameter radiant).\n"
        "The parameter 'splineOrder' indicates the order of the splines used for interpolation.\n"
        "If the 'out' parameter is given, the image is cropped for it's dimensions. If the 'out'\n"
        "parameter is not given, an output image with the same dimensions as the input image is created.\n\n"
        "For more details, see GeometricTransformations.rotationMatrix2DRadians_ in the vigra C++ documentation.\n" 
        );
    def("rotateImageDegree",
        registerConverters(&pythonFreeRotateImageDegree<float>),
        (arg("image"), arg("degree"), arg("direction")=ROTATE_CW, arg("splineOrder")=0, arg("out")=object()),
        "Rotate an image by an arbitrary angle using splines for interpolation around its center.\n"
        "\n"
        "The angle may be given in degree (parameter degree).\n"
        "The parameter 'splineOrder' indicates the order of the splines used for interpolation.\n"
        "If the 'out' parameter is given, the image is cropped for it's dimensions. If the 'out'\n"
        "parameter is not given, an output image with the same dimensions as the input image is created.\n\n"
        "For more details, see GeometricTransformations.rotationMatrix2DDegrees_ in the vigra C++ documentation.\n"  
        );
    def("rotateImageSimple",
        registerConverters(&pythonFixedRotateImage<float>),
        (arg("image"), arg("orientation")=ROTATE_CW,arg("out")=object()),
        "Rotate an image by a multiple of 90 degrees.\n"
        "\n"
        "The 'orientation' parameter (which must be one of CLOCKWISE, COUNTER_CLOCKWISE and UPSIDE_DOWN\n"
        "indicates the rotation direction. The 'out' parameter must, if given, have the according dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see rotateImage_ in the vigra C++ documentation.\n" 
        );

//    def("rotateImageAboutCenter",
//       &DUMMY_FUNCTION,               // also multiband
//        (arg("image"), arg("degrees"), arg("center"), arg("splineOrder")));

    def("resampleImage",
        registerConverters(&pythonResampleImage<float>),               // also multiband
        (arg("image"), arg("factor"),arg("out")=object()),
        "Resample an image by the given 'factor'\n"
        "\n"
        "The 'out' parameter must have, if given, the according dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see resampleImage_ in the vigra C++ documentation.\n" 
        );

    def("resamplingGaussian", registerConverters(&resamplingGaussian2D<float>),
          (arg("image"), 
           arg("sigmaX")=1.0, arg("derivativeOrderX")=0, arg("samplingRatioX")=2.0, arg("offsetX")=0.0, 
           arg("sigmaY")=1.0, arg("derivativeOrderY")=0, arg("samplingRatioY")=2.0, arg("offsetY")=0.0, 
           arg("out") = python::object()),
          "Resample image using a gaussian filter::\n\n"
          "   resamplingGaussian(image,\n"
          "                      sigmaX=1.0, derivativeOrderX=0, samplingRatioX=2.0, offsetX=0.0,\n"
          "                      sigmaY=1.0, derivativeOrderY=0, samplingRatioY=2.0, offsetY=0.0,\n"
          "                      out=None)\n"
          "\n"
          "This function utilizes resamplingConvolveImage_ with a Gaussianfilter\n"
          "(see the vigra C++ documentation for details).\n\n");

    def("resizeImageNoInterpolation",
        registerConverters(&pythonResizeImageNoInterpolation<float>),               // also multiband
        (arg("image"), arg("shape")=object(), arg("out")=object()),
        "Resize image by repeating the nearest pixel values.\n"
        "\n"
        "The desired shape of the output image is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see resizeImageNoInterpolation_ in the vigra C++ documentation.\n"
        );


    def("resizeImageLinearInterpolation",
        registerConverters(&pythonResizeImageLinearInterpolation<float>),               // also multiband>
        (arg("image"), arg("shape")=object(), arg("out")=object()),
        "Resize image using linear interpolation.\n"
        "The function uses the standard separable bilinear interpolation algorithm to obtain a good compromise between quality and speed.\n"
        "\n" 
        "The desired shape of the output image is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see resizeImageLinearInterpolation_ in the vigra C++ documentation.\n"
        );

    def("resizeImageSplineInterpolation",
        registerConverters(&pythonResizeImageSplineInterpolation<float,3>),               // also multiband
        (arg("image"), arg("shape")=object(), arg("order") = 3, arg("out") = object()),
        "Resize image using B-spline interpolation.\n"
        "\n"
        "The spline order is given in the parameter 'order'.\n"
        "The desired shape of the output image is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see resizeImageSplineInterpolation_ in the vigra C++ documentation.\n"
       );

    def("resizeImageCatmullRomInterpolation",
        registerConverters(&pythonResizeImageCatmullRomInterpolation<float>),               // also multiband
        (arg("image"), arg("shape")=object(), arg("out")=object()),
        "Resize image using the Catmull/Rom interpolation function.\n"
        "\n" 
        "The desired shape of the output image is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see resizeImageCatmullRomInterpolation_ in the vigra C++ documentation.\n"
       );

    def("resizeImageCoscotInterpolation",
        registerConverters(&pythonResizeImageCoscotInterpolation<float>),               // also multiband
        (arg("image"), arg("shape")=object(), arg("out")=object()),
        "Resize image using the Coscot interpolation function.\n" 
        "\n" 
        "The desired shape of the output image is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree.\n"
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "For more details, see resizeImageCoscotInterpolation_ in the vigra C++ documentation.\n"
        );

    def("resizeVolumeSplineInterpolation",
        registerConverters(&pythonResizeImageSplineInterpolation<float,4>),               // also multiband
        (arg("image"), arg("shape")=object(), arg("order") = 3, arg("out") = object()),
        "Resize volume using B-spline interpolation.\n"
        "\n"
        "The spline order is given in the parameter 'order'.\n"
        "The dimensions of the output volume is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree.\n"
        "This function also works for multiband volumes, it is then executed on every band.\n\n"
        "For more details, see resizeMultiArraySplineInterpolation_ in the vigra C++ documentation.\n"
       );


    def("resize", registerConverters(&pythonResizeImageSplineInterpolation<float,3>),
        (arg("image"), arg("shape")=object(), arg("order") = 3, arg("out") = object()),
        "Resize image or volume using B-spline interpolation.\n"
        "\n"
        "The spline order is given in the parameter 'order'.\n"
        "The desired shape of the output array is taken either from 'shape' or 'out'.\n"
        "If both are given, they must agree. This function also works for multi-channel "
        "data, it is then executed on every channel independently.\n\n"
        "For more details, see resizeImageSplineInterpolation_ and "
        "resizeMultiArraySplineInterpolation_ in the vigra C++ documentation.\n"
       );

       def("resize",
        registerConverters(&pythonResizeImageSplineInterpolation<float,4>),
        (arg("image"), arg("shape")=object(), arg("order") = 3, arg("out") = object()));


    defSplineView<SplineImageView<0, float> >("SplineImageView0");
    defSplineView<SplineImageView<1, float> >("SplineImageView1");
    defSplineView<SplineImageView<2, float> >("SplineImageView2");
    defSplineView<SplineImageView<3, float> >("SplineImageView3");
    defSplineView<SplineImageView<4, float> >("SplineImageView4");
    defSplineView<SplineImageView<5, float> >("SplineImageView5");

}

} // namespace vigra


using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(sampling)
{
    import_vigranumpy();
    defineSampling();
}

