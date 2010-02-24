/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyconvolution_PyArray_API
//#define NO_IMPORT_ARRAY

#include <Python.h>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/functorexpression.hxx>
#include <vigra/convolution.hxx>
#include <vigra/nonlineardiffusion.hxx>
#include <vigra/resampling_convolution.hxx>
#include <vigra/recursiveconvolution.hxx>
#include <vigra/splineimageview.hxx>
#include "vigranumpykernel.hxx"

namespace python = boost::python;

namespace vigra
{

template < class PixelType>
NumpyAnyArray structureTensor2D(NumpyArray<3, Multiband<PixelType> > image, double innerScale, double outerScale,
                                NumpyArray<2, TinyVector<PixelType, 3> > res=python::object() )
{
    using namespace vigra::functor;
    
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "structureTensor2D(): Output array has wrong shape.");
    
    MultiArrayView<2, PixelType, StridedArrayTag> band = image.bindOuter(0);
    structureTensor(srcImageRange(band), 
                    destImage(res), 
                    innerScale, outerScale);
    
    if(image.shape(2) > 1)
    {
        MultiArray<2, TinyVector<PixelType, 3> > st(res.shape());
        for(int b=1; b<image.shape(2); ++b)
        {
            MultiArrayView<2, PixelType, StridedArrayTag> band = image.bindOuter(b);
            structureTensor(srcImageRange(band), 
                            destImage(st), 
                            innerScale, outerScale);
            combineTwoImages(srcImageRange(res), srcImage(st), destImage(res), Arg1() + Arg2());
        }
    }
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pystructureTensor2D, structureTensor2D)

template <class PixelType>
NumpyAnyArray nonlinearDiffusion2D(NumpyArray<3, Multiband<PixelType> > image, 
                                   double edgeThreshold, double scale,
                                   NumpyArray<3, Multiband<float> > res=python::object())
{
	res.reshapeIfEmpty(image.shape(), "nonlinearDiffusion2D(): Output array has wrong shape.");
    for(int k=0; k<image.shape(2); ++k)
    {
        MultiArrayView<2, float, StridedArrayTag> bres = res.bindOuter(k);
        nonlinearDiffusion(srcImageRange(image.bindOuter(k)), 
                        destImage(bres), 
                        DiffusivityFunctor< double >(edgeThreshold), scale);
    }
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pynonlinearDiffusion2D, nonlinearDiffusion2D)

template <class PixelType>
NumpyAnyArray simpleSharpening2D(NumpyArray<3, Multiband<PixelType> > image, double sharpeningFactor,
								 NumpyArray<3, Multiband<PixelType> > res=python::object() )
{
	res.reshapeIfEmpty(image.shape(), "simpleSharpening2D(): Output array has wrong shape.");
    vigra_precondition(sharpeningFactor >= 0 ,
       "simpleSharpening(): sharpeningFactor must be >= 0.");
    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        simpleSharpening(srcImageRange(bimage), destImage(bres),
                sharpeningFactor);
    }
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pysimpleSharpening2D, simpleSharpening2D)

template <class PixelType>
NumpyAnyArray gaussianSharpening2D(NumpyArray<3, Multiband<PixelType> > image,
    double sharpeningFactor, double scale, NumpyArray<3, Multiband<PixelType> > res=python::object() )
{
	res.reshapeIfEmpty(image.shape(), "gaussianSharpening2D(): Output array has wrong shape.");
	vigra_precondition(sharpeningFactor >= 0 ,
       "gaussianSharpening(): sharpeningFactor must be >= 0.");
    vigra_precondition(sharpeningFactor >= 0 ,
       "gaussianSharpening(): scale must be >= 0.");
    for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    gaussianSharpening(srcImageRange(bimage), destImage(bres),
	            sharpeningFactor, scale);
	}
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pygaussianSharpening2D, gaussianSharpening2D)

template <class PixelType>
NumpyAnyArray laplacianOfGaussian2D(NumpyArray<3, Multiband<PixelType> > image,
    double scale, NumpyArray<3, Multiband<PixelType> > res=python::object() )
{
	res.reshapeIfEmpty(image.shape(), "laplacianOfGaussian2D(): Output array has wrong shape.");
	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    laplacianOfGaussian(srcImageRange(bimage), destImage(bres), scale);
	}
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pylaplacianOfGaussian2D, laplacianOfGaussian2D)

template <class PixelType>
NumpyAnyArray hessianMatrixOfGaussian2D(NumpyArray<2, Singleband<PixelType> > image, double scale,
		NumpyArray<3, Multiband<PixelType> > res = python::object())
{
	res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(0), image.shape(1), 3), "hessianMatrixOfGaussian2D(): Output array has wrong shape.");
	MultiArrayView<2, PixelType, StridedArrayTag> hxx = res.bindOuter(0);
	MultiArrayView<2, PixelType, StridedArrayTag> hxy = res.bindOuter(1);
	MultiArrayView<2, PixelType, StridedArrayTag> hyy = res.bindOuter(2);
	hessianMatrixOfGaussian(srcImageRange(image), destImage(hxx),
	        destImage(hxy), destImage(hyy), scale);
	return res;
}

template <class PixelType>
NumpyAnyArray resamplingGaussian2D(NumpyArray<3, Multiband<PixelType> > image, 
    double sigma, unsigned int derivativeOrder,
    double samplingRatioX, double offsetX,
    double samplingRatioY, double offsetY, 
    NumpyArray<3, Multiband<PixelType> > res = python::object())
{
    vigra_precondition(samplingRatioX > 0 ,
       "resamplingGaussian(): samplingRatioX must be > 0.");
    vigra_precondition(samplingRatioY > 0 ,
       "resamplingGaussian(): samplingRatioY must be > 0.");
    Rational<int> xratio(samplingRatioX), yratio(samplingRatioY),
                  xoffset(offsetX), yoffset(offsetY);
    Gaussian< double > smooth(sigma, derivativeOrder);

	res.reshapeIfEmpty(MultiArrayShape<3>::type(rational_cast< int >(image.shape(0)*xratio), rational_cast< int >(image.shape(1)*yratio), 
	                                            image.shape(2)), "resamplingGaussian2D(): Output array has wrong shape.");

	for(int k=0; k<image.shape(2); ++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    resamplingConvolveImage(srcImageRange(bimage), destImageRange(bres),
	            smooth, xratio, xoffset, smooth, yratio, yoffset);
	}
    return res;
}

template <class PixelType>
NumpyAnyArray pythonConvolveImage(NumpyArray<3, Multiband<PixelType> > image,
                                  TwoDKernel const & kernel, 
                                  NumpyArray<3, Multiband<PixelType> > res = python::object())
{
	res.reshapeIfEmpty(image.shape(), "convolve2D(): Output array has wrong shape.");

	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    convolveImage(srcImageRange(bimage), destImage(bres),
	                  kernel2d(kernel));
	}
    return res;
}

template <class PixelType>
NumpyAnyArray pythonRecursiveFilter1(NumpyArray<3, Multiband<PixelType> > image,
                                     double b, BorderTreatmentMode borderTreatment, 
                                     NumpyArray<3, Multiband<PixelType> > res = python::object())
{
	res.reshapeIfEmpty(image.shape(), "recursiveFilter2D(): Output array has wrong shape.");

	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    recursiveFilterX(srcImageRange(bimage), destImage(bres), b, borderTreatment);
	    recursiveFilterY(srcImageRange(bres), destImage(bres), b, borderTreatment);
	}
    return res;
}

template <class PixelType>
NumpyAnyArray pythonRecursiveFilter2(NumpyArray<3, Multiband<PixelType> > image,
                                     double b1, double b2, 
                                     NumpyArray<3, Multiband<PixelType> > res = python::object())
{
	res.reshapeIfEmpty(image.shape(), "recursiveFilter2D(): Output array has wrong shape.");

	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    recursiveFilterX(srcImageRange(bimage), destImage(bres), b1, b2);
	    recursiveFilterY(srcImageRange(bres), destImage(bres), b1, b2);
	}
    return res;
}


template <class PixelType>
NumpyAnyArray pythonRecursiveSmooth(NumpyArray<3, Multiband<PixelType> > image,
                                    double scale, BorderTreatmentMode borderTreatment, 
                                    NumpyArray<3, Multiband<PixelType> > res = python::object())
{
    return pythonRecursiveFilter1(image, std::exp(-1.0/scale), borderTreatment, res);
}

template <class PixelType>
NumpyAnyArray pythonRecursiveGradient(NumpyArray<2, Singleband<PixelType> > image,
                                      double scale, 
                                      NumpyArray<2, TinyVector<PixelType, 2> > res = python::object())
{
	res.reshapeIfEmpty(image.shape(), "recursiveGradient2D(): Output array has wrong shape.");

    VectorComponentValueAccessor<TinyVector<PixelType, 2> > band(0);
    recursiveFirstDerivativeX(srcImageRange(image), destImage(res, band), scale);
    recursiveSmoothY(srcImageRange(res, band), destImage(res, band), scale);

    band.setIndex(1);
    recursiveSmoothX(srcImageRange(image), destImage(res, band), scale);
    recursiveFirstDerivativeY(srcImageRange(res, band), destImage(res, band), scale);
    
    return res;
}

template <class PixelType>
NumpyAnyArray pythonRecursiveLaplacian(NumpyArray<3, Multiband<PixelType> > image,
                                     double scale, 
                                     NumpyArray<3, Multiband<PixelType> > res = python::object())
{
	using namespace vigra::functor;
	
	res.reshapeIfEmpty(image.shape(), "recursiveLaplacian2D(): Output array has wrong shape.");

    MultiArrayShape<2>::type tmpShape(image.shape().begin());
    MultiArray<2, PixelType > tmp(tmpShape);
	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);

        recursiveSecondDerivativeX(srcImageRange(bimage), destImage(bres), scale);
        recursiveSmoothY(srcImageRange(bres), destImage(bres), scale);

        recursiveSmoothX(srcImageRange(bimage), destImage(tmp), scale);
        recursiveSecondDerivativeY(srcImageRange(tmp), destImage(tmp), scale);
        
        combineTwoImages(srcImageRange(bres), srcImage(tmp), destImage(bres), Arg1()+Arg2());
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
    typedef typename SplineView::value_type Value;
    typedef typename SplineView::difference_type Shape;
    
    Value (SplineView::*callfct)(double, double) const = &SplineView::operator();
    Value (SplineView::*callfct2)(double, double, unsigned int, unsigned int) const = &SplineView::operator();

    static python::class_<SplineView> theclass(name, python::no_init);
    theclass
        .def("__init__", python::make_constructor(registerConverters(&pySplineView<SplineView, UInt8>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView<SplineView, npy_int32>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView<SplineView, float>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView1<SplineView, UInt8>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView1<SplineView, npy_int32>)))
        .def("__init__", python::make_constructor(registerConverters(&pySplineView1<SplineView, float>)))
        .def("size", &SplineView::shape)
        .def("shape", &SplineView::shape)
        .def("width", &SplineView::width)
        .def("height", &SplineView::height)
        .def("isInside", &SplineView::isInside)
        .def("isValid", &SplineView::isValid)
        .def("__getitem__", (Value (SplineView::*)(Shape const &) const)&SplineView::operator())
        .def("__call__", callfct)
        .def("__call__", callfct2)
        .def("dx", (Value (SplineView::*)(double, double) const)&SplineView::dx)
        .def("dy", (Value (SplineView::*)(double, double) const)&SplineView::dy)
        .def("dxx", (Value (SplineView::*)(double, double) const)&SplineView::dxx)
        .def("dxy", (Value (SplineView::*)(double, double) const)&SplineView::dxy)
        .def("dyy", (Value (SplineView::*)(double, double) const)&SplineView::dyy)
        .def("dx3", (Value (SplineView::*)(double, double) const)&SplineView::dx3)
        .def("dxxy", (Value (SplineView::*)(double, double) const)&SplineView::dxxy)
        .def("dxyy", (Value (SplineView::*)(double, double) const)&SplineView::dxyy)
        .def("dy3", (Value (SplineView::*)(double, double) const)&SplineView::dy3)
        .def("g2", (Value (SplineView::*)(double, double) const)&SplineView::g2)
        .def("g2x", (Value (SplineView::*)(double, double) const)&SplineView::g2x)
        .def("g2y", (Value (SplineView::*)(double, double) const)&SplineView::g2y)
        .def("dxImage", &SplineView_dxImage<SplineView>)
        .def("dyImage", &SplineView_dyImage<SplineView>)
        .def("dxxImage", &SplineView_dxxImage<SplineView>)
        .def("dxyImage", &SplineView_dxyImage<SplineView>)
        .def("dyyImage", &SplineView_dyyImage<SplineView>)
        .def("dx3Image", &SplineView_dx3Image<SplineView>)
        .def("dxxyImage", &SplineView_dxxyImage<SplineView>)
        .def("dxyyImage", &SplineView_dxyyImage<SplineView>)
        .def("dy3Image", &SplineView_dy3Image<SplineView>)
        .def("g2Image", &SplineView_g2Image<SplineView>)
        .def("g2xImage", &SplineView_g2xImage<SplineView>)
        .def("g2yImage", &SplineView_g2yImage<SplineView>)
        .def("coefficientImage", &SplineView_coefficientImage<SplineView>)
        .def("interpolatedImage", &SplineView_interpolatedImage<SplineView>)
        .def("facetCoefficients", &SplineView_facetCoefficients<SplineView>,
             "SplineImageView::facetCoefficients(x, y)\n\n"
             "Return the facet coefficient matrix so that spline values can be computed\n"
             "explicitly. The matrix has size (order+1)x(order+1), where order is \n"
             "the order of the spline. The matrix must be multiplied from left and right\n"
             "with the powers of the local facet x- and y-coordinates respectively\n"
             "(note that local facet coordinates are in the range [0,1] for odd order\n"
             "splines and [-0.5, 0.5] for even order splines).\n\n"
             "Usage:\n\n"
             "s = SplineImageView3(image)  # odd order\n"
             "c = s.coefficients(10.1, 10.7)\n"
             "x = matrix([1, 0.1, 0.1**2, 0.1**3])\n"
             "y = matrix([1, 0.7, 0.7**2, 0.7**3])\n"
             "assert abs(x * c * y.T - s[10.1, 10.7]) < smallNumber\n"
             "\n"
             "s = SplineImageView2(image)  # even order\n"
             "c = s.coefficients(10.1, 10.7)\n"
             "x = matrix([1, 0.1, 0.1**2])\n"
             "y = matrix([1, -0.3, (-0.3)**2])\n"
             "assert abs(x * c * y.T - s[10.1, 10.7]) < smallNumber\n")
        ;

    return theclass;
}

void defineConvolutionFunctions()
{
    using namespace python;

    multidef("structureTensor", pystructureTensor2D<float>(),
      (arg("image"), arg("innerScale")=1.0, arg("outerScale")=2.0, arg("out") = python::object()),
      "Calculate the structure tensor for each pixel of an image, using Gaussian (derivative) filters at the given scales.\n"
      "\n"
      "For details see structureTensor_ in the vigra C++ documentation.");
    
    multidef("nonlinearDiffusion", pynonlinearDiffusion2D<float, UInt8>(),
        (arg("image"), arg("edgeThreshold"), arg("scale")),
        "Perform edge-preserving smoothing at the given scale."
        "\n\n"
        "For details see nonlinearDiffusion_ in the vigra C++ documentation.");

    multidef("simpleSharpening", pysimpleSharpening2D<float>(),
              (arg("image"), arg("sharpeningFactor")=1.0, arg("out") = python::object()),
              "Perform simple sharpening function.\n"
              "\n"
              "For details see simpleSharpening_ in the vigra C++ documentation.");

    multidef("gaussianSharpening", pygaussianSharpening2D<float>(),
                  (arg("image"), arg("sharpeningFactor")=1.0, arg("scale")=1.0,
                  arg("out") = python::object()),
                  "Perform sharpening function with gaussian filter."
                  "\n\n"
                  "For details see gaussianSharpening_ in the vigra C++ documentation.");

    multidef("laplacianOfGaussian", pylaplacianOfGaussian2D<float>(),
                          (arg("image"), arg("scale") = 1.0, arg("out") = python::object()),
                          "Filter image with the Laplacian of Gaussian operator at the given scale.\n"
                          "\n"
                          "For details see laplacianOfGaussian_ in the vigra C++ documentation.");

    def("hessianMatrixOfGaussian", registerConverters(&hessianMatrixOfGaussian2D<float>),
      (arg("image"), arg("scale")=1.0, arg("out") = python::object()),
      "Filter image with the 2nd derivatives of the Gaussian at the given scale to get the Hessian matrix.\n"
      "\n"
      "For details see hessianMatrixOfGaussian_ in the vigra C++ documentation.");

    def("resamplingGaussian", registerConverters(&resamplingGaussian2D<float>),
          (arg("image"), arg("sigma")=1.0, arg("derivativeOrder")=1,
           arg("samplingRatioX")=1, arg("samplingRatioY")=1, arg("out") = python::object()),
          "Resample image using a gaussian filter.\n"
          "\n"
          "This function utilize resamplingConvolveImage_ (see the vigra C++ documentation for details).");

    def("convolve2D", registerConverters(&pythonConvolveImage<float>),
              (arg("image"), arg("kernel")=python::object(), arg("out") = python::object()),
              "Perform 2D convolution with a 2D kernel (useful if the kernel is non-separable).\n"
              "\n"
              "For details see StandardConvolution.convolveImage_ in the vigra C++ documentation.");

    def("recursiveFilter2D", registerConverters(&pythonRecursiveFilter1<float>),
              (arg("image"), arg("b"), arg("borderTreament") = BORDER_TREATMENT_REFLECT, arg("out") = python::object()),
              "Perform 2D convolution with a first-order recursive filter with "
              "parameter 'b' and given 'borderTreatment'. 'b' must be between -1 and 1.\n"
              "\n"
              "For details see recursiveFilterX_ and recursiveFilterY_ (which "
              "this function calls in succession) in the vigra C++ documentation.");

    def("recursiveFilter2D", registerConverters(&pythonRecursiveFilter2<float>),
              (arg("image"), arg("b1"), arg("b2"), arg("out") = python::object()),
              "Perform 2D convolution with a second-order recursive filter with "
              "parameters 'b1' and 'b2'. Border treatment is always BORDER_TREATMENT_REFLECT.\n"
              "\n"
              "For details see recursiveFilterX_ and recursiveFilterY_ (which "
              "this function calls in succession) in the vigra C++ documentation.");

    def("recursiveSmooth2D", registerConverters(&pythonRecursiveSmooth<float>),
              (arg("image"), arg("scale"), arg("borderTreament") = BORDER_TREATMENT_REFLECT, arg("out") = python::object()),
              "Calls recursiveFilter2D() with b = exp(-1/scale), which "
              "corresponds to smoothing with an exponential filter exp(-abs(x)/scale).\n"
              "\n"
              "For details see recursiveSmoothLine_ in the vigra C++ documentation.");

    def("recursiveGradient2D", registerConverters(&pythonRecursiveSmooth<float>),
              (arg("image"), arg("scale"), arg("out") = python::object()),
              "Compute the gradient of a scalar image using a recursive (exponential) filter "
              "at the given 'scale'. The output image (if given) must have two channels.\n"
              "\n"
              "For details see recursiveSmoothLine_ and recursiveFirstDerivativeLine_ (which "
              "this function calls internally) in the vigra C++ documentation.");

    def("recursiveLaplacian2D", registerConverters(&pythonRecursiveLaplacian<float>),
              (arg("image"), arg("scale"), arg("out") = python::object()),
              "Compute the gradient of a 2D scalar or multiband image using a recursive (exponential) filter "
              "at the given 'scale'. The output image (if given) must have as many channels as the input.\n"
              "\n"
              "For details see recursiveSmoothLine_ and recursiveSecondDerivativeLine_ (which "
              "this function calls internally) in the vigra C++ documentation.");

    defSplineView<SplineImageView<0, float> >("SplineImageView0");
    defSplineView<SplineImageView<1, float> >("SplineImageView1");
    defSplineView<SplineImageView<2, float> >("SplineImageView2");
    defSplineView<SplineImageView<3, float> >("SplineImageView3");
    defSplineView<SplineImageView<4, float> >("SplineImageView4");
    defSplineView<SplineImageView<5, float> >("SplineImageView5");

}

void defineKernels();
void defineMultiConvolutionFunctions();

} // namespace vigra

using namespace vigra;

BOOST_PYTHON_MODULE_INIT(convolution)
{
	import_vigranumpy();
    defineKernels();
    defineConvolutionFunctions();
    defineMultiConvolutionFunctions();
}
