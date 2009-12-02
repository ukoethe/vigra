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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycmodule_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/convolution.hxx>
#include <vigra/nonlineardiffusion.hxx>
#include <vigra/resampling_convolution.hxx>
#include <vigra/recursiveconvolution.hxx>
#include <vigra/splineimageview.hxx>
#include "vigranumpykernel.hxx"

namespace python = boost::python;

namespace vigra
{

// TODO: Multiband??
template < class PixelType>
NumpyAnyArray structureTensor2D(NumpyArray<2, Singleband<PixelType> > image, double innerScale, double outerScale,
                                NumpyArray<2, TinyVector<float, 3> > res=python::object() )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "structureTensor2D(): Output array has wrong shape.");
    structureTensor(srcImageRange(image), 
                    destImage(res), 
                    innerScale, outerScale);
    return res;
}
//VIGRA_PYTHON_MULTITYPE_FUNCTOR(pystructureTensor2D, structureTensor2D)

// Multiband: addieren...

template <class PixelType>
NumpyAnyArray nonlinearDiffusion2D(NumpyArray<3, Multiband<PixelType> > image, double edgeThreshold, double scale)
{
    NumpyArray<3, Multiband<float> > res(image.shape());
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
NumpyAnyArray gaussianGradientMagnitude2D(NumpyArray<3, Multiband<PixelType> > image,
    double scale, NumpyArray<3, Multiband<PixelType> > res=python::object() )
{
	res.reshapeIfEmpty(image.shape(), "gaussianGradientMagnitude2D(): Output array has wrong shape.");
	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    gaussianGradientMagnitude(srcImageRange(bimage), destImage(bres),
	            scale);
	}
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pygaussianGradientMagnitude2D, gaussianGradientMagnitude2D)

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
		//NumpyArray<2, TinyVector<PixelType, 3> > res = python::object())
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
NumpyAnyArray resamplingGaussian2D(NumpyArray<3, Multiband<PixelType> > image, double sigma, unsigned int derivativeOrder,
    int samplingRatioX, int samplingRatioY, int offsetRatio, NumpyArray<3, Multiband<PixelType> > res = python::object())
{
    vigra_precondition(samplingRatioX > 0 ,
       "resamplingGaussian(): samplingRatioX must be > 0.");
    vigra_precondition(samplingRatioY > 0 ,
       "resamplingGaussian(): samplingRatioY must be > 0.");
    Rational<int> xratio(samplingRatioX), yratio(samplingRatioY),
        offset(offsetRatio);
    Gaussian< double > smooth(sigma, derivativeOrder);

	res.reshapeIfEmpty(MultiArrayShape<3>::type(rational_cast< int >(image.shape(0)*xratio), rational_cast< int >(image.shape(1)*yratio), image.shape(2)), "resamplingGaussian2D(): Output array has wrong shape.");

	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    resamplingConvolveImage(srcImageRange(bimage), destImageRange(bres),
	            smooth, xratio, offset, smooth, yratio, offset);
	}
    return res;
}

template <class PixelType>
NumpyAnyArray pythonConvolveImage(NumpyArray<3, Multiband<PixelType> > image,
    TwoDKernel & kernel, NumpyArray<3, Multiband<PixelType> > res = python::object())
{
	res.reshapeIfEmpty(image.shape(), "pythonConvolveImage(): Output array has wrong shape.");

	for(int k=0;k<image.shape(2);++k)
	{
	    MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
	    MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
	    convolveImage(srcImageRange(bimage), destImage(bres),
	            kernel2d(kernel));
	}
    return res;
}

#if 0
template <class PixelType>
NumpyAnyArray pythonRecursiveFilterLineFirstOrder(NumpyArray<2, Singleband<PixelType> > image,
    double b1, BorderTreatmentMode border, NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "recursiveFilterLineFirstOrder(): Output array has wrong shape.");

    ConstStridedImageIterator<PixelType>
        src_beg(image.data(), 1, image.stride(0), image.stride(1)),
        src_end = src_beg + Size2D(image.shape(0), image.shape(1));
    MultibandVectorAccessor<PixelType> src_acc(image.shape(2), image.stride(2));

    ConstStridedImageIterator<PixelType>
            dest_beg(res.data(), 1, res.stride(0), res.stride(1)),
            dest_end = dest_beg + Size2D(res.shape(0), res.shape(1));

    /*recursiveFilterLine(image.begin(), image.end(), image.accessor(),
        res.begin(), res.accessor(), b1, border);*/

    recursiveFilterLine(src_beg, src_end, src_acc,
            dest_beg, dest_end, b1, border);
    return res;
}
#endif

#if 0
template < class PixelType, int pyArrayTypeConstant >
PyObject* nonlinearDiffusion2D(BasicImageView<
    PixelType > const & image, double edgeThreshold, double scale)
{
    vigra_precondition(scale >= 0 ,
       "nonlinearDiffusion(): scale must be >= 0.");
    vigra_precondition(edgeThreshold >= 0 ,
       "nonlinearDiffusion(): edgeThreshold must be >= 0.");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    nonlinearDiffusion(srcImageRange(image), destImage(*imageView),
        DiffusivityFunctor< double >(edgeThreshold), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* simpleSharpening2D(BasicImageView< PixelType > const & image,
    double sharpeningFactor)
{
    vigra_precondition(sharpeningFactor >= 0 ,
       "simpleSharpening(): sharpeningFactor must be >= 0.");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    simpleSharpening(srcImageRange(image), destImage(*imageView),
        sharpeningFactor);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* gaussianSharpening2D(BasicImageView< PixelType > const & image,
    double sharpeningFactor, double scale)
{
    vigra_precondition(sharpeningFactor >= 0 ,
       "gaussianSharpening(): sharpeningFactor must be >= 0.");
    vigra_precondition(sharpeningFactor >= 0 ,
       "gaussianSharpening(): scale must be >= 0.");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    gaussianSharpening(srcImageRange(image), destImage(*imageView),
        sharpeningFactor, scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class SrcPixelType, class DestPixelType, int pyArrayTypeConstant >
PyObject* gaussianGradientMagnitude2D(BasicImageView< SrcPixelType >
    const & image, double scale)
{
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< DestPixelType >::isScalar());
    BasicImageView< DestPixelType >* imageView = createView2D<
        DestPixelType >(array,
        typename NumericTraits< DestPixelType >::isScalar());
    gaussianGradientMagnitude(srcImageRange(image), destImage(*imageView),
        scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* laplacianOfGaussian2D(BasicImageView< PixelType > const & image,
    double scale)
{
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    laplacianOfGaussian(srcImageRange(image), destImage(*imageView), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* hessianMatrixOfGaussian2D(BasicImageView<
    PixelType > const & image, double scale)
{
    PyObject* arrayxx = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageViewxx = createView2D<
        PixelType >(arrayxx, typename NumericTraits< PixelType >::isScalar());
    PyObject* arrayxy = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageViewxy = createView2D<
        PixelType >(arrayxy, typename NumericTraits< PixelType >::isScalar());
    PyObject* arrayyy = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageViewyy = createView2D<
        PixelType >(arrayyy, typename NumericTraits< PixelType >::isScalar());

    hessianMatrixOfGaussian(srcImageRange(image), destImage(*imageViewxx),
        destImage(*imageViewxy), destImage(*imageViewyy), scale);
    delete imageViewxx;
    delete imageViewxy;
    delete imageViewyy;
    return PyTuple_Pack(3, PyArray_Return((PyArrayObject*) arrayxx),
        PyArray_Return((PyArrayObject*) arrayxy),
        PyArray_Return((PyArrayObject*) arrayyy));
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* resamplingGaussian2D(BasicImageView< PixelType >
    const & image, double sigma, unsigned int derivativeOrder,
    int samplingRatioX, int samplingRatioY, int offsetRatio)
{
    vigra_precondition(samplingRatioX > 0 ,
       "resamplingGaussian(): samplingRatioX must be > 0.");
    vigra_precondition(samplingRatioY > 0 ,
       "resamplingGaussian(): samplingRatioY must be > 0.");
    Rational< int > xratio(samplingRatioX), yratio(samplingRatioY),
        offset(offsetRatio);
    Gaussian< double > smooth(sigma, derivativeOrder);
    PyObject* array = createNumpyArray(Size2D(
        rational_cast< int >(xratio*image.width()),
        rational_cast< int >(yratio*image.height())), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    resamplingConvolveImage(srcImageRange(image), destImageRange(*imageView),
        smooth, xratio, offset, smooth, yratio, offset);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

void py2DKernel_initExplicitly(TwoDKernel & self, int upperleftX,
    int upperleftY, int lowerrightX, int lowerrightY,
    python::object const & args)
{
    vigra_precondition(upperleftX <= 0 ,
       "initExplicitly(): upperleftX must be <= 0.");
    vigra_precondition(upperleftY <= 0 ,
       "initExplicitly(): upperleftY must be <= 0.");
    vigra_precondition(lowerrightX >= 0 ,
       "initExplicitly(): lowerrightX must be >= 0.");
    vigra_precondition(lowerrightY >= 0 ,
       "initExplicitly(): lowerrightY must be >= 0.");
    Diff2D upperleft(upperleftX, upperleftY);
    Diff2D lowerright(lowerrightX, lowerrightY);
    if(! PySequence_Check(args.ptr()))
    {
        KernelValueType value = python::extract<KernelValueType>(args);
        self.initExplicitly(upperleft, lowerright) = value;
    }
    else
    {
        KernelValueType value = python::extract<KernelValueType>(args[0]);
        TwoDKernel::InitProxy ip = self.initExplicitly(upperleft, lowerright) =
            value;
        if(python::len(args) != (self.width() * self.height()))
        {
            std::stringstream str;
            str << "Wrong number of init values. The number must be ";
            str << self.width() * self.height();
            PyErr_SetString(PyExc_ValueError, str.str().c_str());
            python::throw_error_already_set();
        }
        else
        {
            int size = self.width() * self.height();
            for(int i=1; i<size; ++i)
            {
                ip,(python::extract<KernelValueType>(args[i]));
            }
        }
    }
}

void py2DKernel_copy(TwoDKernel & self, TwoDKernel const & kernel)
{
    self=kernel;
}

KernelValueType py2DKernel_getitem(TwoDKernel const & self,
    python::tuple const & pos)
{
    if(python::len(pos) != 2)
    {
        std::stringstream str;
        str << "Wrong number of positional values. The number must be 2.";
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
    }
    int x = python::extract<int>(pos[0]);
    int y = python::extract<int>(pos[1]);
    Diff2D position(x, y);
    Point2D ul(self.upperLeft());
    Point2D lr(self.lowerRight());
    if((ul.px() <= position.x) &&
        (ul.py() <= position.y) &&
        (lr.px() >= position.x) &&
        (lr.py() >= position.y))
    {
        return self[position];
    }
    else
    {
        std::stringstream str;
        str << "Bad position.";
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
        return 0;
    }
}

void py2DKernel_setitem(TwoDKernel & self, python::tuple const & pos,
    KernelValueType value)
{
    if(python::len(pos) != 2)
    {
        std::stringstream str;
        str << "Wrong number of positional values. The number must be 2.";
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
    }
    int x = python::extract<int>(pos[0]);
    int y = python::extract<int>(pos[1]);
    Diff2D position(x, y);
    Point2D ul(self.upperLeft());
    Point2D lr(self.lowerRight());
    if((ul.px() <= position.x) &&
        (ul.py() <= position.y) &&
        (lr.px() >= position.x) &&
        (lr.py() >= position.y))
    {
        self[position] = value;
    }
    else
    {
        std::stringstream str;
        str << "Bad position.";
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
    }
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* pythonConvolveImage(BasicImageView< PixelType > const & image,
    TwoDKernel & kernel)
{
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    convolveImage(srcImageRange(image), destImage(*imageView),
        kernel2d(kernel));
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterLineFirstOrder(BasicImageView< PixelType >
    const & image, double b1, BorderTreatmentMode border)
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFilterLine(image.begin(), image.end(), image.accessor(),
        imageView->begin(), imageView->accessor(), b1, border);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterLineSecondOrder(BasicImageView< PixelType >
    const & image, double b1, double b2)
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    vigra_precondition((-1 < b2) && (b2 < 1),
       "-1 < b2 < 1");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFilterLine(image.begin(), image.end(), image.accessor(),
        imageView->begin(), imageView->accessor(), b1, b2);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterSmoothLine2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveSmoothLine(image.begin(), image.end(), image.accessor(),
        imageView->begin(), imageView->accessor(), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFirstDerivativeLine2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFirstDerivativeLine(image.begin(), image.end(), image.accessor(),
        imageView->begin(), imageView->accessor(), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveSecondDerivativeLine2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveSecondDerivativeLine(image.begin(), image.end(), image.accessor(),
        imageView->begin(), imageView->accessor(), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterXFirstOrder(BasicImageView< PixelType >
    const & image, double b1, BorderTreatmentMode border)
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFilterX(srcImageRange(image), destImage(*imageView), b1, border);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterXSecondOrder(BasicImageView< PixelType >
    const & image, double b1, double b2)
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    vigra_precondition((-1 < b2) && (b2 < 1),
       "-1 < b2 < 1");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFilterX(srcImageRange(image), destImage(*imageView), b1, b2);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterYFirstOrder(BasicImageView< PixelType >
    const & image, double b1, BorderTreatmentMode border)
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFilterY(srcImageRange(image), destImage(*imageView), b1, border);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFilterYSecondOrder(BasicImageView< PixelType >
    const & image, double b1, double b2)
{
    vigra_precondition((-1 < b1) && (b1 < 1),
       "-1 < b1 < 1");
    vigra_precondition((-1 < b2) && (b2 < 1),
       "-1 < b2 < 1");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFilterY(srcImageRange(image), destImage(*imageView), b1, b2);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveSmoothX2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveSmoothX(srcImageRange(image), destImage(*imageView), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveSmoothY2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveSmoothY(srcImageRange(image), destImage(*imageView), scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFirstDerivativeX2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFirstDerivativeX(srcImageRange(image), destImage(*imageView),
        scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveSecondDerivativeX2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveSecondDerivativeX(srcImageRange(image), destImage(*imageView),
        scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveFirstDerivativeY2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveFirstDerivativeY(srcImageRange(image), destImage(*imageView),
        scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class PixelType, int pyArrayTypeConstant >
PyObject* recursiveSecondDerivativeY2D(BasicImageView< PixelType >
    const & image, double scale)
{
    vigra_precondition(scale > 0, "scale > 0");
    PyObject* array = createNumpyArray(image.size(), pyArrayTypeConstant,
        typename NumericTraits< PixelType >::isScalar());
    BasicImageView< PixelType >* imageView = createView2D<
        PixelType >(array, typename NumericTraits< PixelType >::isScalar());
    recursiveSecondDerivativeY(srcImageRange(image), destImage(*imageView),
        scale);
    delete imageView;
    return PyArray_Return((PyArrayObject*) array);
}
#endif

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

namespace { void DUMMY_FUNCTION(int, int, int, int, int, int) {} }

void defineConvolutionFunctions()
{
    using namespace python;
    
    enum_<BorderTreatmentMode>("BorderTreatmentMode")
        .value("BORDER_TREATMENT_CLIP", BORDER_TREATMENT_CLIP)
        .value("BORDER_TREATMENT_AVOID", BORDER_TREATMENT_AVOID)
        .value("BORDER_TREATMENT_REFLECT", BORDER_TREATMENT_REFLECT)
        .value("BORDER_TREATMENT_REPEAT", BORDER_TREATMENT_REPEAT)
        .value("BORDER_TREATMENT_WRAP", BORDER_TREATMENT_WRAP)
    ;

    def("structureTensor", registerConverters(&structureTensor2D<float>),
      (arg("image"), arg("innerScale")=1.0, arg("outerScale")=2.0, arg("out") = python::object()));
    
    multidef("nonlinearDiffusion", pynonlinearDiffusion2D<float, UInt8>(),
        (arg("image"), arg("edgeThreshold"), arg("scale")));

    multidef("simpleSharpening", pysimpleSharpening2D<float>(),
              (arg("image"), arg("sharpeningFactor")=1.0, arg("out") = python::object()));

    multidef("gaussianSharpening", pygaussianSharpening2D<float>(),
                  (arg("image"), arg("sharpeningFactor")=1.0, arg("scale")=1.0,
                  arg("out") = python::object()));

    multidef("gaussianGradientMagnitude", pygaussianGradientMagnitude2D<float>(),
                      (arg("image"), arg("scale") = 1.0, arg("out") = python::object()));

    multidef("laplacianOfGaussian", pylaplacianOfGaussian2D<float>(),
                          (arg("image"), arg("scale") = 1.0, arg("out") = python::object()));

    def("hessianMatrixOfGaussian", registerConverters(&hessianMatrixOfGaussian2D<float>),
      (arg("image"), arg("scale")=1.0, arg("out") = python::object()));

    def("resamplingGaussian", registerConverters(&resamplingGaussian2D<float>),
          (arg("image"), arg("sigma")=1.0, arg("derivativeOrder")=1,
           arg("samplingRatioX")=1, arg("samplingRatioY")=1, arg("out") = python::object()));

    def("convolveImage", registerConverters(&pythonConvolveImage<float>),
              (arg("image"), arg("kernel")=python::object(), arg("out") = python::object()));

    /*def("recursiveFilterLineFirstOrder", registerConverters(&pythonRecursiveFilterLineFirstOrder<float>),
                  (arg("image"), arg("b1")=1.0, arg("borderTreatment") = BORDER_TREATMENT_REPEAT,
                   arg("out") = python::object()));*/

#if 0
    class_< Kernel >("Kernel1D",
        "Generic 1 dimensional convolution kernel.\n\n"
        "This kernel may be used for convolution of 1 dimensional signals or "
        "for separable convolution of multidimensional signals.\n\n"
        "The kernel's size is given by its left() (<=0) and right() (>= 0) "
        "methods. The desired border treatment mode is returned by "
        "getBorderTreatment().\n\n"
        "The different init functions create a kernel with the specified "
        "properties.\n\n",
        init<>("The standard constructor.")
        )
        .def(init< Kernel >(args("kernel"),
            "The copy constructor."))
        .def("initGaussian",
            (void (Kernel::*)(double))&Kernel::initGaussian,
            (arg("std_dev")),
            "Init as a Gaussian function with norm 1.")
        .def("initGaussian",
            (void (Kernel::*)(double, KernelValueType))&Kernel::initGaussian,
            (arg("std_dev"), arg("norm")),
            "Init as a sampled Gaussian function. The radius of the kernel is "
            "always 3*std_dev. 'norm' denotes the sum of all bins of the "
            "kernel (i.e. the kernel is corrected for the normalization error "
            "introduced by windowing the Gaussian to a finite interval). "
            "However, if norm is 0.0, the kernel is normalized to 1 by the "
            "analytic expression for the Gaussian, and no correction for the "
            "windowing error is performed.")
        .def("initDiscreteGaussian",
            (void (Kernel::*)(double))&Kernel::initDiscreteGaussian,
            (arg("std_dev")),
            "Init as a Lindeberg's discrete analog of the Gaussian function "
            "with norm 1.")
        .def("initDiscreteGaussian",
            (void (Kernel::*)(double,
            KernelValueType))&Kernel::initDiscreteGaussian,
            (arg("std_dev"), arg("norm")),
            "Init as Lindeberg's discrete analog of the Gaussian function. "
            "The radius of the kernel is always 3*std_dev. 'norm' denotes "
            "the sum of all bins of the kernel.")
        .def("initGaussianDerivative",
            (void (Kernel::*)(double, int))&Kernel::initGaussianDerivative,
            (arg("std_dev"), arg("order")),
            "Init as a Gaussian derivative with norm 1.")
        .def("initGaussianDerivative",
            (void (Kernel::*)(double, int,
            KernelValueType))&Kernel::initGaussianDerivative,
            (arg("std_dev"), arg("order"), arg("norm")),
            "Init as a Gaussian derivative of order 'order'. The radius of "
            "the kernel is always 3*std_dev + 0.5*order. 'norm' denotes "
            "the norm of the kernel. Thus, the kernel will be corrected for "
            "the error introduced by windowing the Gaussian to a finite "
            "interval. However, if norm is 0.0, the kernel is normalized to 1 "
            "by the analytic expression for the Gaussian derivative, and no "
            "correction for the windowing error is performed.")
        .def("initBinomial",
            (void (Kernel::*)(int))&Kernel::initBinomial,
            (arg("radius")),
            "Init as a Binomial filter with norm 1.")
        .def("initBinomial",
            (void (Kernel::*)(int, KernelValueType))&Kernel::initBinomial,
            (arg("radius"), arg("norm")),
            "Init as a Binomial filter. 'norm' denotes the sum of all bins "
            "of the kernel.")
        .def("initAveraging",
            (void (Kernel::*)(int))&Kernel::initAveraging,
            (arg("radius")),
            "Init as a Averaging filter with norm 1.")
        .def("initAveraging",
            (void (Kernel::*)(int, KernelValueType))&Kernel::initAveraging,
            (arg("radius"), arg("norm")),
            "Init as an Averaging filter. 'norm' denotes the sum of all "
            "bins of the kernel. The window size is (2*radius+1) * "
            "(2*radius+1)")
        .def("initSymmetricGradient",
            (void (Kernel::*)())&Kernel::initSymmetricGradient,
            "Init as a symmetric gradient filter with norm 1.")
        .def("initSymmetricGradient",
            (void (Kernel::*)(KernelValueType))&Kernel::initSymmetricGradient,
            (arg("norm")),
            "Init as a symmetric gradient filter of the form [ 0.5 * norm, "
            "0.0 * norm, -0.5 * norm]")
        .def("left", &Kernel::left,
            "Left border of kernel.")
        .def("right", &Kernel::right,
            "Right border of kernel.")
        .def("size", &Kernel::size,
            "Size of kernel.")
        .def("norm", &Kernel::norm,
            "Norm of kernel.")
        .def("normalize", (void (Kernel::*)())&Kernel::normalize,
            "Normalize kernel to norm 1.")
        .def("normalize", (void (Kernel::*)(KernelValueType, unsigned int,
            double))&Kernel::normalize,
            (arg("norm"), arg("derivativeOrder") = 0, arg("offset") = 0.0),
            "Set a new norm and normalize kernel, use the normalization "
            "formula for the given derivativeOrder.")
        .def("initExplicitly", &pyKernel_initExplicitly,
            (arg("left"), arg("right"), arg("initializersSequence")),
            "Init the kernel by an explicit initializer list. The left and "
            "right boundaries of the kernel must be passed. An initializer "
            "tuple is passed. The norm is set to the sum of the initialzer "
            "values. If the wrong number of values is given, a run-time error "
            "results. It is, however, possible to give just one initializer. "
            "This creates an averaging filter with the given constant.")
        .def("copy", pyKernel_copy,
            (arg("kernel1D")),
            "Copy parameter kernel \"kernel1D\" into the current kernel.")
        .def("__getitem__", pyKernel_getitem,
            (arg("position")),
            "Access kernel value at the specified location.")
        .def("__setitem__", pyKernel_setitem,
            (arg("position"), arg("value")),
            "Assign value at the specified location.")
        .def("setBorderTreatment", &Kernel::setBorderTreatment,
            (arg("borderTreatmentMode")),
            "Set border treatment mode.")
        .def("getBorderTreatment", &Kernel::borderTreatment,
            "Current border treatment mode.")
    ;

    class_< TwoDKernel >("Kernel2D",
        "Generic 2 dimensional convolution kernel.\n\n"
        "This kernel may be used for convolution of 2 dimensional signals.\n\n"
        "The desired border treatment mode is returned by borderTreatment()."
        "(Note that the 2D convolution functions don't currently support all "
        "modes.)\n\n"
        "The different init functions create a kernel with the specified "
        "properties.\n\n",
        init<>("The standard constructor.")
        )
        .def(init< TwoDKernel >(args("kernel"),
            "The copy constructor."))
        .def("initSeparable",
            (void (TwoDKernel::*)(Kernel &, Kernel &))
            &TwoDKernel::initSeparable,
            (arg("kernelX"), arg("kernelY")),
            "Init the 2D kernel as the cartesian product of two 1D kernels of "
            "type Kernel1D. The norm becomes the product of the two original "
            "norms.")
        .def("initDisk",
            (void (TwoDKernel::*)(int))&TwoDKernel::initDisk,
            (arg("radius")),
            "Init the 2D kernel as a circular averaging filter. The norm will "
            "be calculated as 1 / (number of non-zero kernel values).\n\n"
            "Precondition:\n\n"
            "   radius > 0")
        .def("width", &TwoDKernel::width,
            "Width of the kernel.")
        .def("height", &TwoDKernel::height,
            "Height of the kernel.")
        .def("norm", &TwoDKernel::norm,
            "Norm of kernel (i.e. sum of its elements).")
        .def("normalize", (void (TwoDKernel::*)())&TwoDKernel::normalize,
            "Normalize kernel to norm 1.")
        .def("normalize",
            (void (TwoDKernel::*)(KernelValueType))&TwoDKernel::normalize,
            (arg("norm")),
            "Normalize the kernel to the given value. (The norm is the sum of "
            "all kernel elements.)")
        .def("initExplicitly", &py2DKernel_initExplicitly,
            (arg("upperleftX"), arg("upperleftY"), arg("lowerrightY"),
            arg("lowerrightY"), arg("initializersSequence")),
            "Init the kernel by an explicit initializer list. The left and "
            "right boundaries of the kernel must be passed. An initializer "
            "tuple is passed. The norm is set to the sum of the initialzer "
            "values. If the wrong number of values is given, a run-time error "
            "results. It is, however, possible to give just one initializer. "
            "This creates an averaging filter with the given constant.")
        .def("copy", py2DKernel_copy,
            (arg("kernel2D")),
            "Copy parameter kernel \"kernel2D\" into the current kernel.")
        .def("__getitem__", py2DKernel_getitem,
            (arg("positionX"), arg("positionY")),
            "Access kernel value at the specified location.")
        .def("__setitem__", py2DKernel_setitem,
            (arg("positionX"), arg("positionY"), arg("value")),
            "Assign value at the specified location.")
        .def("setBorderTreatment", &TwoDKernel::setBorderTreatment,
            (arg("borderTreatmentMode")),
            "Set border treatment mode.")
        .def("getBorderTreatment", &TwoDKernel::borderTreatment,
            "Current border treatment mode.")
    ;
#endif


#if 0
    def("convolveOneDimension2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("dim"), arg("kernel1D")));

    def("separableConvolve2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("kernel"), arg("kernel2") = python::object()));

    def("gaussianSmoothing2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("sigma")));

    def("simpleSharpening2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("sharpeningFactor")));

    def("gaussianSharpening2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("sharpeningFactor"), arg("scale")));

    def("symmetricGradient2D",
        &DUMMY_FUNCTION,               // also multiband?
        (arg("image")));

    def("gaussianGradient2D",
        &DUMMY_FUNCTION,               // also multiband?
        (arg("image"), arg("sigma")));

    def("gaussianGradientMagnitude2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("laplacianOfGaussian2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("hessianMatrixOfGaussian2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("resamplingGaussian2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("sigma"), arg("derivativeOrder"),
        arg("samplingRatioX"), arg("samplingRatioY"), arg("offset")));


    def("convolveImage",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("kernel2D")));

    def("recursiveFilterLineFirstOrder",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("b1"), arg("borderTreatment")));

    def("recursiveFilterLineSecondOrder",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("b1"), arg("b2")));

    def("recursiveFilterSmoothLine2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveFirstDerivativeLine2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveSecondDerivativeLine2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveFilterXFirstOrder",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("b1"), arg("borderTreatment")));

    def("recursiveFilterXSecondOrder",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("b1"), arg("b2")));

    def("recursiveFilterYFirstOrder",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("b1"), arg("borderTreatment")));

    def("recursiveFilterYSecondOrder",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("b1"), arg("b2")));

    def("recursiveSmoothX2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveSmoothY2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveFirstDerivativeX2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveSecondDerivativeX2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveFirstDerivativeY2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));

    def("recursiveSecondDerivativeY2D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("image"), arg("scale")));
#endif

    defSplineView<SplineImageView<0, float> >("SplineImageView0");
    defSplineView<SplineImageView<1, float> >("SplineImageView1");
    defSplineView<SplineImageView<2, float> >("SplineImageView2");
    defSplineView<SplineImageView<3, float> >("SplineImageView3");
    defSplineView<SplineImageView<4, float> >("SplineImageView4");
    defSplineView<SplineImageView<5, float> >("SplineImageView5");

}

} // namespace vigra

