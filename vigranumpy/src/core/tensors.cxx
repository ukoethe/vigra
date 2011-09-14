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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyfilters_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/boundarytensor.hxx>
#include <vigra/orientedtensorfilters.hxx>
#include <vigra/tensorutilities.hxx>
#include <vigra/multi_tensorutilities.hxx>
#include "vigranumpyscaleparam.hxx"

namespace python = boost::python;

namespace vigra
{

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientND(NumpyArray<ndim, Singleband<VoxelType> > array,
                         python::object sigma,
                         NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res = NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> >(),
                         python::object sigma_d = python::object(0.0), 
                         python::object step_size = python::object(1.0),
                         double window_size = 0.0, 
                         python::object roi = python::object())
{
    pythonScaleParam<ndim> params(sigma, sigma_d, step_size, "gaussianGradient");
    params.permuteLikewise(array);
    std::string description("Gaussian gradient, scale=");
    description += asString(sigma);
    
    ConvolutionOptions<ndim> opt(params().filterWindowSize(window_size));
    
    if(roi != python::object())
    {
        typedef typename MultiArrayShape<ndim>::type Shape;
        Shape start = array.permuteLikewise(python::extract<Shape>(roi[0])());
        Shape stop  = array.permuteLikewise(python::extract<Shape>(roi[1])());
        opt.subarray(start, stop);
        res.reshapeIfEmpty(array.taggedShape().resize(stop-start).setChannelDescription(description), 
                       "gaussianGradient(): Output array has wrong shape.");
    }
    else
    {
        res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
                       "gaussianGradient(): Output array has wrong shape.");
    }

    PyAllowThreads _pythread;
    gaussianGradientMultiArray(srcMultiArrayRange(array), destMultiArray(res), opt);

    return res;
}


template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > array,
                                  const ConvolutionOptions<ndim-1> & opt,
                                  NumpyArray<ndim-1, Singleband<VoxelType> > res = NumpyArray<ndim-1, Singleband<VoxelType> >())
{
    using namespace vigra::functor;
    static const int sdim = ndim - 1;
    
    std::string description("Gaussian gradient magnitude");
    typedef typename MultiArrayShape<sdim>::type Shape;
    Shape tmpShape(array.shape().begin());
    if(opt.to_point != Shape())
        tmpShape = opt.to_point-opt.from_point;
    
    res.reshapeIfEmpty(array.taggedShape().resize(tmpShape).setChannelDescription(description), 
          "gaussianGradientMagnitude(): Output array has wrong shape.");
    res.init(VoxelType());
    
    PyAllowThreads _pythread;
    MultiArray<sdim, TinyVector<VoxelType, sdim> > grad(tmpShape);
    
    for(int k=0; k<array.shape(sdim); ++k)
    {
        MultiArrayView<sdim, VoxelType, StridedArrayTag> barray = array.bindOuter(k);
    
        gaussianGradientMultiArray(srcMultiArrayRange(barray), destMultiArray(grad), opt);
        combineTwoMultiArrays(srcMultiArrayRange(grad), srcMultiArray(res), destMultiArray(res), 
                              squaredNorm(Arg1())+Arg2());
    }
    transformMultiArray(srcMultiArrayRange(res), destMultiArray(res), sqrt(Arg1()));
    
    return res;
}


template < class PixelType>
NumpyAnyArray 
pythonRieszTransformOfLOG2D(NumpyArray<2, Singleband<PixelType> > image,
                            double scale, 
                            unsigned int xorder, unsigned int yorder,
                            NumpyArray<2, Singleband<PixelType> > res = NumpyArray<2, Singleband<PixelType> >())
{
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription("Riesz transform"), 
              "rieszTransformOfLOG2D(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    rieszTransformOfLOG(srcImageRange(image), destImage(res), scale, xorder, yorder);
     
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                  const ConvolutionOptions<ndim-1> & opt,
                                  NumpyArray<ndim, Multiband<VoxelType> > res = NumpyArray<ndim, Multiband<VoxelType> >())
{
    using namespace vigra::functor;
    static const int sdim = ndim - 1;
    
    std::string description("channel-wise Gaussian gradient magnitude");
    
    typedef typename MultiArrayShape<sdim>::type Shape;
    Shape tmpShape(volume.shape().begin());
    if(opt.to_point != Shape())
        tmpShape = opt.to_point-opt.from_point;
    
    res.reshapeIfEmpty(volume.taggedShape().resize(tmpShape).setChannelDescription(description), 
             "gaussianGradientMagnitude(): Output array has wrong shape.");
    
    PyAllowThreads _pythread;
    MultiArray<sdim, TinyVector<VoxelType, sdim> > grad(tmpShape);
    
    for(int k=0; k<volume.shape(sdim); ++k)
    {
        MultiArrayView<sdim, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
        MultiArrayView<sdim, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    
        gaussianGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(grad), opt);
        transformMultiArray(srcMultiArrayRange(grad), destMultiArray(bres), norm(Arg1()));
    }
    
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitude(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                python::object sigma, bool accumulate,
                                NumpyAnyArray res,
                                python::object sigma_d, 
                                python::object step_size,
                                double window_size = 0.0, 
                                python::object roi = python::object())
{
    pythonScaleParam<ndim - 1> params(sigma, sigma_d, step_size, "gaussianGradientMagnitude");
    params.permuteLikewise(volume);
    ConvolutionOptions<ndim-1> opt(params().filterWindowSize(window_size));
    
    typedef typename MultiArrayShape<ndim - 1>::type Shape;
    if(roi != python::object())
    {
        opt.subarray(volume.permuteLikewise(python::extract<Shape>(roi[0])()), 
                     volume.permuteLikewise(python::extract<Shape>(roi[1])()));
    }
    else
    {
        opt.subarray(Shape(), Shape(volume.shape().begin()));
    }
    
    return accumulate
              ? pythonGaussianGradientMagnitudeND(volume, opt, NumpyArray<ndim-1, Singleband<VoxelType> >(res))
              : pythonGaussianGradientMagnitudeND(volume, opt, NumpyArray<ndim, Multiband<VoxelType> >(res));
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonSymmetricGradientND(NumpyArray<ndim, Singleband<VoxelType> > volume,
                          double sigma,
                          NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res=python::object(),
                          python::object step_size = python::object(1.0), 
                          python::object roi = python::object())
{
    pythonScaleParam<ndim> params(python::object(0.0), python::object(0.0),
                                 step_size, "symmetricGradient");
    params.permuteLikewise(volume);
    ConvolutionOptions<ndim> opt(params());
    
    if(roi != python::object())
    {
        typedef typename MultiArrayShape<ndim>::type Shape;
        Shape start = volume.permuteLikewise(python::extract<Shape>(roi[0])());
        Shape stop  = volume.permuteLikewise(python::extract<Shape>(roi[1])());
        opt.subarray(start, stop);
        res.reshapeIfEmpty(volume.taggedShape().resize(stop-start).setChannelDescription("symmetric gradient"), 
                 "symmetricGradient(): Output array has wrong shape.");
    }
    else
    {
        res.reshapeIfEmpty(volume.taggedShape().setChannelDescription("symmetric gradient"), 
                 "symmetricGradient(): Output array has wrong shape.");
    }
    
    PyAllowThreads _pythread;
    symmetricGradientMultiArray(srcMultiArrayRange(volume), destMultiArray(res), opt);
    return res;
}

template < class VoxelType, unsigned int N >
NumpyAnyArray 
pythonHessianOfGaussianND(NumpyArray<N, Singleband<VoxelType> > array,
                          python::object sigma,
                          NumpyArray<N, TinyVector<VoxelType, int(N*(N+1)/2)> > res= NumpyArray<N, TinyVector<VoxelType, int(N*(N+1)/2)> >(),
                          python::object sigma_d = python::object(0.0), 
                          python::object step_size = python::object(1.0),
                          double window_size = 0.0, 
                          python::object roi = python::object())
{
    std::string description("Hessian of Gaussian (flattened upper triangular matrix), scale=");
    description += asString(sigma);
    
    pythonScaleParam<N> params(sigma, sigma_d, step_size, "hessianOfGaussian");
    params.permuteLikewise(array);
    ConvolutionOptions<N> opt(params().filterWindowSize(window_size));
    
    if(roi != python::object())
    {
        typedef typename MultiArrayShape<N>::type Shape;
        Shape start = array.permuteLikewise(python::extract<Shape>(roi[0])());
        Shape stop  = array.permuteLikewise(python::extract<Shape>(roi[1])());
        opt.subarray(start, stop);
        res.reshapeIfEmpty(array.taggedShape().resize(stop-start).setChannelDescription(description), 
               "hessianOfGaussian(): Output array has wrong shape.");
    }
    else
    {
        res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
               "hessianOfGaussian(): Output array has wrong shape.");
    }
    
    PyAllowThreads _pythread;
    hessianOfGaussianMultiArray(srcMultiArrayRange(array), destMultiArray(res), opt);
    return res;
}

#if 0 // FIXME: this is probably no longer needed thanks to axistags
template < class VoxelType>
NumpyAnyArray 
pythonHessianOfGaussian3D(NumpyArray<3, Singleband<VoxelType> > volume,
                          python::object sigma,
                          NumpyArray<3, TinyVector<VoxelType, 6> > res=NumpyArray<3, TinyVector<VoxelType, 6> >(),
                          python::object sigma_d = python::object(0.0), python::object step_size = python::object(1.0))
{
    pythonScaleParam<3> params(sigma, sigma_d, step_size, "hessianOfGaussian");
    params.permuteLikewise(volume);
    std::string description("Hessian of Gaussian (flattened upper triangular matrix), scale=");
    description += asString(sigma);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
          "hessianOfGaussian(): Output array has wrong shape.");
    
    PyAllowThreads _pythread;
    hessianOfGaussianMultiArray(srcMultiArrayRange(volume), destMultiArray(res), params());
    
    return res;
}

template < class PixelType>
NumpyAnyArray 
pythonHessianOfGaussian2D(NumpyArray<2, Singleband<PixelType> > image,
                          python::object sigma,
                          NumpyArray<2, TinyVector<PixelType, 3> > res=NumpyArray<2, TinyVector<PixelType, 3> >(),
                          python::object sigma_d = python::object(0.0), python::object step_size = python::object(1.0))
{
    pythonScaleParam<2> params(sigma, sigma_d, step_size, "hessianOfGaussian");
    params.permuteLikewise(image);
    std::string description("Hessian of Gaussian (flattened upper triangular matrix), scale=");
    description += asString(sigma);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
             "hessianOfGaussian(): Output array has wrong shape.");
    
    PyAllowThreads _pythread;
    hessianOfGaussianMultiArray(srcMultiArrayRange(image), destMultiArray(res), params());
    
    return res;
}
#endif

template <class PixelType, unsigned int N>
NumpyAnyArray 
pythonStructureTensor(NumpyArray<N, Multiband<PixelType> > array, 
                      python::object innerScale, python::object outerScale,
                      NumpyArray<N-1, TinyVector<PixelType, int(N*(N-1)/2)> > res=NumpyArray<N-1, TinyVector<PixelType, int(N*(N-1)/2)> >(),
                      python::object sigma_d = python::object(0.0), 
                      python::object step_size = python::object(1.0),
                      double window_size = 0.0, 
                      python::object roi = python::object())
{
    using namespace vigra::functor;
    static const int sdim = N - 1;
    
    std::string description("structure tensor (flattened upper triangular matrix), inner scale=");
    description += asString(innerScale) + ", outer scale=" + asString(outerScale);
    
    pythonScaleParam<N-1> params(innerScale, sigma_d, step_size, outerScale, "structureTensor");
    params.permuteLikewise(array);
    ConvolutionOptions<N-1> opt(params().filterWindowSize(window_size));
    
    if(roi != python::object())
    {
        typedef typename MultiArrayShape<N-1>::type Shape;
        Shape start = array.permuteLikewise(python::extract<Shape>(roi[0])());
        Shape stop  = array.permuteLikewise(python::extract<Shape>(roi[1])());
        opt.subarray(start, stop);
        res.reshapeIfEmpty(array.taggedShape().resize(stop-start).setChannelDescription(description), 
                     "structureTensor(): Output array has wrong shape.");
    }
    else
    {
        res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
                     "structureTensor(): Output array has wrong shape.");
    }
    
    PyAllowThreads _pythread;

    MultiArrayView<sdim, PixelType, StridedArrayTag> band = array.bindOuter(0);	
    structureTensorMultiArray(srcMultiArrayRange(band), destMultiArray(res), opt);
    
    if(array.shape(sdim) > 1)
    {
        MultiArray<sdim, TinyVector<PixelType, int(N*(N-1)/2)> > st(res.shape());
        
        for(int b=1; b<array.shape(sdim); ++b)
        {
            MultiArrayView<sdim, PixelType, StridedArrayTag> band = array.bindOuter(b);
            structureTensorMultiArray(srcMultiArrayRange(band), destMultiArray(st), opt);
            combineTwoMultiArrays(srcMultiArrayRange(res), srcMultiArray(st), 
                                  destMultiArray(res), Arg1() + Arg2());
        }
    }
    
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray 
pythonBoundaryTensor2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                       double scale,
                       NumpyArray<2, TinyVector<DestPixelType, 3> > res = NumpyArray<2, TinyVector<DestPixelType, 3> >())
{
    std::string description("boundary tensor (flattened upper triangular matrix), scale=");
    description += asString(scale);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
           "boundaryTensor2D(): Output array has wrong shape.");    

    PyAllowThreads _pythread;

    boundaryTensor(srcImageRange(image), destImage(res), scale);
     
    return res;
}


template < class SrcPixelType, typename DestPixelType  >
NumpyAnyArray 
pythonTensorEigenRepresentation2D(NumpyArray<2, TinyVector<SrcPixelType, 3> > image,
                                  NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    std::string description("tensor eigen representation (ev1, ev2, angle)");
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
                    "tensorEigenRepresentation2D(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    tensorEigenRepresentation(srcImageRange(image), destImage(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonVectorToTensor(NumpyArray<N, TinyVector<PixelType, int(N)> > array,
                     NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > res = python::object())
{
    std::string description("outer product tensor (flattened upper triangular matrix)");

    res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
            "vectorToTensor(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    vectorToTensorMultiArray(srcMultiArrayRange(array), destMultiArray(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonTensorTrace(NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > array,
                  NumpyArray<N, Singleband<PixelType> > res = python::object())
{
    std::string description("tensor trace");

    res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
           "tensorTrace(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    tensorTraceMultiArray(srcMultiArrayRange(array), destMultiArray(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonTensorDeterminant(NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > array,
                        NumpyArray<N, Singleband<PixelType> > res = python::object())
{
    std::string description("tensor determinant");

    res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
                "tensorDeterminant(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    tensorDeterminantMultiArray(srcMultiArrayRange(array), destMultiArray(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonTensorEigenvalues(NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > array,
                        NumpyArray<N, TinyVector<PixelType, int(N)> > res = python::object())
{
    std::string description("tensor eigenvalues");

    res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
                 "tensorEigenvalues(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    tensorEigenvaluesMultiArray(srcMultiArrayRange(array), destMultiArray(res));
    
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray 
pythonHourGlassFilter2D(NumpyArray<2, TinyVector<SrcPixelType, 3> > image,
                        double sigma, 
                        double rho,
                        NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    std::string description("hourglass tensor (flattened upper triangular matrix), scale=");
    description += asString(sigma) + ", rho=" + asString(rho);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "hourGlassFilter2D(): Output array has wrong shape.");    
    
    PyAllowThreads _pythread;
    hourGlassFilter(srcImageRange(image), destImage(res), sigma, rho);
     
    return res;
}

void defineTensor()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);
    
    def("gaussianGradient",
        registerConverters(&pythonGaussianGradientND<float,2>),
        (arg("image"), arg("sigma"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Calculate the gradient vector by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D scalar image.\n\n"
        "If 'sigma' is a single value, an isotropic filter at this scale is "
        "applied (i.e., each dimension is filtered in the same way). "
        "If 'sigma' is a tuple or list of values, the amount of smoothing "
        "will be different for each spatial dimension.\n"
        "The optional 'sigma_d' (single, tuple, or list) denotes the resolution standard deviation "
        "per axis, the optional 'step_size' (single, tuple, or list) the distance between two adjacent "
        "pixels for each dimension. "
        "The length of the tuples or lists must be equal to the "
        "number of spatial dimensions.\n\n"        
        "'window_size' and 'roi' have the same meaning as in :func:`gaussianSmoothing`.\n\n"
        "For details see gaussianGradientMultiArray_ and ConvolutionOptions_ in the vigra C++ documentation.\n");

    def("gaussianGradient",
        registerConverters(&pythonGaussianGradientND<float,3>),
        (arg("volume"), arg("sigma"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Likewise for a 3D scalar volume.\n");

    def("rieszTransformOfLOG2D",
        registerConverters(&pythonRieszTransformOfLOG2D<float>),
        (arg("image"), arg("scale"), arg("xorder"), arg("yorder"),arg("out")=python::object()),
        "Calculate Riesz transforms of the Laplacian of Gaussian.\n\n"
        "For details see rieszTransformOfLOG_ in the vigra C++ documentation.\n");

    def("gaussianGradientMagnitude",
        registerConverters(&pythonGaussianGradientMagnitude<float,3>),
        (arg("image"), arg("sigma"), arg("accumulate")=true, arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Calculate the gradient magnitude by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D scalar or multiband image.\n"
        "If 'accumulate' is True (the default), the gradients are accumulated (in the "
        "L2-norm sense) over all  channels of a multi-channel array. Otherwise, "
        "a separate gradient magnitude is computed for each channel.\n\n"
        "If 'sigma' is a single value, an isotropic filter at this scale is "
        "applied (i.e., each dimension is filtered in the same way). "
        "If 'sigma' is a tuple or list of values, the amount of smoothing "
        "will be different for each spatial dimension.\n"
        "The optional 'sigma_d' (single, tuple, or list) denotes the resolution standard deviation "
        "per axis, the optional 'step_size' (single, tuple, or list) the distance between two adjacent "
        "pixels for each dimension. "
        "The length of the tuples or lists must be equal to the "
        "number of spatial dimensions.\n\n"        
        "'window_size' and 'roi' have the same meaning as in :func:`gaussianSmoothing`.\n\n"
        "For details see gaussianGradientMultiArray_ and ConvolutionOptions_ in the vigra C++ documentation.\n");

    def("gaussianGradientMagnitude",
        registerConverters(&pythonGaussianGradientMagnitude<float,4>),
        (arg("volume"), arg("sigma"), arg("accumulate")=true, arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Likewise for a 3D scalar or multiband volume.\n");

    def("symmetricGradient",
        registerConverters(&pythonSymmetricGradientND<float,2>),
        (arg("image"), arg("out")=python::object(), arg("step_size")=1.0, arg("roi")=python::object()),
        "Calculate gradient of a scalar 2D image using symmetric difference filters."
        "\n"
        "The optional tuple or list 'step_size' denotes the distance between two "
        "adjacent pixels for each dimension; its length must be equal to the "
        "number of spatial dimensions.\n\n"        
        "'roi' has the same meaning as in :func:`gaussianSmoothing`.\n\n"
        "For details see symmetricGradientMultiArray_ and ConvolutionOptions_ in the vigra C++ documentation.\n");

    def("symmetricGradient",
        registerConverters(&pythonSymmetricGradientND<float,3>), 
        (arg("volume"), arg("out")=python::object(), arg("step_size")=1.0, arg("roi")=python::object()),
        "Likewise for a 3D scalar volume.\n");
    
    // FIXME: is this function still needed?
    def("hessianOfGaussian2D",
        registerConverters(&pythonHessianOfGaussianND<float, 2>),
        (arg("image"), arg("sigma"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Calculate the Hessian matrix by means of a derivative of "
        "Gaussian filters at the given scale for a 2D scalar image.\n"
        "\n"
        "If 'sigma' is a single value, an isotropic filter at this scale is "
        "applied (i.e., each dimension is filtered in the same way). "
        "If 'sigma' is a tuple or list of values, the amount of smoothing "
        "will be different for each spatial dimension.\n"
        "The optional 'sigma_d' (single, tuple, or list) denotes the resolution standard deviation "
        "per axis, the optional 'step_size' (single, tuple, or list) the distance between two adjacent "
        "pixels for each dimension. "
        "The length of the tuples or lists must be equal to the "
        "number of spatial dimensions.\n\n"        
        "'window_size' and 'roi' have the same meaning as in :func:`gaussianSmoothing`.\n\n"
        "For details see hessianOfGaussianMultiArray_ and ConvolutionOptions_ in the vigra C++ documentation.\n");

    // FIXME: is this function still needed?
    def("hessianOfGaussian3D",
        registerConverters(&pythonHessianOfGaussianND<float, 3>),
        (arg("volume"), arg("sigma"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Calculate the Hessian matrix by means of a derivative of "
        "Gaussian filters at the given scale for a 3D scalar image.\n"
        "\n"
        "For details see hessianOfGaussianMultiArray_ in the vigra C++ documentation.\n");

    def("hessianOfGaussian",
        registerConverters(&pythonHessianOfGaussianND<float,2>),
        (arg("image"), arg("sigma"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Calculate the Hessian matrix by means of a derivative of "
        "Gaussian filters at the given scale for a 2D scalar image.\n"
        "\n"
        "If 'sigma' is a single value, an isotropic filter at this scale is "
        "applied (i.e., each dimension is filtered in the same way). "
        "If 'sigma' is a tuple or list of values, the amount of smoothing "
        "will be different for each spatial dimension.\n"
        "The optional 'sigma_d' (single, tuple, or list) denotes the resolution standard deviation "
        "per axis, the optional 'step_size' (single, tuple, or list) the distance between two adjacent "
        "pixels for each dimension. "
        "The length of the tuples or lists must be equal to the "
        "number of spatial dimensions.\n\n"        
        "'window_size' and 'roi' have the same meaning as in :func:`gaussianSmoothing`.\n\n"
        "For details see hessianOfGaussianMultiArray_ in the vigra C++ documentation.\n");

    def("hessianOfGaussian",
        registerConverters(&pythonHessianOfGaussianND<float,3>),
        (arg("volume"), arg("sigma"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Likewise for a 3D scalar or multiband volume.\n");

    def("structureTensor",
        registerConverters(&pythonStructureTensor<float,3>),
        (arg("image"), arg("innerScale"), arg("outerScale"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Calculate the structure tensor of an image by means of Gaussian "
        "(derivative) filters at the given scales. If the input has multiple channels, "
        "the structure tensors of each channel are added to get the result.\n\n"
        "If 'innerScale' and 'outerScale' are single values, "
        "isotropic filters at these scales are "
        "applied (i.e., each dimension is filtered in the same way). "
        "If 'innerScale' and / or 'outerScale' are are tuples or lists of "
        "values, the amount of smoothing "
        "will be different for each spatial dimension.\n"
        "The optional 'sigma_d' (single, tuple, or list) denotes the resolution standard deviation "
        "per axis, the optional 'step_size' (single, tuple, or list) the distance between two adjacent "
        "pixels for each dimension. "
        "The length of the tuples or lists must be equal to the "
        "number of spatial dimensions.\n\n"        
        "'window_size' and 'roi' have the same meaning as in :func:`gaussianSmoothing`.\n\n"
        "For details see structureTensorMultiArray_ and ConvolutionOptions_ in the vigra C++ documentation.\n");

    def("structureTensor",
        registerConverters(&pythonStructureTensor<float,4>),
        (arg("volume"), arg("innerScale"), arg("outerScale"), arg("out")=python::object(), 
         arg("sigma_d")=0.0, arg("step_size")=1.0, arg("window_size")=0.0, arg("roi")=python::object()),
        "Likewise for a 3D scalar or multiband volume.\n");

    def("boundaryTensor2D",
        registerConverters(&pythonBoundaryTensor2D<float, float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Calculate the boundary tensor for a scalar valued 2D image."
        "For details see boundaryTensor_ in the vigra C++ documentation.\n");
        
    /** FIXME: Export of Kernel2D before
    def("gradientEnergyTensor2D",
        registerConverters(&gradientEnergyTensor2D<float,float>),
        (arg("image"), arg("derivKernel"), arg("smoothKernel"),arg("out")=python::object()));
        
    */

    def("tensorEigenRepresentation2D",
        registerConverters(&pythonTensorEigenRepresentation2D<float,float>),
        (arg("image"),arg("out")=python::object()),
        "Calculate eigen representation of a symmetric 2x2 tensor.\n\n"
        "For details see tensorEigenRepresentation_ in the vigra C++ documentation.\n"
        );

    def("vectorToTensor",
        registerConverters(&pythonVectorToTensor<float,2>),
        (arg("image"),arg("out")=python::object()),
        "Turn a 2D vector valued image (e.g. the gradient image) into "
        "a tensor image by computing the outer product in every pixel.\n\n"
        "For details see vectorToTensorMultiArray_ in the vigra C++ documentation.\n");

    def("vectorToTensor",
        registerConverters(&pythonVectorToTensor<float,3>),
        (arg("volume"),arg("out")=python::object()),
        "Likewise for a 3D vector-valued volume.\n");

    def("tensorTrace",
        registerConverters(&pythonTensorTrace<float,2>),
        (arg("image"),arg("out")=python::object()),
        "Calculate the trace of a 2x2 tensor image.\n\n"
        "For details see tensorTraceMultiArray_ in the vigra C++ documentation.\n");

    def("tensorTrace",
        registerConverters(&pythonTensorTrace<float,3>),
        (arg("volume"),arg("out")=python::object()),
        "Likewise for a 3x3 tensor volume.\n");

    def("tensorDeterminant",
        registerConverters(&pythonTensorDeterminant<float,2>),
        (arg("image"),arg("out")=python::object()),
        "Calculate the determinant of a 2x2 tensor image.\n\n"
        "For details see tensorDeterminantMultiArray_ in the vigra C++ documentation.\n");

    def("tensorDeterminant",
        registerConverters(&pythonTensorDeterminant<float,3>),
        (arg("volume"),arg("out")=python::object()),
        "Likewise for a 3x3 tensor volume.\n");

    def("tensorEigenvalues",
        registerConverters(&pythonTensorEigenvalues<float,2>),
        (arg("image"),arg("out")=python::object()),
        "Calculate the eigenvalues in each pixel/voxel of a 2x2 tensor image.\n\n"
        "For details see tensorEigenvaluesMultiArray_ in the vigra C++ documentation.\n");

    def("tensorEigenvalues",
        registerConverters(&pythonTensorEigenvalues<float,3>),
        (arg("volume"),arg("out")=python::object()),
        "Likewise for a 3x3 tensor volume.\n");

    def("hourGlassFilter2D",
        registerConverters(&pythonHourGlassFilter2D<float,float>),
        (arg("image"), arg("sigma"), arg("rho"),arg("out")=python::object()),
        "Anisotropic tensor smoothing with the hourglass filter. \n\n"
        "For details see hourGlassFilter_ in the vigra C++ documentation.\n");
 
 /* Wee, tons of errors here
    def("ellipticGaussian2D",
        registerConverters(&ellipticGaussian2D<float,float>),
        (arg("image"), arg("sigmamax"), arg("sigmamin"),arg("out")=python::object()));
    def("ellipticGaussian2D",
        registerConverters(&ellipticGaussian2D<float,float>),
        (arg("image"), arg("sigmamax"), arg("sigmamin"),arg("out")=python::object()));
  */
}

} // namespace vigra
