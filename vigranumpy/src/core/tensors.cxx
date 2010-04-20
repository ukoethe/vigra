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

namespace python = boost::python;

namespace vigra
{

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonGaussianGradientND(NumpyArray<ndim, Singleband<VoxelType> > volume,
                                       double sigma,
                                       NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "gaussianGradient(): Output array has wrong shape.");
    gaussianGradientMultiArray(srcMultiArrayRange(volume), destMultiArray(res), sigma);
    return res;
}


template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                  double sigma, 
                                  NumpyArray<ndim-1, Singleband<VoxelType> > res=python::object())
{
    using namespace vigra::functor;
    
    typename MultiArrayShape<ndim-1>::type tmpShape(volume.shape().begin());
    res.reshapeIfEmpty(tmpShape, "gaussianGradientMagnitude(): Output array has wrong shape.");
    res.init(NumericTraits<VoxelType>::zero());
    MultiArray<ndim-1, TinyVector<VoxelType, (int)(ndim-1)> > grad(tmpShape);
    for(int k=0; k<volume.shape(ndim-1); ++k)
    {
        MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    
        gaussianGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(grad), sigma);
        combineTwoMultiArrays(srcMultiArrayRange(grad), srcMultiArray(res), destMultiArray(res), 
                              squaredNorm(Arg1())+Arg2());
    }
    transformMultiArray(srcMultiArrayRange(res), destMultiArray(res), sqrt(Arg1()));
    return res;
}


template < class PixelType>
NumpyAnyArray 
pythonRieszTransformOfLOG2D(NumpyArray<2, Singleband<PixelType> > image,
                            double scale, unsigned int xorder,
                            unsigned int yorder,
                            NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "rieszTransformOfLOG2D(): Output array has wrong shape.");    
    
    rieszTransformOfLOG(srcImageRange(image), destImage(res),
        scale, xorder, yorder);
     
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                  double sigma, 
                                  NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    using namespace vigra::functor;
    
    res.reshapeIfEmpty(volume.shape(), "gaussianGradientMagnitude(): Output array has wrong shape.");
    
    typename MultiArrayShape<ndim-1>::type tmpShape(volume.shape().begin());
    MultiArray<ndim-1, TinyVector<VoxelType, (int)(ndim-1)> > grad(tmpShape);
    for(int k=0; k<volume.shape(ndim-1); ++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    
        gaussianGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(grad), sigma);
        transformMultiArray(srcMultiArrayRange(grad), destMultiArray(bres), norm(Arg1()));
    }
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitude(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                double sigma, bool accumulate, 
                                NumpyAnyArray res=python::object())
{
    return accumulate
              ? pythonGaussianGradientMagnitudeND(volume, sigma, NumpyArray<ndim-1, Singleband<VoxelType> >(res))
              : pythonGaussianGradientMagnitudeND(volume, sigma, NumpyArray<ndim, Multiband<VoxelType> >(res));
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonSymmetricGradientND(NumpyArray<ndim, Singleband<VoxelType> > volume,
                                        double sigma,
                                        NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "symmetricGradient(): Output array has wrong shape.");
    symmetricGradientMultiArray(srcMultiArrayRange(volume), destMultiArray(res));
    return res;
}

template < class VoxelType, unsigned int N >
NumpyAnyArray 
pythonHessianOfGaussianND(NumpyArray<N, Singleband<VoxelType> > volume,
                          double sigma,
                          NumpyArray<N, TinyVector<VoxelType, int(N*(N-1)/2)> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "hessianOfGaussian(): Output array has wrong shape.");
    hessianOfGaussianMultiArray(srcMultiArrayRange(volume), destMultiArray(res), sigma);
    return res;
}

template <class PixelType, unsigned int N>
NumpyAnyArray 
pythonStructureTensor(NumpyArray<N, Multiband<PixelType> > image, 
                      double innerScale, double outerScale,
                      NumpyArray<N-1, TinyVector<PixelType, int(N*(N-1)/2)> > res=python::object() )
{
    using namespace vigra::functor;
    
    res.reshapeIfEmpty(typename MultiArrayShape<N-1>::type(image.shape().begin()), 
                 "structureTensor(): Output array has wrong shape.");
    
    MultiArrayView<N-1, PixelType, StridedArrayTag> band = image.bindOuter(0);
    structureTensorMultiArray(srcMultiArrayRange(band), destMultiArray(res), 
                              innerScale, outerScale);
    
    if(image.shape(N-1) > 1)
    {
        MultiArray<N-1, TinyVector<PixelType, int(N*(N-1)/2)> > st(res.shape());
        for(int b=1; b<image.shape(N-1); ++b)
        {
            MultiArrayView<N-1, PixelType, StridedArrayTag> band = image.bindOuter(b);
            structureTensorMultiArray(srcMultiArrayRange(band), destMultiArray(st), 
                                      innerScale, outerScale);
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
                       NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "boundaryTensor2D(): Output array has wrong shape.");    

    boundaryTensor(srcImageRange(image), destImage(res), scale);
     
    return res;
}


template < class SrcPixelType, typename DestPixelType  >
NumpyAnyArray 
pythonTensorEigenRepresentation2D(NumpyArray<2, TinyVector<SrcPixelType, 3> >image,
                                  NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "tensorEigenRepresentation2D(): Output array has wrong shape.");    
    
    tensorEigenRepresentation(srcImageRange(image), destImage(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonVectorToTensor(NumpyArray<N, TinyVector<PixelType, int(N)> > image,
                     NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "vectorToTensor(): Output array has wrong shape.");    
    
    vectorToTensorMultiArray(srcMultiArrayRange(image), destMultiArray(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonTensorTrace(NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > image,
                  NumpyArray<N, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "tensorTrace(): Output array has wrong shape.");    
    
    tensorTraceMultiArray(srcMultiArrayRange(image), destMultiArray(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonTensorDeterminant(NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > image,
                        NumpyArray<N, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "tensorDeterminant(): Output array has wrong shape.");    
    
    tensorDeterminantMultiArray(srcMultiArrayRange(image), destMultiArray(res));
     
    return res;
}

// FIXME: generalize to handle non-interleaved representations
template < class PixelType, unsigned int N >
NumpyAnyArray 
pythonTensorEigenvalues(NumpyArray<N, TinyVector<PixelType, int(N*(N+1)/2)> > image,
                        NumpyArray<N, TinyVector<PixelType, int(N)> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "tensorEigenvalues(): Output array has wrong shape.");    
    
    tensorEigenvaluesMultiArray(srcMultiArrayRange(image), destMultiArray(res));
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonHourGlassFilter2D(NumpyArray<2, TinyVector<SrcPixelType, 3> >image,
                                double sigma, 
                                double rho,
                                NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "hourGlassFilter2D(): Output array has wrong shape.");    
    
    hourGlassFilter(srcImageRange(image), destImage(res), sigma, rho);
     
    return res;
}

void defineTensor()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);
    
    def("gaussianGradient",
    	registerConverters(&pythonGaussianGradientND<float,2>),
    	(arg("image"), arg("sigma"), arg("out")=python::object()),
        "Calculate the gradient vector by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D scalar image.\n"
        "\n"
        "For details see gaussianGradientMultiArray_ in the vigra C++ documentation.\n");

    def("gaussianGradient",
    	registerConverters(&pythonGaussianGradientND<float,3>),
    	(arg("volume"), arg("sigma"), arg("out")=python::object()),
        "Likewise for a 3D scalar volume.\n");

    def("rieszTransformOfLOG2D",
        registerConverters(&pythonRieszTransformOfLOG2D<float>),
        (arg("image"), arg("scale"), arg("xorder"), arg("yorder"),arg("out")=python::object()),
        "Calculate Riesz transforms of the Laplacian of Gaussian.\n\n"
        "For details see rieszTransformOfLOG_ in the vigra C++ documentation.\n");

    def("gaussianGradientMagnitude",
    	registerConverters(&pythonGaussianGradientMagnitude<float,3>),
    	(arg("image"), arg("sigma"), arg("accumulate")=true, arg("out")=python::object()),
        "Calculate the gradient magnitude by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D scalar or multiband image.\n"
        "If 'accumulate' is True (the default), the gradients are accumulated (in the "
        "L2-norm sense) over all  channels of a multi-channel array. Otherwise, "
        "a separate gradient magnitude is computed for each channel.\n"
        "\n"
        "For details see gaussianGradientMultiArray_ in the vigra C++ documentation.\n");

    def("gaussianGradientMagnitude",
    	registerConverters(&pythonGaussianGradientMagnitude<float,4>),
    	(arg("volume"), arg("sigma"), arg("accumulate")=true, arg("out")=python::object()),
        "Likewise for a 3D scalar or multiband volume.\n");

    def("symmetricGradient",
        registerConverters(&pythonSymmetricGradientND<float,2>),
        (arg("image"), arg("out")=python::object()),
        "Calculate gradient of a scalar 2D image using symmetric difference filters."
        "\n"
        "For details see symmetricGradientMultiArray_ in the vigra C++ documentation.\n");

    def("symmetricGradient",
        registerConverters(&pythonSymmetricGradientND<float,3>), 
        (arg("volume"), arg("out")=python::object()),
        "Likewise for a 3D scalar volume.\n");
    
    def("hessianOfGaussian",
    	registerConverters(&pythonHessianOfGaussianND<float,2>),
    	(arg("image"), arg("sigma"), arg("out")=python::object()),
        "Calculate the Hessian matrix by means of a derivative of "
        "Gaussian filters at the given scale for a 2D scalar image.\n"
        "\n"
        "For details see hessianOfGaussianMultiArray_ in the vigra C++ documentation.\n");

    def("hessianOfGaussian",
    	registerConverters(&pythonHessianOfGaussianND<float,3>),
    	(arg("volume"), arg("sigma"), arg("out")=python::object()),
        "Likewise for a 3D scalar volume.\n");

    def("structureTensor",
    	registerConverters(&pythonStructureTensor<float,3>),
    	(arg("image"), arg("innerScale"), arg("outerScale"), arg("out")=python::object()),
        "Calculate the structure tensor of an image by means of Gaussian "
        "(derivative) filters at the given scales. If the input has multiple channels, "
        "the structure tensors of each channel are added to get the result.\n"
        "\n"
        "For details see structureTensorMultiArray_ in the vigra C++ documentation.\n");

    def("structureTensor",
    	registerConverters(&pythonStructureTensor<float,4>),
    	(arg("volume"), arg("innerScale"), arg("outerScale"), arg("out")=python::object()),
        "Likewise for a 3D scalar or multiband volume.\n");

    def("boundaryTensor2D",
        registerConverters(&pythonBoundaryTensor2D<float, float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Calculate the boundary tensor for a scalar valued 2D image."
        "For details see boundaryTensor_ in the vigra C++ documentation.\n");
        
    /** Export of Kernel2D before
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
