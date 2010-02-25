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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyfilters_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/boundarytensor.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>
#include <vigra/symmetry.hxx>

#include <cmath>

#include "tensors.hxx"
namespace python = boost::python;

namespace vigra
{
namespace { void DUMMY_FUNCTION(int, int, int, int, int, int, int, int, int, int) {} }


template < class PixelType>
NumpyAnyArray pythonEigHessian2d(NumpyArray<2,Singleband<PixelType> > image, double scale,
                                 NumpyArray<3,Singleband<PixelType>  > res)
{
    //reshape if empty
    res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(0),image.shape(1),2));

    MultiArray<3,PixelType> hessian(MultiArrayShape<3>::type(image.shape(0),image.shape(1),3));

    hessianOfGaussian(image, hessian, scale);

    eigenValuesPerPixel(hessian,res);

    return res;
}

template<class PixelType>
NumpyAnyArray pythonEigStructureTensor2d(NumpyArray<2,Singleband<PixelType> > image, double innerScale,double outerScale,
                                         NumpyArray<3,Singleband<PixelType> > res)
{
    //reshape output
    res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(0),image.shape(1),2));

    MultiArray<3,PixelType> st(MultiArrayShape<3>::type(image.shape(0),image.shape(1), 3));

    structureTensor(image,st,innerScale,outerScale);

    eigenValuesPerPixel(st,res);
    return res;
}

template<class PixelType>
NumpyAnyArray pythonGaussianSmooth2d(NumpyArray<2,Singleband<PixelType> > image, double scale,
                                     NumpyArray<2,Singleband<PixelType> > res)
{
    //reshape output
    res.reshapeIfEmpty(image.shape());

    gaussianSmoothMultiArray(srcMultiArrayRange(image), destMultiArray(res), scale);
   

    return res;
}

template < class VoxelType >
NumpyAnyArray pythonDistanceTransform3D(NumpyArray<3, Singleband<VoxelType> > volume, 
                                  bool background,
                                  NumpyArray<3, Singleband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "distanceTransform3D(): Output array has wrong shape.");
    
    separableMultiDistance(srcMultiArrayRange(volume),
        destMultiArray(res), background);
    return res;
}

template < class PixelType>
NumpyAnyArray pythonRieszTransformOfLOG2D(NumpyArray<2, Singleband<PixelType> > image,
                                    double scale, unsigned int xorder,
                                    unsigned int yorder,
                                    NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "rieszTransformOfLOG2D(): Output array has wrong shape.");    
    
    rieszTransformOfLOG(srcImageRange(image), destImage(res),
        scale, xorder, yorder);
     
    return res;
}

template < class PixelType, typename DestPixelType >
NumpyAnyArray pythonDistanceTransform2D(NumpyArray<2, Singleband<PixelType> > image,
                                  PixelType background, 
                                  int norm,
                                  NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "distanceTransform2D(): Output array has wrong shape.");
    distanceTransform(srcImageRange(image), destImage(res), background,
        norm);
    return res;
}



template < class SrcPixelType >
NumpyAnyArray pythonRadialSymmetryTransform2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                        double scale = 1,
                                        NumpyArray<2, Singleband<SrcPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "radialSymmetryTransform2D(): Output array has wrong shape.");    
    radialSymmetryTransform(srcImageRange(image), destImage(res), scale);
    return res;
}



void defineFilters2D()
{
    using namespace python;

    def("eigHessian2d",
        registerConverters(&pythonEigHessian2d<float>),
        (arg("image"), arg("scale"), arg("out")=object()),
        "computes the eigenvalues of the Hessian in each pixel of an image."
        "\n"
        "The 'scale' parameter gives the scale of 2nd derivative filter.\n"
       );
    def("eigStructureTensor2d",
        registerConverters(&pythonEigStructureTensor2d<float>),
        (arg("image"), arg("scale"), arg("out")=object()),
        "computes the eigenvalues of the structure tensor in each pixel of an image."
        "\n"
        "The innerScale parameter gives the scale of the gradient filter."
        "The outerScale parameter gives the outer scale of the smoothing filter\n"
       );
    def("gaussianSmooth2d",
        registerConverters(&pythonGaussianSmooth2d<float>),
        (arg("image"), arg("scale"), arg("out")=object()),
        "computes the gradient magnitude in each pixel of an image.\n"
        "\n"
        "The scale parameter gives the scale of the derivative filter.\n"
       );

    def("distanceTransform3D",
        registerConverters(&pythonDistanceTransform3D<npy_int32>),
        (arg("array"), arg("background"),
         arg("out")=python::object()),
        "For all background pixels, calculate the distance to the nearest object or contour."
        "The label of the pixels to be considered background in the source image is passed "
        "in the parameter 'background'."
        "Source pixels with other labels will be considered objects."
        "In the destination image, all pixels corresponding to background will be assigned "
        "the their distance value, all pixels corresponding to objects will be assigned 0.\n"
        "\n"
        "For more details see separableMultiDistance_ in the vigra C++ documentation.\n");

    def("rieszTransformOfLOG2D",
        registerConverters(&pythonRieszTransformOfLOG2D<float>),        // also multiband
        (arg("image"), arg("scale"), arg("xorder"), arg("yorder"),arg("out")=python::object()),
        "Calculate Riesz transforms of the Laplacian of Gaussian.\n\n"
        "For details see rieszTransformOfLOG_ in the vigra C++ documentation.\n");

    def("distanceTransform2D",
        registerConverters(&pythonDistanceTransform2D<npy_int32, float>),
        (arg("image"), 
         arg("background")=0, 
         arg("norm")=2,
         arg("out")=python::object()),
        "For all background pixels, calculate the distance to the nearest object or contour. "
        "The label of the pixels to be considered background in the source image is passed "
        "in the parameter 'background'. "
        "Source pixels with other labels will be considered objects. "
        "In the destination image, all pixels corresponding to background will be assigned "
        "the their distance value, all pixels corresponding to objects will be assigned 0.\n\n"
        "The norm parameter gives the distance norm to use.\n\n"
        "For details see distanceTransform_ in the vigra C++ documentation.\n");
    def("distanceTransform2D",
        registerConverters(&pythonDistanceTransform2D<UInt8,float>),
        (arg("image"), 
         arg("background")=0, 
         arg("norm")=2,
         arg("out")=python::object()));

    def("radialSymmetryTransform2D",
        registerConverters(&pythonRadialSymmetryTransform2D<float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find centers of radial symmetry in an 2D image.\n\n"
        "This algorithm implements the Fast Radial Symmetry Transform according to "
        "[G. Loy, A. Zelinsky: \"A Fast Radial Symmetry Transform for Detecting Points of Interest\", "
        "in: A. Heyden et al. (Eds.): Proc. of 7th European Conf. on Computer Vision, Part 1, pp. 358-368, Springer LNCS 2350, 2002]\n\n"
        "For details see radialSymmetryTransform_ in the vigra C++ documentation.\n");
}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(filters)
{
    import_vigranumpy();
    defineFilters2D();
}
