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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyanalysis_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/cornerdetection.hxx>
#include <vigra/boundarytensor.hxx>
#include <vigra/mathutil.hxx>

namespace python = boost::python;

namespace vigra
{

template < class PixelType >
NumpyAnyArray 
pythonCornerResponseFunction2D(NumpyArray<2, Singleband<PixelType> > image,
                               double scale=1.0,
                               NumpyArray<2, Singleband<PixelType> > res = python::object() )
{
    res.reshapeIfEmpty(image.shape(), "cornernessHarris(): Output array has wrong shape.");    
    
    cornerResponseFunction(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonFoerstnerCornerDetector2D(NumpyArray<2, Singleband<PixelType> > image,
                                double scale=1.0,
                                NumpyArray<2, Singleband<PixelType> > res = python::object() )
{
    res.reshapeIfEmpty(image.shape(), "cornernessFoerstner(): Output array has wrong shape.");    
    
    foerstnerCornerDetector(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonRohrCornerDetector2D(NumpyArray<2, Singleband<PixelType> > image,
                           double scale = 1.0,
                           NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "cornernessRohr(): Output array has wrong shape.");    
    
    rohrCornerDetector(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonBeaudetCornerDetector2D(NumpyArray<2, Singleband<PixelType> > image,
                              double scale=1.0,
                              NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "cornernessBeaudet(): Output array has wrong shape.");    
    
    beaudetCornerDetector(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonBoundaryTensorCornerDetector2D(NumpyArray<2, Singleband<PixelType> > image,
                                     double scale=1.0,
                                     NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    MultiArrayShape<2>::type shape(image.shape());
    
    res.reshapeIfEmpty(shape, "cornernessBoundaryTensor(): Output array has wrong shape.");    
    
    MultiArray<2, TinyVector<PixelType, 3> > bt(shape);
    boundaryTensor(srcImageRange(image), destImage(bt), scale);
    
    PixelType ev1, ev2;
    for(int y=0; y<shape[1]; ++y)
    {
        for(int x=0; x<shape[0]; ++x)
        {
            symmetric2x2Eigenvalues(bt(x,y)[0], bt(x,y)[1], bt(x,y)[2], &ev1, &ev2);
            res(x,y) = PixelType(2.0)*ev2;
        }
    }
    return res;
}


void defineInterestpoints()
{
    using namespace python;

    def("cornernessHarris",
        registerConverters(&pythonCornerResponseFunction2D<float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in a scalar 2D image using the method of Harris at the given 'scale'.\n\n"
        "For details see cornerResponseFunction_ in the vigra C++ documentation.\n"
        );

    def("cornernessFoerstner",
        registerConverters(&pythonFoerstnerCornerDetector2D<float>), 
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in a scalar 2D image using the method of Foerstner at the given 'scale'.\n\n"
        "For details see foerstnerCornerDetector_ in the vigra C++ documentation.\n");

    def("cornernessRohr",
        registerConverters(&pythonRohrCornerDetector2D<float>), 
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in a scalar 2D image using the method of Rohr at the given 'scale'.\n\n"
        "For details see rohrCornerDetector_ in the vigra C++ documentation.\n");

    def("cornernessBeaudet",
        registerConverters(&pythonBeaudetCornerDetector2D<float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in a scalar 2D image using the method of Beaudet at the given 'scale'.\n\n"
        "For details see beaudetCornerDetector_ in the vigra C++ documentation.\n");

    def("cornernessBoundaryTensor",
        registerConverters(&pythonBoundaryTensorCornerDetector2D<float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in a scalar 2D image using the boundary tensor at the given 'scale'.\n\n"
        "Specifically, the cornerness is defined as twice the small eigenvalue "
        "of the boundary tensor.\n\n"
        "For details see boundaryTensor_ in the vigra C++ documentation.\n");
}

} // namespace vigra
