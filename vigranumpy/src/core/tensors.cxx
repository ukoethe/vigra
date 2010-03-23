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
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/orientedtensorfilters.hxx>
#include <vigra/tensorutilities.hxx>
#include <vigra/boundarytensor.hxx>

#include <cmath>

namespace python = boost::python;

namespace vigra
{

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonBoundaryTensor2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                               double scale,
                               NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "boundaryTensor2D(): Output array has wrong shape.");    

    boundaryTensor(srcImageRange(image), destImage(res), scale);
     
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



template < class SrcPixelType, typename DestPixelType  >
NumpyAnyArray pythonTensorEigenRepresentation2D(NumpyArray<2, TinyVector<SrcPixelType, 3> >image,
                                          NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "tensorEigenRepresentation2D(): Output array has wrong shape.");    
    
    tensorEigenRepresentation(srcImageRange(image), destImage(res));
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonTensorTrace2D(NumpyArray<2, TinyVector<SrcPixelType, 3> > image,
                            NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "tensorTrace2D(): Output array has wrong shape.");    
    
    tensorTrace(srcImageRange(image), destImage(res));
     
    return res;
}



void defineTensor()
{
    using namespace python;
    

    
    def("boundaryTensor2D",
        registerConverters(&pythonBoundaryTensor2D<float, float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Calculate the boundary tensor for a scalar valued 2D image."
        "For details see boundaryTensor_ in the vigra C++ documentation.");
    /** Export of Kernel2D before
  def("gradientEnergyTensor2D",
        registerConverters(&gradientEnergyTensor2D<float,float>),
        (arg("image"), arg("derivKernel"), arg("smoothKernel"),arg("out")=python::object()));
        
    */


    def("hourGlassFilter2D",
        registerConverters(&pythonHourGlassFilter2D<float,float>),
        (arg("image"), arg("sigma"), arg("rho"),arg("out")=python::object()),
        "Anisotropic tensor smoothing with the hourglass filter. \n\n"
        "For details see hourGlassFilter_ in the vigra C++ documentation.");
 
    def("tensorEigenRepresentation2D",
        registerConverters(&pythonTensorEigenRepresentation2D<float,float>),                 // change coordinate system!
        (arg("image"),arg("out")=python::object()),
        "Calculate eigen representation of a symmetric 2x2 tensor.\n\n"
        "For details see tensorEigenRepresentation_ in the vigra C++ documentation."
        );

    def("tensorTrace2D",
        registerConverters(&pythonTensorTrace2D<float,float>),
        (arg("image"),arg("out")=python::object()),
        "Calculate the trace of a 2x2 tensor.\n\n"
        "For details see tensorTrace_ in the vigra C++ documentation.");
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
