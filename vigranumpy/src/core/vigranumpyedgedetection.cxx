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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyedgedetection_PyArray_API
//#define NO_IMPORT_ARRAY
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/cornerdetection.hxx>
#include <vigra/edgedetection.hxx>

#include <cmath>

namespace python = boost::python;

namespace vigra
{

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonCornerResponseFunction2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                       double scale=1,
                                       NumpyArray<2, Singleband<DestPixelType> > res = python::object() )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "cornerResponseFunction2D(): Output array has wrong shape.");    
    
    cornerResponseFunction(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonFoerstnerCornerDetector2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                        double scale=1,
                                        NumpyArray<2, Singleband<DestPixelType> > res = python::object() )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "foerstnerCornerDetector2D(): Output array has wrong shape.");    
    
    foerstnerCornerDetector(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonRohrCornerDetector2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                   double scale = 1,
                                   NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "rohrCornerDetector2D(): Output array has wrong shape.");    
    
    rohrCornerDetector(srcImageRange(image), destImage(res), scale);
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonBeaudetCornerDetector2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                      double scale=1,
                                      NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "beaudetCornerDetector2D(): Output array has wrong shape.");    
    
    beaudetCornerDetector(srcImageRange(image), destImage(res), scale);
    return res;
}



template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonCannyEdgeImage(NumpyArray<2, Singleband<SrcPixelType> > image,
    double scale, double threshold, DestPixelType edgeMarker, NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "cannyEdgeImage(): Output array has wrong shape.");    
    
    cannyEdgeImage(srcImageRange(image), destImage(res), scale,
        threshold, edgeMarker);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonCannyEdgeImageWithThinning(NumpyArray<2, Singleband<SrcPixelType> > image,
       double scale, double threshold, DestPixelType edgeMarker, bool addBorder = true, NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "cannyEdgeImageWithThinning(): Output array has wrong shape.");    
    
    cannyEdgeImageWithThinning(srcImageRange(image), destImage(res),
        scale, threshold, edgeMarker, addBorder);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonShenCastanEdgeImage(NumpyArray<2, Singleband<SrcPixelType> > image,
      double scale, double threshold, DestPixelType edgeMarker, NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "pythonShenCastanEdgeImage(): Output array has wrong shape.");    
    
    differenceOfExponentialEdgeImage(srcImageRange(image),
        destImage(res), scale, threshold, edgeMarker);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonShenCastanCrackEdgeImage(NumpyArray<2, Singleband<SrcPixelType> > image,
      double scale, double threshold, DestPixelType edgeMarker, NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0)*2-1, image.shape(1)*2-1), "shenCastanCrackEdgeImage(): Output array has wrong shape.");    
    
    differenceOfExponentialCrackEdgeImage(srcImageRange(image),
        destImage(res), scale, threshold, edgeMarker);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonRemoveShortEdges(NumpyArray<2, Singleband<SrcPixelType> > image,
      int minEdgeLength, DestPixelType nonEdgeMarker, NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "shenCastanCrackEdgeImage(): Output array has wrong shape.");    
    
    copyImage(srcImageRange(image), destImage(res));
    removeShortEdges(destImageRange(res), minEdgeLength, nonEdgeMarker);
     
    return res;
}

template < class PixelType >
NumpyAnyArray pythonBeautifyCrackEdgeImage(NumpyArray<2, Singleband<PixelType> > image,
                                           PixelType edgeMarker, 
                                           PixelType backgroundMarker,
                                           NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "beautifyCrackEdgeImage(): Output array has wrong shape.");    
    
    copyImage(srcImageRange(image), destImage(res));
    beautifyCrackEdgeImage(destImageRange(res), edgeMarker,
        backgroundMarker);
     
    return res;
}

template < class PixelType >
NumpyAnyArray pythonCloseGapsInCrackEdgeImage(NumpyArray<2, Singleband<PixelType> > image,
                                              PixelType edgeMarker,
                                              NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "closeGapsInCrackEdgeImage(): Output array has wrong shape.");    
    
    copyImage(srcImageRange(image), destImage(res));
    closeGapsInCrackEdgeImage(destImageRange(res), edgeMarker);
     
    return res;
}

template < class PixelType >
NumpyAnyArray pythonRegionImageToCrackEdgeImage2D(NumpyArray<2, Singleband<PixelType> > image,
                                            PixelType edgeLabel = 0,
                                            NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0)*2-1, image.shape(1)*2-1), "regionImageToCrackEdgeImage2D(): Output array has wrong shape. Needs to be (w,h)*2 -1");


    regionImageToCrackEdgeImage(srcImageRange(image), destImage(res),
        edgeLabel);
    return res;
}

template < class PixelType >
NumpyAnyArray pythonRegionImageToEdgeImage2D(NumpyArray<2, Singleband<PixelType> > image,
                                       PixelType edgeLabel = 0,
                                       NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "regionImageToEdgeImage2D(): Output array has wrong shape.");

    regionImageToEdgeImage(srcImageRange(image), destImage(res),
        edgeLabel);
    return res;
}


void defineEdgedetection()
{
    using namespace python;
    

    def("cornerResponseFunction2D",
        registerConverters(&pythonCornerResponseFunction2D<float, float>),               // also multiband!
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in an 2D image.\n\n"
        "This algorithm implements the so called 'corner response function' to measure the 'cornerness' of each pixel in the image,"
        "according to [C.G. Harris and M.J. Stevens: \"A Combined Corner and Edge Detector\","
        "Proc. of 4th Alvey Vision Conference, 1988].\n"
        "Several studies have found this to be a very robust corner detector, although it moves the corners somewhat into one region,"
        " depending on the scale.\n\n"
        "For details see cornerResponseFunction_ in the vigra C++ documentation."
        );

    def("foerstnerCornerDetector2D",
        registerConverters(&pythonFoerstnerCornerDetector2D<float, float>),               // also multiband
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in an 2D image.\n\n"
        "This algorithm implements the so called 'Foerstner Corner Detector' to measure the 'cornerness' of each pixel in the image,"
        " according to [W. Foerstner: \"A feature based correspondence algorithms for image matching\", Intl."
        " Arch. Photogrammetry and Remote Sensing, vol. 24, pp 160-166, 1986]."
        " It is also known as the \"Plessey Detector\" by Harris. However, it should not be confused with the cornerResponseFunction_"
        ", another detector invented by Harris.\n\n"
        "For details see foerstnerCornerDetector_ in the vigra C++ documentation.");

    def("rohrCornerDetector2D",
        registerConverters(&pythonRohrCornerDetector2D<float, float>),               // also multiband
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in an 2D image.\n\n"
        "This algorithm implements yet another structure tensor-based corner detector,"
        " according to [K. Rohr:"
        "\"Untersuchung von grauwertabhaengigen Transformationen zur Ermittlung der optischen Flusses in Bildfolgen\""
        ", Diploma thesis, Inst. fuer Nachrichtensysteme, Univ. Karlsruhe, 1987, see also K. Rohr: "
        "\"Modelling and Identification of Characteristic Intensity Variations\""
        ", Image and Vision Computing 10:2 (1992) 66-76 and K. Rohr: "
        "\"Localization Properties of Direct Corner Detectors\", J. of Mathematical Imaging and Vision 4:2 (1994) 139-150].\n\n"
        "For details see rohrCornerDetector_ in the vigra C++ documentation.");

    def("beaudetCornerDetector2D",
        registerConverters(&pythonBeaudetCornerDetector2D<float, float>),               // also multiband
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find corners in an 2D image.\n\n"
        "This algorithm implements a corner detector according to [P.R. Beaudet:  \"Rotationally Invariant Image Operators\", Proc. Intl. Joint Conf. on Pattern Recognition, Kyoto, Japan, 1978, pp. 579-583].\n\n"
        "For details look into beaudetCornerDetector_ in vigra.");

    def("cannyEdgeImage",
        registerConverters(&pythonCannyEdgeImage<float,float>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),arg("out")=python::object()),
        "Detect and mark edges in an edge image using Canny's algorithm.\n\n"
        "For details see cannyEdgeImage_ in the vigra C++ documentation.");

    def("cannyEdgeImageWithThinning",
        registerConverters(&pythonCannyEdgeImageWithThinning<float,float>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),
        arg("addBorder")=true,arg("out")=python::object()),
        "Detect and mark edges in an edge image using Canny's algorithm.\n\n"
        "For details see cannyEdgeImageWithThinning_ in the vigra C++ documentation.");

    def("shenCastanEdgeImage",
        registerConverters(&pythonShenCastanEdgeImage<float,float>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),arg("out")=python::object()),
        "Detect and mark edges in an edge image using the Shen/Castan zero-crossing detector.\n\n"
        "For details see differenceOfExponentialEdgeImage_ in the vigra C++ documentation.");

    def("shenCastanCrackEdgeImage",
        registerConverters(&pythonShenCastanCrackEdgeImage<float,float>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),arg("out")=python::object()),
        "Detect and mark edges in a crack edge image using the Shen/Castan zero-crossing detector.\n\n"
        "For details see differenceOfExponentialCrackEdgeImage_ in the vigra C++ documentation.");

    def("removeShortEdges",
        registerConverters(&pythonRemoveShortEdges<Int32,Int32>),
        (arg("image"), arg("minEdgeLength"), arg("nonEdgeMarker"),arg("out")=python::object()),
        "Remove short edges from an edge image.\n\n"
        "For details see removeShortEdges_ in the vigra C++ documentation.");

    def("beautifyCrackEdgeImage",
        registerConverters(&pythonBeautifyCrackEdgeImage<Int32>),
        (arg("image"), arg("edgeMarker"), arg("backgroundMarker"),arg("out")=python::object()),
        "Beautify crack edge image for visualization.\n\n"
        "For details see beautifyCrackEdgeImage_ in the vigra C++ documentation.");

    def("closeGapsInCrackEdgeImage",
        registerConverters(&pythonCloseGapsInCrackEdgeImage<Int32>),
        (arg("image"), arg("edgeMarker"),arg("out")=python::object()),
        "Close one-pixel wide gaps in a cell grid edge image.\n\n"
        "For details see closeGapsInCrackEdgeImage_ in the vigra C++ documentation.");


    def("regionImageToEdgeImage",
        registerConverters(&pythonRegionImageToEdgeImage2D<Int32>),
        (arg("image"), 
         arg("edgeLabel") = 0,
         arg("out")=python::object()),
        "Transform a labeled image into an edge image.\n\n"
        "For details see regionImageToEdgeImage_ in the vigra C++ documentation.");

    def("regionImageToEdgeImage",
        registerConverters(&pythonRegionImageToEdgeImage2D<Int64>),
        (arg("image"), 
         arg("edgeLabel") = 0,
         arg("out")=python::object()));

    def("regionImageToCrackEdgeImage",
         registerConverters(&pythonRegionImageToCrackEdgeImage2D<Int64>),
         (arg("image"), 
          arg("edgeLabel") = 0, 
          arg("out")=python::object()));

    def("regionImageToCrackEdgeImage",
         registerConverters(&pythonRegionImageToCrackEdgeImage2D<Int32>),
         (arg("image"), 
          arg("edgeLabel") = 0, 
          arg("out")=python::object()),
         "Transform a labeled image into a crack edge image. \n\n"
         "For details see regionImageToCrackEdgeImage_ in the vigra C++ documentation.");

}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(edgedetection)
{
    import_vigranumpy();
    defineEdgedetection();
}
