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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyanalysis_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/edgedetection.hxx>

namespace python = boost::python;

namespace vigra
{

double Edgel__getitem__(Edgel const & e, unsigned int i)
{
    if(i > 1)
    {
        PyErr_SetString(PyExc_IndexError,
            "Edgel.__getitem__(): index out of bounds.");
        python::throw_error_already_set();
    }
    return i == 0 ? e.x : e.y;
}

void Edgel__setitem__(Edgel & e, unsigned int i, double v)
{
    if(i > 1)
    {
        PyErr_SetString(PyExc_IndexError,
            "Edgel.__setitem__(): index out of bounds.");
        python::throw_error_already_set();
    }
    if(i==0)
        e.x = Edgel::value_type(v);
    else
        e.y = Edgel::value_type(v);
}

unsigned int Edgel__len__(Edgel const & e)
{
    return 2;
}

PyObject * Edgel__repr__(Edgel const & e)
{
        std::stringstream s;
        s << std::setprecision(14)
          << "Edgel(x=" << e.x << ", y=" << e.y << ", strength=" << e.strength << ", angle=" << e.orientation << ")";
        return PyString_FromString(s.str().c_str());
}

template < class PixelType>
python::list
pythonFindEdgelsFromGrad(NumpyArray<2, TinyVector<PixelType, 2> > grad,
                         double threshold) 
{
    std::vector<Edgel> edgels;
    cannyEdgelList(srcImageRange(grad), edgels);

    python::list pyEdgels;
    for(unsigned int i = 0; i < edgels.size(); ++i)
    {
        if(edgels[i].strength >= threshold)
            pyEdgels.append(edgels[i]);
    }
    return pyEdgels;
}

template < class PixelType>
python::list
pythonFindEdgels(NumpyArray<2, PixelType> image,
                 double scale, double threshold)
{
    std::vector<Edgel> edgels;
    cannyEdgelList(srcImageRange(image), edgels, scale);

    python::list pyEdgels;
    for(unsigned int i = 0; i < edgels.size(); ++i)
    {
        if(edgels[i].strength >= threshold)
            pyEdgels.append(edgels[i]);
    }
    return pyEdgels;
}

template < class PixelType>
python::list
pythonFindEdgels3x3FromGrad(NumpyArray<2, TinyVector<PixelType, 2> > grad,
                            double threshold) 
{
    std::vector<Edgel> edgels;
    cannyEdgelList3x3(srcImageRange(grad), edgels);

    python::list pyEdgels;
    for(unsigned int i = 0; i < edgels.size(); ++i)
    {
        if(edgels[i].strength >= threshold)
            pyEdgels.append(edgels[i]);
    }
    return pyEdgels;
}

template < class PixelType>
python::list
pythonFindEdgels3x3(NumpyArray<2, PixelType> image,
                    double scale, double threshold)
{
    std::vector<Edgel> edgels;
    cannyEdgelList3x3(srcImageRange(image), edgels, scale);

    python::list pyEdgels;
    for(unsigned int i = 0; i < edgels.size(); ++i)
    {
        if(edgels[i].strength >= threshold)
            pyEdgels.append(edgels[i]);
    }
    return pyEdgels;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray 
pythonCannyEdgeImage(NumpyArray<2, Singleband<SrcPixelType> > image,
                     double scale, double threshold, DestPixelType edgeMarker, 
                     NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "cannyEdgeImage(): Output array has wrong shape.");    
    
    cannyEdgeImage(srcImageRange(image), destImage(res), 
                   scale, threshold, edgeMarker);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray 
pythonCannyEdgeImageWithThinning(NumpyArray<2, Singleband<SrcPixelType> > image,
                                 double scale, double threshold, DestPixelType edgeMarker, bool addBorder = true,
                                 NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "cannyEdgeImageWithThinning(): Output array has wrong shape.");    
    
    cannyEdgeImageWithThinning(srcImageRange(image), destImage(res),
                               scale, threshold, edgeMarker, addBorder);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray 
pythonShenCastanEdgeImage(NumpyArray<2, Singleband<SrcPixelType> > image,
                          double scale, double threshold, DestPixelType edgeMarker, 
                          NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "shenCastanEdgeImage(): Output array has wrong shape.");    
    
    differenceOfExponentialEdgeImage(srcImageRange(image), destImage(res), 
                                     scale, threshold, edgeMarker);
     
    return res;
}

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray 
pythonShenCastanCrackEdgeImage(NumpyArray<2, Singleband<SrcPixelType> > image,
                               double scale, double threshold, DestPixelType edgeMarker, 
                               NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(2*image.shape() - MultiArrayShape<2>::type(1,1), 
                       "shenCastanCrackEdgeImage(): Output array has wrong shape.");    
    
    differenceOfExponentialCrackEdgeImage(srcImageRange(image), destImage(res), 
                                          scale, threshold, edgeMarker);
     
    return res;
}

template < class PixelType>
NumpyAnyArray 
pythonRemoveShortEdges(NumpyArray<2, Singleband<PixelType> > image,
                       int minEdgeLength, PixelType nonEdgeMarker, 
                       NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "removeShortEdges(): Output array has wrong shape.");    
    
    copyImage(srcImageRange(image), destImage(res));
    removeShortEdges(destImageRange(res), minEdgeLength, nonEdgeMarker);
     
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonBeautifyCrackEdgeImage(NumpyArray<2, Singleband<PixelType> > image,
                             PixelType edgeMarker, 
                             PixelType backgroundMarker,
                             NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "beautifyCrackEdgeImage(): Output array has wrong shape.");    
    
    copyImage(srcImageRange(image), destImage(res));
    beautifyCrackEdgeImage(destImageRange(res), edgeMarker, backgroundMarker);
     
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonCloseGapsInCrackEdgeImage(NumpyArray<2, Singleband<PixelType> > image,
                                PixelType edgeMarker,
                                NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "closeGapsInCrackEdgeImage(): Output array has wrong shape.");    
    
    copyImage(srcImageRange(image), destImage(res));
    closeGapsInCrackEdgeImage(destImageRange(res), edgeMarker);
     
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonRegionImageToCrackEdgeImage(NumpyArray<2, Singleband<PixelType> > image,
                                  PixelType edgeLabel = 0,
                                  NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(2*image.shape() - MultiArrayShape<2>::type(1,1), 
                 "regionImageToCrackEdgeImage(): Output array has wrong shape. Needs to be (w,h)*2 -1");

    regionImageToCrackEdgeImage(srcImageRange(image), destImage(res), edgeLabel);
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonRegionImageToEdgeImage(NumpyArray<2, Singleband<PixelType> > image,
                             PixelType edgeLabel = 0,
                             NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "regionImageToEdgeImage2D(): Output array has wrong shape.");

    regionImageToEdgeImage(srcImageRange(image), destImage(res), edgeLabel);
    return res;
}


void defineEdgedetection()
{
    using namespace python;

    class_<Edgel> edgel("Edgel", "Represent an Edgel at a particular subpixel position (x, y), having "
                                  "given 'strength' and 'orientation'.\n\n"
                                  "For details, see Edgel_ in the vigra C++ documentation.\n",
                         init<>("Standard constructor::\n\n   Edgel()\n\n"));
    edgel
       .def(init<float, float, float, float>(args("x", "y", "strength", "orientation"), 
            "Constructor::\n\n    Edgel(x, y, strength, orientation)\n\n"))
       .def_readwrite("x", &Edgel::x, "The edgel's x position.")
       .def_readwrite("y", &Edgel::y, "The edgel's y position.")
       .def_readwrite("strength", &Edgel::strength, "The edgel's strength.")
       .def_readwrite("orientation", &Edgel::orientation, "The edgel's orientation.")
       .def("__getitem__", &Edgel__getitem__)
       .def("__setitem__", &Edgel__setitem__)
       .def("__repr__", &Edgel__repr__)
       .def("__len__", &Edgel__len__)
       ;

    def("cannyEdgelList", 
        registerConverters(&pythonFindEdgelsFromGrad<float>),
        args("gradient", "threshold"),
        "Return a list of :class:`Edgel` objects whose strength is at least 'threshold'.\n\n"
        "The function comes in two forms::\n\n"
        "    cannyEdgelList(gradient, threshold) -> list\n"
        "    cannyEdgelList(image, scale, threshold) -> list\n\n"
        "The first form expects a gradient image (i.e. with two channels) to compute "
        "edgels, whereas the second form expects a scalar image and computes the "
        "gradient internally at 'scale'.\n\n"
        "For details see cannyEdgelList_ in the vigra C++ documentation.\n");

    def("cannyEdgelList",  
        registerConverters(&pythonFindEdgels<float>),
        args("image", "scale", "threshold"));

    def("cannyEdgelList3x3", 
        registerConverters(&pythonFindEdgels3x3FromGrad<float>),
        args("gradient", "threshold"),
        "Return a list of :class:`Edgel` objects whose strength is at least 'threshold'.\n\n"
        "The function comes in two forms::\n\n"
        "    cannyEdgelList3x3(gradient, threshold) -> list\n"
        "    cannyEdgelList3x3(image, scale, threshold) -> list\n\n"
        "The first form expects a gradient image (i.e. with two channels) to compute "
        "edgels, whereas the second form expects a scalar image and computes the "
        "gradient internally at 'scale'. The results are slightly better than "
        "those of :func:`cannyEdgelList`.\n\n"
        "For details see cannyEdgelList3x3_ in the vigra C++ documentation.\n");

    def("cannyEdgelList3x3",  
        registerConverters(&pythonFindEdgels3x3<float>),
        args("image", "scale", "threshold"));

    def("cannyEdgeImage",
        registerConverters(&pythonCannyEdgeImage<float, UInt8>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),arg("out")=python::object()),
        "Detect and mark edges in an edge image using Canny's algorithm.\n\n"
        "For details see cannyEdgeImage_ in the vigra C++ documentation.\n");

    def("cannyEdgeImageWithThinning",
        registerConverters(&pythonCannyEdgeImageWithThinning<float, UInt8>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),
        arg("addBorder")=true,arg("out")=python::object()),
        "Detect and mark edges in an edge image using Canny's algorithm.\n\n"
        "For details see cannyEdgeImageWithThinning_ in the vigra C++ documentation.\n");

    def("shenCastanEdgeImage",
        registerConverters(&pythonShenCastanEdgeImage<float, UInt8>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),arg("out")=python::object()),
        "Detect and mark edges in an edge image using the Shen/Castan zero-crossing detector.\n\n"
        "For details see differenceOfExponentialEdgeImage_ in the vigra C++ documentation.\n");

    def("shenCastanCrackEdgeImage",
        registerConverters(&pythonShenCastanCrackEdgeImage<float, UInt8>),
        (arg("image"), arg("scale"), arg("threshold"), arg("edgeMarker"),arg("out")=python::object()),
        "Detect and mark edges in a crack edge image using the Shen/Castan zero-crossing detector.\n\n"
        "For details see differenceOfExponentialCrackEdgeImage_ in the vigra C++ documentation.\n");

    def("removeShortEdges",
        registerConverters(&pythonRemoveShortEdges<UInt8>),
        (arg("image"), arg("minEdgeLength"), arg("nonEdgeMarker"),arg("out")=python::object()),
        "Remove short edges from an edge image.\n\n"
        "For details see removeShortEdges_ in the vigra C++ documentation.\n");

    def("beautifyCrackEdgeImage",
        registerConverters(&pythonBeautifyCrackEdgeImage<UInt8>),
        (arg("image"), arg("edgeMarker"), arg("backgroundMarker"),arg("out")=python::object()),
        "Beautify crack edge image for visualization.\n\n"
        "For details see beautifyCrackEdgeImage_ in the vigra C++ documentation.\n");

    def("closeGapsInCrackEdgeImage",
        registerConverters(&pythonCloseGapsInCrackEdgeImage<UInt8>),
        (arg("image"), arg("edgeMarker"),arg("out")=python::object()),
        "Close one-pixel wide gaps in a cell grid edge image.\n\n"
        "For details see closeGapsInCrackEdgeImage_ in the vigra C++ documentation.\n");

    def("regionImageToEdgeImage",
        registerConverters(&pythonRegionImageToEdgeImage<npy_uint32>),
        (arg("image"), 
         arg("edgeLabel") = 0,
         arg("out")=python::object()),
        "Transform a labeled image into an edge image.\n\n"
        "For details see regionImageToEdgeImage_ in the vigra C++ documentation.\n");

    def("regionImageToEdgeImage",
        registerConverters(&pythonRegionImageToEdgeImage<npy_uint64>),
        (arg("image"), 
         arg("edgeLabel") = 0,
         arg("out")=python::object()));

    def("regionImageToCrackEdgeImage",
         registerConverters(&pythonRegionImageToCrackEdgeImage<npy_uint32>),
         (arg("image"), 
          arg("edgeLabel") = 0, 
          arg("out")=python::object()),
         "Transform a labeled image into a crack edge image. \n\n"
         "For details see regionImageToCrackEdgeImage_ in the vigra C++ documentation.\n");

    def("regionImageToCrackEdgeImage",
         registerConverters(&pythonRegionImageToCrackEdgeImage<npy_uint64>),
         (arg("image"), 
          arg("edgeLabel") = 0, 
          arg("out")=python::object()));

}

} // namespace vigra
