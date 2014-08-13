/*
 * morpho.cxx
 *
 *  Created on: Aug 21, 2013
 *      Author: Thomas Walter
 */


// define PY_ARRAY_UNIQUE_SYMBOL (required by the numpy C-API)
#define PY_ARRAY_UNIQUE_SYMBOL vigranumpymorpho_PyArray_API

// include the vigranumpy C++ API
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

// specific includes
#include <vigra/morpho_basic.hxx>
#include <vigra/morpho_utilities.hxx>
#include <vigra/morpho_geodesy.hxx>
#include <vigra/morpho_criteria.hxx>
#include <vigra/morpho_watershed.hxx>
#include <vigra/morpho_dynamic.hxx>

// implementation of your wrapper functions and classes
//template <class PixelType>
//NumpyAnyArray
//pythonConvolveImage(NumpyArray<3, Multiband<PixelType> > image,
//                    TwoDKernel const & kernel,
//                    NumpyArray<3, Multiband<PixelType> > res = python::object())
//{
//    res.reshapeIfEmpty(image.taggedShape(),
//            "convolve(): Output array has wrong shape.");
//
//    {
//        PyAllowThreads _pythread;
//        for(int k=0;k<image.shape(2);++k)
//        {
//            MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
//            MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
//            convolveImage(srcImageRange(bimage), destImage(bres),
//                          kernel2d(kernel));
//        }
//    }
//    return res;
//}

namespace python = boost::python;

namespace vigra
{

//template<class T>
//void pythonInitExplicitlyKernel1D(Kernel1D<T> &k, int left, int right, NumpyArray<1,T> contents)
//{
//    vigra_precondition(contents.size() == 1 || right-left+1 == contents.shape(0),
//              "Kernel1D::initExplicitly(): 'contents' must contain as many elements as the kernel (or just one element).");
//
//    k.initExplicitly(left,right);
//    for(int i=left; i<=right; ++i)
//    {
//        k[i] = (contents.size() == 1)
//                     ? contents(0)
//                     : contents(i-left);
//    }
//}

//NumpyArray<2, Singleband<PixelType> > image,

class structuringElement2D_Conv: public morpho::structuringElement2D
{
    public:
        structuringElement2D_Conv(): structuringElement2D() {}
        structuringElement2D_Conv(boost::python::list a, int s=1) {
            size = s;
            boost::python::ssize_t len = boost::python::len(a);
            for(int i=0; i<len;i++){
                Diff2D from_python((int)boost::python::extract<int>(a[i][0]),
                                   (int)boost::python::extract<int>(a[i][1]));
                support.push_back(from_python);

            }
            CalculateExtension();
        }
};


template <class PixelType>
NumpyAnyArray
pythonMorphoReconstructionByDilation(NumpyArray<2, Singleband<PixelType> > marker,
                                     NumpyArray<2, Singleband<PixelType> > mask,
                                     int neighborhood_graph
                                   )
{
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(marker.width(), marker.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoReconstructionByDilation(destImageRange(marker), srcImage(mask), nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoReconstructionByDilation(destImageRange(marker), srcImage(mask), nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoOpeningByReconstruction(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return marker;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoReconstructionByErosion(NumpyArray<2, Singleband<PixelType> > marker,
                                     NumpyArray<2, Singleband<PixelType> > mask,
                                     int neighborhood_graph
                                   )
{
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(marker.width(), marker.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoReconstructionByErosion(destImageRange(marker), srcImage(mask), nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoReconstructionByErosion(destImageRange(marker), srcImage(mask), nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoOpeningByReconstruction(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return marker;
}


template <class PixelType>
NumpyAnyArray
pythonMorphoOpeningByReconstruction(NumpyArray<2, Singleband<PixelType> > image,
                                    structuringElement2D_Conv se,
                                    int neighborhood_graph,
                                    NumpyArray<2, Singleband<PixelType> > res
                                   )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoOpeningByReconstruction(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoOpeningByReconstruction(image, res, se, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoOpeningByReconstruction(image, res, se, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoOpeningByReconstruction(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            //std::cout << "This neighborhood is not defined." << std::endl;
            break;
        }
    };

    //morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
    //morpho::morphoOpeningByReconstruction(image, res, se, nb);

    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoDiameterOpening(NumpyArray<2, Singleband<PixelType> > image,
                            int diameter,
                            int neighborhood_graph,
                            NumpyArray<2, Singleband<PixelType> > res
                            )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoDiameterOpening(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoDiameterOpening(image, res, diameter, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoDiameterOpening(image, res, diameter, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoDiameterOpening(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}

template <class PixelType, class LabelType>
NumpyAnyArray
pythonMorphoWatershed(NumpyArray<2, Singleband<PixelType> > image,
                      int neighborhood_graph,
                      NumpyArray<2, Singleband<LabelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoWatershed(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoWatershed(image, res, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoWatershed(image, res, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoWatershed(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}


template <class PixelType, class LabelType>
NumpyAnyArray
pythonMorphoSelectiveWatershed(NumpyArray<2, Singleband<PixelType> > image,
                               int dyn_thresh,
                               int neighborhood_graph,
                               NumpyArray<2, Singleband<LabelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoSelectiveWatershed(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoSelectiveWatershed(image, res, dyn_thresh, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoSelectiveWatershed(image, res, dyn_thresh, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoSelectiveWatershed(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoDynMinima(NumpyArray<2, Singleband<PixelType> > image,
                      int neighborhood_graph,
                      NumpyArray<2, Singleband<PixelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
                       "pythonMorphoDynMinima(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());
//
//    def("morphoDynMinima", registerConverters(&pythonMorphoDynMinima<UInt8>),
//        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
//        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
//        "(default is 8 neighborhood) as inputs. "
//        "The function calculates the morphological dynamics for each minimum"
//        "and returns an image which is only non-zero on the local minima and ");
//    def("morphoDynMinima", registerConverters(&pythonMorphoDynMinima<UInt16>),
//        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
//        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
//        "(default is 8 neighborhood) as inputs. "
//        "The function calculates the morphological dynamics for each minimum"
//        "and returns an image which is only non-zero on the local minima and ");

    //std::vector<PixelType> dynamics;

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoDynMinima(image, res, nb);
            //vigra::NumpyArray<2, vigra::Singleband<unsigned short>, vigra::StridedArrayTag>::value_type
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            //morpho::morphoDynMinima(srcImageRange(image), destImageRange(res), nb);
            morpho::morphoDynMinima(image, res, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoDynMinima(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}

template <class PixelType, class LabelType>
NumpyAnyArray
pythonMorphoBifurcationPoints(NumpyArray<2, Singleband<PixelType> > image,
                              int neighborhood_graph,
                              NumpyArray<2, Singleband<LabelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
                       "morphoBifurcationPoints(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoBifurcationPoints(image, res, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoBifurcationPoints(image, res, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoBifurcationPoints(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}

template <class PixelType, class MarkerType, class LabelType>
NumpyAnyArray
pythonMorphoConstrainedWatershed(NumpyArray<2, Singleband<PixelType> > image,
                                 NumpyArray<2, Singleband<MarkerType> > marker_image,
                                 int neighborhood_graph,
                                 NumpyArray<2, Singleband<LabelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoConstrainedWatershed(): Output array has wrong shape.");

    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoConstrainedWatershed(image, marker_image, res, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoConstrainedWatershed(image, marker_image, res, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoConstrainedWatershed(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}


template <class PixelType>
NumpyAnyArray
pythonMorphoDiameterClosing(NumpyArray<2, Singleband<PixelType> > image,
                            int diameter,
                            int neighborhood_graph,
                            NumpyArray<2, Singleband<PixelType> > res
                            )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoDiameterClosing(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoDiameterClosing(image, res, diameter, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoDiameterClosing(image, res, diameter, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoDiameterClosing(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoAreaOpening(NumpyArray<2, Singleband<PixelType> > image,
                                    int area,
                                    int neighborhood_graph,
                                    NumpyArray<2, Singleband<PixelType> > res
                                   )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoAreaOpening(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoAreaOpening(image, res, area, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoAreaOpening(image, res, area, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoAreaOpening(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            //std::cout << "This neighborhood is not defined." << std::endl;
            break;
        }
    };

    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoAreaClosing(NumpyArray<2, Singleband<PixelType> > image,
                                    int area,
                                    int neighborhood_graph,
                                    NumpyArray<2, Singleband<PixelType> > res
                                   )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoAreaClosing(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoAreaClosing(image, res, area, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoAreaClosing(image, res, area, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoAreaClosing(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            //std::cout << "This neighborhood is not defined." << std::endl;
            break;
        }
    };

    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoClosingByReconstruction(NumpyArray<2, Singleband<PixelType> > image,
                                    structuringElement2D_Conv se,
                                    int neighborhood_graph,
                                    NumpyArray<2, Singleband<PixelType> > res
                                   )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoClosingByReconstruction(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoClosingByReconstruction(image, res, se, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoClosingByReconstruction(image, res, se, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoClosingByReconstruction(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}



template <class PixelType>
NumpyAnyArray
pythonMorphoOpeningByDynamics(NumpyArray<2, Singleband<PixelType> > image,
                              int h,
                              int neighborhood_graph,
                              NumpyArray<2, Singleband<PixelType> > res
)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoClosingByReconstruction(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoOpeningByDynamics(image, res, h, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoOpeningByDynamics(image, res, h, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoOpeningByDynamics(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}


template <class PixelType>
NumpyAnyArray
pythonMorphoClosingByDynamics(NumpyArray<2, Singleband<PixelType> > image,
                              int h,
                              int neighborhood_graph,
                              NumpyArray<2, Singleband<PixelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoClosingByReconstruction(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    Diff2D image_size = Diff2D(image.width(), image.height());

    switch(neighborhood_graph) {
        case 8:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER8, image_size);
            morpho::morphoClosingByDynamics(image, res, h, nb);
            break;
        }
        case 4:
        {
            morpho::neighborhood2D nb(morpho::WITHOUTCENTER4, image_size);
            morpho::morphoClosingByDynamics(image, res, h, nb);
            break;
        }
        default:
        {
            vigra_precondition(false,
                         "pythonMorphoClosingByDynamics(): Unknown neighborhood graph.\n"
                         "This function is only defined for 4- and 8-neighborhoods.");
            break;
        }
    };

    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoDilation(NumpyArray<2, Singleband<PixelType> > image,
                    structuringElement2D_Conv se,
                    NumpyArray<2, Singleband<PixelType> > res // = python::object(),
                    )
{

    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoDilation(srcImageRange(image), destImageRange(res), se);
    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoOpening(NumpyArray<2, Singleband<PixelType> > image,
                    structuringElement2D_Conv se,
                    NumpyArray<2, Singleband<PixelType> > res // = python::object(),
                    )
{

    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoOpening(srcImageRange(image), destImageRange(res), se);
    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoClosing(NumpyArray<2, Singleband<PixelType> > image,
                    structuringElement2D_Conv se,
                    NumpyArray<2, Singleband<PixelType> > res // = python::object(),
                    )
{

    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoClosing(srcImageRange(image), destImageRange(res), se);
    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoInternalGradient(NumpyArray<2, Singleband<PixelType> > image,
                             structuringElement2D_Conv se,
                             NumpyArray<2, Singleband<PixelType> > res // = python::object(),
)
{

    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoInternalGradient(srcImageRange(image), destImageRange(res), se);
    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoExternalGradient(NumpyArray<2, Singleband<PixelType> > image,
                             structuringElement2D_Conv se,
                    NumpyArray<2, Singleband<PixelType> > res // = python::object(),
                    )
{

    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoExternalGradient(srcImageRange(image), destImageRange(res), se);
    return res;
}

template <class PixelType>
NumpyAnyArray
pythonMorphoGradient(NumpyArray<2, Singleband<PixelType> > image,
                     structuringElement2D_Conv se,
                     NumpyArray<2, Singleband<PixelType> > res
                    )
{

    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoGradient(srcImageRange(image), destImageRange(res), se);
    return res;
}


template <class PixelType>
NumpyAnyArray
pythonMorphoErosion(NumpyArray<2, Singleband<PixelType> > image,
                    structuringElement2D_Conv se,
                    NumpyArray<2, Singleband<PixelType> > res // = python::object(),
                    )
{
    res.reshapeIfEmpty(image.taggedShape(),
            "morphoErosion(): Output array has wrong shape.");
    PyAllowThreads _pythread;

    morpho::morphoErosion(srcImageRange(image), destImageRange(res), se);
    return res;
}

//struct Diff2Dlist_from_python_pointlist
//{
//
//  // Determine if obj_ptr can be converted in a QString
//  static void* convertible(PyObject* obj_ptr)
//    {
//      //if (!PyString_Check(obj_ptr)) return 0;
//      //return obj_ptr;
//    }
//
//  // Convert obj_ptr into a list of points
//   static void construct(
//     PyObject* obj_ptr,
//     boost::python::converter::rvalue_from_python_stage1_data* data)
//     {
//       // Extract the character data from the python string
//       const char* value = PyString_AsString(obj_ptr);
//
//       // Verify that obj_ptr is a string (should be ensured by convertible())
//       assert(value);
//
//       // Grab pointer to memory into which to construct the new QString
//       void* storage = (
//         (boost::python::converter::rvalue_from_python_storage<QString>*)
//         data)->storage.bytes;
//
//       // in-place construct the new QString using the character data
//       // extraced from the python object
//       new (storage) QString(value);
//
//       // Stash the memory chunk pointer for later use by boost.python
//       data->convertible = storage;
//     }
//
//};


//class VecTest {
//
//    protected:
//        std::vector<int> myvec;
//
//    public:
//        //VecTest(std::vector<int> a): myvec(a) {}
//        VecTest(boost::python::list a)
//        {
//
//            boost::python::ssize_t len = boost::python::len(a);
//            for(int i=0; i<len;i++){
//                myvec.push_back((int)boost::python::extract<int>(a[i]));
//            }
//        }
//
//        void output() {
//            for(std::vector<int>::iterator iter = myvec.begin();
//              iter != myvec.end();
//              ++iter)
//              {
//                std::cout << (*iter) << std::endl;
//              }
//            std::cout << std::endl;
//        };
//};
//
//class PointList {
//
//    protected:
//        std::vector<int> x;
//        std::vector<int> y;
//
//    public:
//        //VecTest(std::vector<int> a): myvec(a) {}
//        PointList(boost::python::list a)
//        {
//
//            boost::python::ssize_t len = boost::python::len(a);
//            std::cout << "length of list: " << len << std::endl;
//            for(int i=0; i<len;i++){
//                x.push_back((int)boost::python::extract<int>(a[i][0]));
//            }
//        }
//
//        void output() {
//            for(std::vector<int>::iterator iter = x.begin();
//              iter != x.end();
//              ++iter)
//              {
//                std::cout << (*iter) << std::endl;
//              }
//            std::cout << std::endl;
//        };
//};
//

using namespace boost::python;

// the argument of the init macro must be the module name
BOOST_PYTHON_MODULE_INIT(morpho)
{
    // initialize numpy and vigranumpy
    vigra::import_vigranumpy();

    // Structuring elements
    //class_<morpho::structuringElement2D>("structuringElement2D", init<std::vector<vigra::Diff2D> >())
    //    .def("output", &morpho::structuringElement2D::output);
//    class_<std::vector<int> >("VectorOfInt")
//            .def(vector_indexing_suite<std::vector<int> >() )
//        ;

    //typedef std::vector<int> VectorOfInt;

//    class_<VecTest>("VecTest", init<boost::python::list>())
//        .def("output", &VecTest::output)
//        ;
//
//    class_<PointList>("PointList", init<boost::python::list>())
//        .def("output", &PointList::output)
//        ;

//    class_<TCurrency>( "TCurrency" )
//        .def( init<long> )
//        .def( init<const std::string&> )
//        ...
//        ;
//
    class_<structuringElement2D_Conv>("structuringElement2D",
                                      "\nA Structuring Element in two dimensions for morphological operations. "
                                      "\nThe constructor takes a a list of points (python list of lists)."
                                      "\nThe structuring element is needed for morphological operations such as erosion, dilation"
                                      "opening, closing. "
                                      "\n\nExample: simple 4 neighborhood structuring element:"
                                      "\n>>> se = vigra.morpho.structuringElement2D([[0,0], [1,0], [0,1], [-1, 0], [0, -1]])"
                                      "\n\nExample: simple 8 neighborhood structuring element:"
                                      "\n>>> se = vigra.morpho.structuringElement2D([[0,0], [1,0], [0,1], [-1, 0], [0, -1], [1, 1], [1, -1], [-1, -1], [-1, 1]])"
                                      )
        .def(init<boost::python::list>())//("structuringElement2D", init<boost::python::list>())
        .def(init<boost::python::list, int>())
        .def("output", &structuringElement2D_Conv::output)
        .def("numberOfPixels", &structuringElement2D_Conv::numberOfPixels)
        .def_readwrite("size", &structuringElement2D_Conv::size)
        .def("transpose", &structuringElement2D_Conv::transpose)
        ;

    // export morphological functions
//    def("morphoDilation", registerConverters(&pythonMorphoDilation<UInt8>),
//        (arg("image"), arg("size"), arg("res")=object()),
//        "Morphological Dilation of the image with a structuring element."
//        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)"
//        "At the moment we only accept size as a parameter, and the structuring"
//        "element is a square with radius of this size "
//        "(i.e. with each side of length 2*size + 1)");

    def("morphoDilation", registerConverters(&pythonMorphoDilation<UInt8>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Dilation of the image with a structuring element."
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)"
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoDilation", registerConverters(&pythonMorphoDilation<UInt16>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Dilation of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoErosion", registerConverters(&pythonMorphoErosion<UInt8>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Erosion of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoErosion", registerConverters(&pythonMorphoErosion<UInt16>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Erosion of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoOpening", registerConverters(&pythonMorphoOpening<UInt8>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Opening of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoOpening", registerConverters(&pythonMorphoOpening<UInt16>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Opening of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoClosing", registerConverters(&pythonMorphoClosing<UInt8>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Closing of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoClosing", registerConverters(&pythonMorphoClosing<UInt16>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Closing of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoInternalGradient", registerConverters(&pythonMorphoInternalGradient<UInt8>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Internal gradient (f - ero(f)) of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoInternalGradient", registerConverters(&pythonMorphoInternalGradient<UInt16>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological Internal gradient (f - ero(f)) of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoExternalGradient", registerConverters(&pythonMorphoExternalGradient<UInt8>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological external gradient (dil(f) - f) of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoExternalGradient", registerConverters(&pythonMorphoExternalGradient<UInt16>),
        (arg("image"), arg("se"), arg("res")=object()),
        "Morphological external gradient (dil(f) - f) of the image with a structuring element."
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoOpeningByReconstruction", registerConverters(&pythonMorphoOpeningByReconstruction<UInt8>),
        (arg("image"), arg("se"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Morphological opening, followed by a reconstruction of the opened image under the original image."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)"
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoOpeningByReconstruction", registerConverters(&pythonMorphoOpeningByReconstruction<UInt16>),
        (arg("image"), arg("se"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Morphological opening, followed by a reconstruction of the opened image under the original image."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)"
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoClosingByReconstruction", registerConverters(&pythonMorphoClosingByReconstruction<UInt8>),
        (arg("image"), arg("se"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Morphological closing, followed by a reconstruction of the closed image over the original image."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)"
        "The structuring element is of type vigra.morpho.structuringElement2D");
    def("morphoClosingByReconstruction", registerConverters(&pythonMorphoClosingByReconstruction<UInt16>),
        (arg("image"), arg("se"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Morphological closing, followed by a reconstruction of the closed image over the original image."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)"
        "The structuring element is of type vigra.morpho.structuringElement2D");

    def("morphoClosingByDynamics", registerConverters(&pythonMorphoClosingByDynamics<UInt8>),
        (arg("image"), arg("h"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Let f be an image and h a positive constant. This function calculates then the "
        "morphological reconstruction of f + h over f."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)");
    def("morphoClosingByDynamics", registerConverters(&pythonMorphoClosingByDynamics<UInt16>),
        (arg("image"), arg("h"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Let f be an image and h a positive constant. This function calculates then the "
        "morphological reconstruction of f + h over f."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)");

    def("morphoOpeningByDynamics", registerConverters(&pythonMorphoOpeningByDynamics<UInt8>),
        (arg("image"), arg("h"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Let f be an image and h a positive constant. This function calculates then the "
        "morphological reconstruction of f - h under f."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)");
    def("morphoOpeningByDynamics", registerConverters(&pythonMorphoOpeningByDynamics<UInt16>),
        (arg("image"), arg("h"), arg("neighborhood_graph")=8, arg("res")=object()),
        "Let f be an image and h a positive constant. This function calculates then the "
        "morphological reconstruction of f - h under f."
        "This function is only implemented for Uint8 and Uint16 images"
        "(for other images, you can use img_mod = img.astype(np.dtype('uint16')) for conversion)");

    def("morphoReconstructionByDilation", registerConverters(&pythonMorphoReconstructionByDilation<UInt8>),
        (arg("marker"), arg("mask"), arg("neighborhood_graph")=8),
        "Takes two images and calculates the morphological reconstruction of the"
        "first under the second. Typically the first image is a marker image and the second is the original"
        "image to be analyzed.");
    def("morphoReconstructionByDilation", registerConverters(&pythonMorphoReconstructionByDilation<UInt16>),
        (arg("marker"), arg("mask"), arg("neighborhood_graph")=8),
        "Takes two images and calculates the morphological reconstruction of the"
        "first under the second. Typically the first image is a marker image and the second is the original"
        "image to be analyzed.");

    def("morphoReconstructionByErosion", registerConverters(&pythonMorphoReconstructionByErosion<UInt8>),
        (arg("marker"), arg("mask"), arg("neighborhood_graph")=8),
        "Takes two images and calculates the morphological reconstruction of the"
        "first over the second. For this, the first image has to be greater than the second (pointwise)."
        "Typically the first image is a marker image and the second is the original"
        "image to be analyzed.");
    def("morphoReconstructionByErosion", registerConverters(&pythonMorphoReconstructionByErosion<UInt16>),
        (arg("marker"), arg("mask"), arg("neighborhood_graph")=8),
        "Takes two images and calculates the morphological reconstruction of the"
        "first over the second. For this, the first image has to be greater than the second (pointwise)."
        "Typically the first image is a marker image and the second is the original"
        "image to be analyzed.");


    def("morphoAreaOpening", registerConverters(&pythonMorphoAreaOpening<UInt8>),
        (arg("image"), arg("area"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, an area parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the area opening of the input image, i.e."
        "the largest output image smaller than the input image for which any "
        "connected component resulting from a threshold operation would have at least <area> pixels.");
    def("morphoAreaOpening", registerConverters(&pythonMorphoAreaOpening<UInt16>),
        (arg("image"), arg("area"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, an area parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the area opening of the input image, i.e."
        "the largest output image smaller than the input image for which any "
        "connected component resulting from a threshold operation would have at least <area> pixels.");

    def("morphoAreaClosing", registerConverters(&pythonMorphoAreaClosing<UInt8>),
        (arg("image"), arg("area"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, an area parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the area closing of the input image, i.e."
        "the smallest output image larger than the input image for which any "
        "connected component resulting from a threshold operation would have at least <area> pixels.");
    def("morphoAreaClosing", registerConverters(&pythonMorphoAreaClosing<UInt16>),
        (arg("image"), arg("area"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, an area parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the area closing of the input image, i.e."
        "the smallest output image larger than the input image for which any "
        "connected component resulting from a threshold operation would have at least <area> pixels.");

    def("morphoWatershed", registerConverters(&pythonMorphoWatershed<UInt8, UInt16>),
        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates a watershed from the local minima.");
    def("morphoWatershed", registerConverters(&pythonMorphoWatershed<UInt16, UInt16>),
        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates a watershed from the local minima.");

    def("morphoSelectiveWatershed", registerConverters(&pythonMorphoSelectiveWatershed<UInt8, UInt16>),
        (arg("image"), arg("dyn_thresh"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image a dynamic threshold and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates a watershed from the local minima.");
    def("morphoSelectiveWatershed", registerConverters(&pythonMorphoSelectiveWatershed<UInt16, UInt16>),
        (arg("image"), arg("dyn_thresh"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image a dynamic threshold and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates a watershed from the local minima.");

    def("morphoDynMinima", registerConverters(&pythonMorphoDynMinima<UInt8>),
        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates the morphological dynamics for each minimum"
        "and returns an image which is only non-zero on the local minima and ");
    def("morphoDynMinima", registerConverters(&pythonMorphoDynMinima<UInt16>),
        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates the morphological dynamics for each minimum"
        "and returns an image which is only non-zero on the local minima and ");

    def("morphoBifurcationPoints", registerConverters(&pythonMorphoBifurcationPoints<UInt8, UInt8>),
        (arg("image"), arg("neighborhood_graph")=4, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The input image is supposed to contain a skeleton, where all non-zero pixels are taken as"
        "foreground (i.e. belonging to the pixel)."
        "The function finds the bifurcation points (all points with more than 2 neighbors).");
    def("morphoBifurcationPoints", registerConverters(&pythonMorphoBifurcationPoints<UInt16, UInt16>),
        (arg("image"), arg("neighborhood_graph")=4, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The input image is supposed to contain a skeleton, where all non-zero pixels are taken as"
        "foreground (i.e. belonging to the pixel)."
        "The function finds the bifurcation points (all points with more than 2 neighbors).");
    def("morphoBifurcationPoints", registerConverters(&pythonMorphoBifurcationPoints<float, float>),
        (arg("image"), arg("neighborhood_graph")=4, arg("res")=object()),
        "takes an image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The input image is supposed to contain a skeleton, where all non-zero pixels are taken as"
        "foreground (i.e. belonging to the pixel)."
        "The function finds the bifurcation points (all points with more than 2 neighbors).");

    def("morphoConstrainedWatershed", registerConverters(&pythonMorphoConstrainedWatershed<UInt8, UInt8, UInt16>),
        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, a marker image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates calculates a constrained watershed from the markers.");
    def("morphoConstrainedWatershed", registerConverters(&pythonMorphoConstrainedWatershed<UInt16, UInt16, UInt16>),
        (arg("image"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, a marker image and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is 8 neighborhood) as inputs. "
        "The function calculates calculates a constrained watershed from the markers.");

    def("morphoDiameterOpening", registerConverters(&pythonMorphoDiameterOpening<UInt8>),
        (arg("image"), arg("diameter"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, a diameter parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the diameter opening of the input image, i.e."
        "the largest output image smaller than the input image for which any "
        "connected component resulting from a threshold operation would have at least a diameter of <diameter>");
    def("morphoDiameterOpening", registerConverters(&pythonMorphoDiameterOpening<UInt16>),
        (arg("image"), arg("diameter"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, a diameter parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the diameter opening of the input image, i.e."
        "the largest output image smaller than the input image for which any "
        "connected component resulting from a threshold operation would have at least a diameter of <diameter>");

    def("morphoDiameterClosing", registerConverters(&pythonMorphoDiameterClosing<UInt8>),
        (arg("image"), arg("diameter"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, a diameter parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the diameter closing of the input image, i.e."
        "the smallest output image larger than the input image for which any "
        "connected component resulting from a threshold operation would have at least a diameter of <diameter>");
    def("morphoDiameterClosing", registerConverters(&pythonMorphoDiameterClosing<UInt16>),
        (arg("image"), arg("diameter"), arg("neighborhood_graph")=8, arg("res")=object()),
        "takes an image, a diameter parameter and (optionally) a constant (4 or 8) indicating the neighborhood graph"
        "(default is a 8 neighborhood) as inputs. "
        "The function calculates the diameter closing of the input image, i.e."
        "the smallest output image larger than the input image for which any "
        "connected component resulting from a threshold operation would have at least a diameter of <diameter>");


}
} // end of namespace vigra

