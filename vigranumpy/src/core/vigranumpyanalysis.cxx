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
#include <vigra/localminmax.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/symmetry.hxx>
#include <vigra/cornerdetection.hxx>
#include <vigra/edgedetection.hxx>
#include <vigra/gradient_energy_tensor.hxx>
#include <vigra/orientedtensorfilters.hxx>
#include <vigra/tensorutilities.hxx>
#include <vigra/watersheds.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/boundarytensor.hxx>

#include <cmath>

namespace python = boost::python;

namespace vigra
{

template < class PixelType >
NumpyAnyArray pythonWatersheds2D(NumpyArray<2, Singleband<PixelType> > image,
                           int neighborhood = 8,
                           NumpyArray<2, Singleband<Int32> > res=python::object())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
       "watersheds2D(image, neighborhood): neighborhood must be 4 or 8.");
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "Watersheds(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 4:
        {
            watersheds(srcImageRange(image), destImage(res),
                FourNeighborCode());
            break;
        }
        case 8:
        {
            watersheds(srcImageRange(image), destImage(res),
                EightNeighborCode());
            break;
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywatersheds2D, pythonWatersheds2D)


template < class PixelType,  typename DestPixelType >
NumpyAnyArray pythonSeededRegionGrowing2D(NumpyArray<2, Singleband<PixelType> > image,
                                    SRGType srgType=CompleteGrow, 
                                    NumpyArray<2, Singleband<DestPixelType> > res=python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "seededRegionGrowing2D(): Output array has wrong shape.");
    //NumpyArray<2, Singleband<DestPixelType> > labels(image.shape());
    extendedLocalMinima(srcImageRange(image), destImage(res));
    DestPixelType max_region_label =
        labelImageWithBackground(srcImageRange(res),
        destImage(res), false, 0);
    ArrayOfRegionStatistics< SeedRgDirectValueFunctor< PixelType > >
      stats(detail::RequiresExplicitCast<unsigned int>::cast(max_region_label));
    seededRegionGrowing(srcImageRange(image), srcImage(res),
        destImage(res), stats, srgType);
    return res;
}



template < class PixelType, typename DestPixelType >
NumpyAnyArray pythonSeededRegionGrowingSeeded2D(NumpyArray<2, Singleband<PixelType> > image, 
                                          NumpyArray<2, Singleband<Int32> > seeds,
                                          SRGType srgType=CompleteGrow,
                                          NumpyArray<2, Singleband<DestPixelType> > res=python::object())
{
    vigra_precondition(image.shape() == seeds.shape(),
         "seededRegionGrowing(): magnitude and seed images must have the same "
         "size.");
    res.reshapeIfEmpty(image.shape(), "seededRegionGrowingSeeded2D(): Output array has wrong shape.");
    
    FindMinMax< DestPixelType > minmax;
    inspectImage(srcImageRange(seeds), minmax);
 
    ArrayOfRegionStatistics< SeedRgDirectValueFunctor< PixelType > >
        stats((unsigned int) std::ceil((double)minmax.max));
    seededRegionGrowing(srcImageRange(image), srcImage(seeds),
        destImage(res), stats, srgType);
    return res;
}



template < class PixelType >
NumpyAnyArray pythonLocalMinima2D(NumpyArray<2, Singleband<PixelType> > image,
    PixelType marker = 1,
    int neighborhood = 8,
    NumpyArray<2, Singleband<PixelType> > res = python::object()
    )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "localMinima2D(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 4:
        {
            localMinima(srcImageRange(image), destImage(res), marker,
                FourNeighborCode());
            break;
        }
        case 8:
        {
            localMinima(srcImageRange(image), destImage(res), marker,
                EightNeighborCode());
            break;
        }
    }

    return res;
}



template < class PixelType >
NumpyAnyArray pythonExtendedLocalMinima2D(NumpyArray<2, Singleband<PixelType> > image,
    PixelType marker = 1,
    int neighborhood = 8,
    NumpyArray<2, Singleband<PixelType> > res = python::object()
    )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "localMinima2D(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 4:
        {
            extendedLocalMinima(srcImageRange(image), destImage(res),
                marker, FourNeighborCode());
            break;
        }
        case 8:
        {
            extendedLocalMinima(srcImageRange(image), destImage(res),
                marker, EightNeighborCode());
            break;
        }
    }
    return res;
}



template < class PixelType >
NumpyAnyArray pythonLocalMaxima2D(NumpyArray<2, Singleband<PixelType> > image,
    PixelType marker = 1,
    int neighborhood = 8,
    NumpyArray<2, Singleband<PixelType> > res=python::object()
    )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "localMinima2D(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 4:
        {
            localMinima(srcImageRange(image), destImage(res), marker,
                FourNeighborCode());
            break;
        }
        case 8:
        {
            localMinima(srcImageRange(image), destImage(res), marker,
                EightNeighborCode());
            break;
        }
    }

    return res;
}



template < class PixelType >
NumpyAnyArray pythonExtendedLocalMaxima2D(NumpyArray<2, Singleband<PixelType> > image,
    PixelType marker = 1,
    int neighborhood = 8,
    NumpyArray<2, Singleband<PixelType> > res = python::object()
    )
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "ExtendedLocalMaxima2D(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 4:
        {
            extendedLocalMaxima(srcImageRange(image), destImage(res),
                marker, FourNeighborCode());
            break;
        }
        case 8:
        {
            extendedLocalMaxima(srcImageRange(image), destImage(res),
                marker, EightNeighborCode());
            break;
        }
    }
    return res;
}


template < class PixelType >
NumpyAnyArray pythonLabelImage2D(NumpyArray<2, Singleband<PixelType> > image,
    int neighborhood = 4,
    NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "labelImage2D(): Output array has wrong shape.");

    switch (neighborhood)
    {
        case 4:
        {
            labelImage(srcImageRange(image), destImage(res), false);
            break;
        }
        case 8:
        {
            labelImage(srcImageRange(image), destImage(res), true);
            break;
        }
    }

    return res;
}


template < class PixelType >
NumpyAnyArray pythonLabelImageWithBackground2D(NumpyArray<2, Singleband<PixelType> > image,
    int neighborhood = 4,
    PixelType background_value = 0,
    NumpyArray<2, Singleband<PixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "labelImageWithBackground2D(): Output array has wrong shape.");

    switch (neighborhood)
    {
        case 4:
        {
            labelImageWithBackground(srcImageRange(image),
                destImage(res), false, background_value);
            break;
        }
        case 8:
        {
            labelImageWithBackground(srcImageRange(image),
                destImage(res), true, background_value);
            break;
        }
    }
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

template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray pythonBoundaryTensor2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                               double scale,
                               NumpyArray<2, TinyVector<DestPixelType, 3> > res = python::object())
{
    res.reshapeIfEmpty(image.shape(), "boundaryTensor2D(): Output array has wrong shape.");    

    boundaryTensor(srcImageRange(image), destImage(res), scale);
     
    return res;
}
/** Wait for Kernel2D export
template < class SrcPixelType, typename DestPixelType >
NumpyAnyArray gradientEnergyTensor2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                     Kernel2D derivKernel, 
                                     Kernel2D smoothKernel,
                                     NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "gradientEnergyTensor2D(): Output array has wrong shape.");    
    
    gradientEnergyTensor(srcImageRange(image), destImage(res),
        derivKernel, smoothKernel);
     
    return res;
}
*/

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

/*template < class SrcPixelType, class DestPixelType>
NumpyAnyArray ellipticGaussian2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                 double sigmamax, double sigmamin,
                                 NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(image.shape(0), image.shape(1)), "ellipticGaussian2D(): Output array has wrong shape.");    
    
    ellipticGaussian(srcImageRange(image), destImage(res), sigmamax,
        sigmamin);
     
    return res;
}*/

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




void defineAnalysisFunctions()
{
    using namespace python;
    
    enum_<vigra::SRGType>("SRGType")
        .value("CompleteGrow", vigra::CompleteGrow)
        .value("CompleteGrowContours", vigra::KeepContours)
        ;

    multidef("watersheds", pywatersheds2D< UInt8, float >(),
      (arg("image"), arg("neighborhood") = 8, arg("out")=python::object()));

    def("seededRegionGrowingSeeded2D",
        registerConverters(&pythonSeededRegionGrowingSeeded2D<float, Int32>),
        (arg("image"), 
         arg("seeds"), 
         arg("srgType")=CompleteGrow,
         arg("out")=python::object()),
        "Region Segmentation by means of Seeded Region Growing.\n"
        "\n"
        "This algorithm implements seeded region growing as described in\n"
        "R. Adams, L. Bischof: \"<em> Seeded Region Growing</em>\", IEEE Trans. on Pattern Analysis and Maschine Intelligence, vol 16, no 6, 1994, and\n"
        "Ullrich Koethe: Primary Image Segmentation, in: G. Sagerer, S. Posch, F. Kummert (eds.): Mustererkennung 1995, Proc. 17. DAGM-Symposium, Springer 1995\n"
        "\n"
        "In differnece to seededReegionGrowing2D, the have to be provided in the seeds image\n"
        "\n"
        "For details see seededRegionGrowing_ in the vigra C++ documentation."
        );
    def("seededRegionGrowingSeeded2D",
        registerConverters(&pythonSeededRegionGrowingSeeded2D<float, Int64>),
        (arg("image"), 
         arg("seeds"), 
         arg("srgType")=CompleteGrow,
         arg("out")=python::object()));
    
    
    def("seededRegionGrowing2D",
        registerConverters(&pythonSeededRegionGrowing2D<float, Int32>),
        (arg("image"),
         arg("srgType")=CompleteGrow,
         arg("out")=python::object()),
        "Region Segmentation by means of Seeded Region Growing.\n"
        "\n"
        "This algorithm implements seeded region growing as described in\n"
        "R. Adams, L. Bischof: \"<em> Seeded Region Growing</em>\", IEEE Trans. on Pattern Analysis and Maschine Intelligence, vol 16, no 6, 1994, and\n"
        "Ullrich Koethe: Primary Image Segmentation, in: G. Sagerer, S. Posch, F. Kummert (eds.): Mustererkennung 1995, Proc. 17. DAGM-Symposium, Springer 1995\n"
        "\n"
        "In differnece to seededReegionGrowing2DSeeded, the seeds are calculated by finding local minimax\n"
        "\n"
        "For details see seededRegionGrowing_ in the vigra C++ documentation."
        );
    def("seededRegionGrowing2D",
        registerConverters(&pythonSeededRegionGrowing2D<float, Int64>),
        (arg("image"),
         arg("srgType")=CompleteGrow,
         arg("out")));

    def("localMinima2D",
        registerConverters(&pythonLocalMinima2D<float>),
        (arg("image"), 
         arg("marker")=1, 
         arg("neighborhood2D") = 8,
         arg("out")=python::object()),
        "Find local minima in an image.\n\n"
        "For details see localMinima_ in the vigra C++ documentation.");

    def("extendedLocalMinima2D",
        registerConverters(&pythonExtendedLocalMinima2D<float>),
        (arg("image"), 
         arg("marker")=1, 
         arg("neighborhood2D") = 8,
         arg("out")=python::object()),
        "Find local minima in an image.\n\n"
        "For details see extendedLocalMinima_ in the vigra C++ documentation."
        );

    def("localMaxima2D",
        registerConverters(&pythonLocalMaxima2D<float>),
        (arg("image"), 
         arg("marker")=1, 
         arg("neighborhood2D") = 8,
         arg("out")=python::object()),
        "Find local maxima in an image.\n\n"
        "For details see localMinima_ in the vigra C++ documentation.");

    def("extendedLocalMaxima2D",
        registerConverters(&pythonExtendedLocalMaxima2D<float>),
        (arg("image"), 
         arg("marker")=1, 
         arg("neighborhood2D") = 8,
         arg("out")=python::object()),
        "Find local maxima in an image.\n\n"
        "For details see localMinima_ in the vigra C++ documentation.");

    def("labelImage",
        registerConverters(&pythonLabelImage2D<float>),
        (arg("image"), 
        arg("neighborhood2D") = 4,
        arg("out")=python::object()),
        "Find the connected components of a segmented image."  );

    def("labelImageWithBackground",
        registerConverters(&pythonLabelImageWithBackground2D<float>),
        (arg("image"), 
        arg("neighborhood2D") = 4,
        arg("background_value") = 0,
        arg("out")=python::object()),
        "Find the connected components of a segmented image, excluding the background from labeling.\n\n"
        "For detailt see labelImageWithBackground int the vigra C++ documentation.");

    def("regionImageToCrackEdgeImage",
         registerConverters(&pythonRegionImageToCrackEdgeImage2D<Int32>),
         (arg("image"), 
          arg("edgeLabel") = 0, 
          arg("out")=python::object()),
         "Transform a labeled image into a crack edge image. \n\n"
         "For details see regionImageToCrackEdgeImage_ in the vigra C++ documentation.");

    def("regionImageToCrackEdgeImage",
         registerConverters(&pythonRegionImageToCrackEdgeImage2D<Int64>),
         (arg("image"), 
          arg("edgeLabel") = 0, 
          arg("out")=python::object()));

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

    def("distanceTransform2D",
        registerConverters(&pythonDistanceTransform2D<Int32,float>),
        (arg("image"), 
         arg("background")=0, 
         arg("norm")=2,
         arg("out")=python::object()),
        "For all background pixels, calculate the distance to the nearest object or contour.\n"
        "The label of the pixels to be considered background in the source image is passed in the parameter 'background'.\n"
        "Source pixels with other labels will be considered objects.\n"
        "In the destination image, all pixels corresponding to background will be assigned the their distance value, all pixels corresponding to objects will be assigned 0.\n\n"
        "The norm parameter gives the distance norm to use.\n\n"
        "For details see distanceTransform_ in the vigra C++ documentation.");
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
        "[G. Loy, A. Zelinsky: \"A Fast Radial Symmetry Transform for Detecting Points of Interest\", in: A. Heyden et al. (Eds.): Proc. of 7th European Conf. on Computer Vision, Part 1, pp. 358-368, Springer LNCS 2350, 2002]\n\n"
        "For details see radialSymmetryTransform_ in the vigra C++ documentation.");

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
    def("rieszTransformOfLOG2D",
        registerConverters(&pythonRieszTransformOfLOG2D<float>),        // also multiband
        (arg("image"), arg("scale"), arg("xorder"), arg("yorder"),arg("out")=python::object()),
        "Calculate Riesz transforms of the Laplacian of Gaussian.\n\n"
        "For details see rieszTransformOfLOG_ in the vigra C++ documentation.");

}

} // namespace vigra

