/************************************************************************/
/*                                                                      */
/*            Copyright 2011-2012 by Ullrich Koethe                     */
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

#include "pythonaccumulator.hxx"
#include <vigra/skeleton.hxx>

namespace python = boost::python;

namespace vigra
{

// workaround for compiler bug in VS 2015 (compiler fails to match the template
// function get_pointer() at line 20 of boost/get_pointer.hpp)
namespace acc
{
	inline PythonRegionFeatureAccumulator const volatile *
    get_pointer(PythonRegionFeatureAccumulator const volatile * p) { return p; }
}

// Helper functions for stringification of macro string constants
#define STR(s) #s
#define XSTR(s) s

#ifdef WITH_LEMON

template <unsigned int N, class T>
python::object
extractConvexHullFeatures(NumpyArray<N, Singleband<T> > const & labels,
                          python::object ignore_label,
                          bool list_features_only=false)
{
    using namespace vigra::acc;

    #define VIGRA_CONVEX_HULL_FEATURE_INPUT_VOLUME "InputVolume"
    #define VIGRA_CONVEX_HULL_FEATURE_HULL_VOLUME "HullVolume"
    #define VIGRA_CONVEX_HULL_FEATURE_CONVEXITY "Convexity"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_MEAN "DefectVolumeMean"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_VARIANCE "DefectVolumeVariance"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_SKEWNESS "DefectVolumeSkewness"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_KURTOSIS "DefectVolumeKurtosis"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_COUNT "DefectCount"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_DISPLACEMENT_MEAN "DefectDisplacementMean"

    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_INPUT_CENTER "InputCenter"
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_HULL_CENTER "HullCenter"
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT_CENTER "DefectCenter"

    if(list_features_only)
    {
        python::list res;
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_INPUT_VOLUME));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_HULL_VOLUME));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_CONVEXITY));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_MEAN));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_VARIANCE));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_SKEWNESS));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_KURTOSIS));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_COUNT));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_DISPLACEMENT_MEAN));

        res.append(XSTR(VIGRA_CONVEX_HULL_VECTOR_FEATURE_INPUT_CENTER));
        res.append(XSTR(VIGRA_CONVEX_HULL_VECTOR_FEATURE_HULL_CENTER));
        res.append(XSTR(VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT_CENTER));

        return res;
    }

    TinyVector<npy_intp, N> permutation = labels.template permuteLikewise<N>();

    AccumulatorChainArray<
            CoupledArrays<N, T>,
            Select<ConvexHullFeatures, DataArg<1>, LabelArg<1> > > acc;

    MultiArrayIndex ignored_label = -1;
    if(ignore_label != python::object())
    {
        ignored_label = python::extract<MultiArrayIndex>(ignore_label)();
        acc.ignoreLabel(ignored_label);
    }

    {
        PyAllowThreads _pythread;
        extractFeatures(labels, acc);
    }

    int size = acc.maxRegionLabel()+1;

    // finalize the calculations
    for (int k = 0; k < size; ++k)
    {
        if (k != ignored_label && get<Count>(acc, k) != 0)
        {
            getAccumulator<ConvexHullFeatures>(acc, k).finalize();
        }
    }

    // initialize return dict
    python::dict res;

    #define VIGRA_CONVEX_HULL_FEATURE(TYPE, NAME, FUNCTION) \
    { \
        NumpyArray<1, TYPE> array((Shape1(size)));        \
        for(int k=0; k<size; ++k) \
        { \
            if(k == ignored_label || get<Count>(acc, k) == 0) \
                continue; \
            array(k) = get<ConvexHullFeatures>(acc, k).FUNCTION(); \
        } \
        res[XSTR(NAME)] = array;                        \
    }

    VIGRA_CONVEX_HULL_FEATURE(npy_uint32, VIGRA_CONVEX_HULL_FEATURE_INPUT_VOLUME, inputVolume)
    VIGRA_CONVEX_HULL_FEATURE(npy_uint32, VIGRA_CONVEX_HULL_FEATURE_HULL_VOLUME, hullVolume)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_CONVEXITY, convexity)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_MEAN, defectVolumeMean)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_VARIANCE, defectVolumeVariance)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_SKEWNESS, defectVolumeSkewness)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_VOLUME_KURTOSIS, defectVolumeKurtosis)
    VIGRA_CONVEX_HULL_FEATURE(npy_uint32, VIGRA_CONVEX_HULL_FEATURE_DEFECT_COUNT, defectCount)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_DISPLACEMENT_MEAN, defectDisplacementMean)

    #undef VIGRA_CONVEX_HULL_FEATURE

    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE(NAME, FUNCTION) \
    { \
        NumpyArray<2, double> array(Shape2(size, N)); \
        for(int k=0; k<size; ++k) \
        { \
            if(k == ignored_label || get<Count>(acc, k) == 0)        \
                continue; \
            for(unsigned j=0; j<N; ++j) \
                array(k, permutation[j]) = get<ConvexHullFeatures>(acc, k).FUNCTION()[j]; \
        } \
        res[XSTR(NAME)] = array; \
    }

    VIGRA_CONVEX_HULL_VECTOR_FEATURE(VIGRA_CONVEX_HULL_VECTOR_FEATURE_INPUT_CENTER, inputCenter)
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(VIGRA_CONVEX_HULL_VECTOR_FEATURE_HULL_CENTER, hullCenter)
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT_CENTER, defectCenter)

    #undef VIGRA_CONVEX_HULL_VECTOR_FEATURE

    return res;
}

#endif // WITH_LEMON

template <unsigned int N, class T>
python::object
pyExtractSkeletonFeatures(NumpyArray<N, Singleband<T> > const & labels,
                          double pruning_threshold,
                          bool list_features_only=false)
{
    using namespace vigra::acc;

    #define VIGRA_SKELETON_FEATURE_DIAMETER "Diameter"
    #define VIGRA_SKELETON_FEATURE_EUCLIDEAN_DIAMETER "Euclidean Diameter"
    #define VIGRA_SKELETON_FEATURE_TOTAL_LENGTH "Total Length"
    #define VIGRA_SKELETON_FEATURE_AVERAGE_LENGTH "Average Length"
    #define VIGRA_SKELETON_FEATURE_BRANCH_COUNT "Branch Count"
    #define VIGRA_SKELETON_FEATURE_HOLE_COUNT "Hole Count"
    #define VIGRA_SKELETON_VECTOR_FEATURE_CENTER "Skeleton Center"
    #define VIGRA_SKELETON_VECTOR_FEATURE_TERMINAL_1 "Terminal 1"
    #define VIGRA_SKELETON_VECTOR_FEATURE_TERMINAL_2 "Terminal 2"

    if(list_features_only)
    {

        python::list res;
        res.append(XSTR(VIGRA_SKELETON_FEATURE_DIAMETER));
        res.append(XSTR(VIGRA_SKELETON_FEATURE_EUCLIDEAN_DIAMETER));
        res.append(XSTR(VIGRA_SKELETON_FEATURE_TOTAL_LENGTH));
        res.append(XSTR(VIGRA_SKELETON_FEATURE_AVERAGE_LENGTH));
        res.append(XSTR(VIGRA_SKELETON_FEATURE_BRANCH_COUNT));
        res.append(XSTR(VIGRA_SKELETON_FEATURE_HOLE_COUNT));
        res.append(XSTR(VIGRA_SKELETON_VECTOR_FEATURE_CENTER));
        res.append(XSTR(VIGRA_SKELETON_VECTOR_FEATURE_TERMINAL_1));
        res.append(XSTR(VIGRA_SKELETON_VECTOR_FEATURE_TERMINAL_2));

        return res;
    }

    TinyVector<npy_intp, N> permutation = labels.template permuteLikewise<N>();
    ArrayVector<SkeletonFeatures> features;

    {
        PyAllowThreads _pythread;

        extractSkeletonFeatures(labels, features,
                                SkeletonOptions().pruneSalienceRelative(pruning_threshold));
    }

    int size = features.size();
    python::dict res;

    #define VIGRA_SKELETON_FEATURE(TYPE, NAME, ATTRIBUTE) \
    { \
        NumpyArray<1, TYPE> array((Shape1(size))); \
        for(int k=0; k<size; ++k) \
        { \
            array(k) = features[k].ATTRIBUTE; \
        } \
        res[XSTR(NAME)] = array; \
    }

    VIGRA_SKELETON_FEATURE(double, VIGRA_SKELETON_FEATURE_DIAMETER, diameter)
    VIGRA_SKELETON_FEATURE(double, VIGRA_SKELETON_FEATURE_EUCLIDEAN_DIAMETER, euclidean_diameter)
    VIGRA_SKELETON_FEATURE(double, VIGRA_SKELETON_FEATURE_TOTAL_LENGTH, total_length)
    VIGRA_SKELETON_FEATURE(double, VIGRA_SKELETON_FEATURE_AVERAGE_LENGTH, average_length)
    VIGRA_SKELETON_FEATURE(npy_uint32, VIGRA_SKELETON_FEATURE_BRANCH_COUNT, branch_count)
    VIGRA_SKELETON_FEATURE(npy_uint32,VIGRA_SKELETON_FEATURE_HOLE_COUNT, hole_count)

    #undef VIGRA_SKELETON_FEATURE

    #define VIGRA_SKELETON_VECTOR_FEATURE(NAME, ATTRIBUTE) \
    { \
        NumpyArray<2, double> array(Shape2(size, N)); \
        for(int k=0; k<size; ++k) \
        { \
            for(unsigned j=0; j<N; ++j) \
                array(k, permutation[j]) = features[k].ATTRIBUTE[j]; \
        } \
        res[XSTR(NAME)] = array; \
    }

    VIGRA_SKELETON_VECTOR_FEATURE(VIGRA_SKELETON_VECTOR_FEATURE_CENTER, center)
    VIGRA_SKELETON_VECTOR_FEATURE(VIGRA_SKELETON_VECTOR_FEATURE_TERMINAL_1, terminal1)
    VIGRA_SKELETON_VECTOR_FEATURE(VIGRA_SKELETON_VECTOR_FEATURE_TERMINAL_2, terminal2)

    #undef VIGRA_SKELETON_VECTOR_FEATURE

    return res;
}

void defineSinglebandRegionAccumulators()
{
    using namespace python;
    using namespace vigra::acc;

    docstring_options doc_options(true, true, false);
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis,
                   Minimum, Maximum, StandardQuantiles<GlobalRangeHistogram<0> >,
                   RegionCenter, RegionRadii, RegionAxes,
                   Weighted<RegionCenter>, Weighted<RegionRadii>, Weighted<RegionAxes>,
                   Select<Coord<Minimum>, Coord<Maximum>, Coord<ArgMinWeight>, Coord<ArgMaxWeight>,
                          Principal<Coord<Skewness> >, Principal<Coord<Kurtosis> >,
                          Principal<Weighted<Coord<Skewness> > >, Principal<Weighted<Coord<Kurtosis> > > >,
                   DataArg<1>, WeightArg<1>, LabelArg<2>
                   > ScalarRegionAccumulators;
    definePythonAccumulatorArraySingleband<2, float, ScalarRegionAccumulators>();
    definePythonAccumulatorArraySingleband<3, float, ScalarRegionAccumulators>();

#ifdef WITH_LEMON
    def("extract2DConvexHullFeatures",
        registerConverters(&extractConvexHullFeatures<2, npy_uint32>),
        (   arg("labels"),
            arg("ignoreLabel")=python::object(),
            arg("list_features_only")=false),
        "\nExtract convex hull features for each region of a labeled 2D image (with\n"
        "dtype=numpy.uint32) and return a dictionary holding the resulting feature\n"
        "arrays. The argument 'ignoreLabel' can be used to specify an optional\n"
        "background label that is to be skipped. Note that the convex hull itself and\n"
        "its features are computed from the interpixel contour around each region. In\n"
        "the following, 'convexity defects' are the connected components of the set\n"
        "difference between the convex hull and the original region.\n"
        "The result dictionary holds the following keys:\n\n"
        "   - InputVolume :  the number of pixels in the original region\n\n"
        "   - HullVolume : the number of pixels in the convex hull\n\n"
        "   - Convexity : the ratio between the convex hull volume and the input\n"
        "     volume\n\n"
        "   - DefectVolumeMean : mean of the volumes of the convexity defects\n\n"
        "   - DefectVolumeVariance : variance of the volumes of the convexity\n"
        "     defects\n\n"
        "   - DefectVolumeSkewness : skewness of the volumes of the convexity\n"
        "     defects\n\n"
        "   - DefectVolumeKurtosis : kurtosis of the volumes of the convexity\n"
        "     defects\n\n"
        "   - DefectCount : number of convexity defects\n\n"
        "   - DefectDisplacementMean : mean distance between the center of the input\n"
        "     region and the center of the defects, weighted by the defect volumes\n\n"
        "   - InputCenter : center of the input region\n\n"
        "   - HullCenter : center of the convex hull\n\n"
        "   - DefectCenter : center of the defects\n\n");

    def("extract3DConvexHullFeatures",
        registerConverters(&extractConvexHullFeatures<3, npy_uint32>),
        (   arg("labels"),
            arg("ignoreLabel")=python::object(),
            arg("list_features_only")=false),
        "\nExtract convex hull features for each region of a labeled 3D image (with\n"
        "dtype=numpy.uint32) and return a dictionary holding the resulting feature\n"
        "arrays. The argument 'ignoreLabel' can be used to specify an optional\n"
        "background label that is to be skipped. Note that the convex hull itself and\n"
        "its features are computed from the interpixel contour around each region. In\n"
        "the following, 'convexity defects' are the connected components of the set\n"
        "difference between the convex hull and the original region.\n"
        "The result dictionary holds the following keys:\n\n"
        "   - InputVolume :  the number of pixels in the original region\n\n"
        "   - HullVolume : the number of pixels in the convex hull\n\n"
        "   - Convexity : the ratio between the convex hull volume and the input\n"
        "     volume\n\n"
        "   - DefectVolumeMean : mean of the volumes of the convexity defects\n\n"
        "   - DefectVolumeVariance : variance of the volumes of the convexity\n"
        "     defects\n\n"
        "   - DefectVolumeSkewness : skewness of the volumes of the convexity\n"
        "     defects\n\n"
        "   - DefectVolumeKurtosis : kurtosis of the volumes of the convexity\n"
        "     defects\n\n"
        "   - DefectCount : number of convexity defects\n\n"
        "   - DefectDisplacementMean : mean distance between the center of the input\n"
        "     region and the center of the defects, weighted by the defect volumes\n\n"
        "   - InputCenter : center of the input region\n\n"
        "   - HullCenter : center of the convex hull\n\n"
        "   - DefectCenter : center of the defects\n\n");
#endif // WITH_LEMON

    def("extractSkeletonFeatures",
         registerConverters(&pyExtractSkeletonFeatures<2, npy_uint32>),
         (arg("labels"),
          arg("pruning_threshold")=0.2,
          arg("list_features_only")=false),
            "\nExtract skeleton features for each region of a labeled 2D image\n"
            "(with dtype=numpy.uint32) and return a dictionary holding the\n"
            "resulting feature arrays. Label 0 is always considered background\n"
            "and therefore skipped. The skeleton is computed using mode\n"
            "'PruneSalienceRelative' with the given 'pruning_threshold'.\n\n"
            "The result dictionary holds the following keys:\n\n"
            "   - 'Diameter':  the longest path between two terminals of the skeleton\n\n"
            "   - 'Center':  the center point of this path\n\n"
            "   - 'Terminal1':  first end point of this path\n\n"
            "   - 'Terminal2':  second end point of this path\n\n"
            "   - 'EuclideanDiameter':  the Euclidean distance between Terminal1 and Terminal2\n\n"
            "   - 'TotalLength':  total length of the (pruned) skeleton\n\n"
            "   - 'AverageLength':  the average length of the skeleton's branches after pruning\n\n"
            "   - 'BranchCount':  the number of skeleton branches (i.e. end points after pruning)\n\n"
            "   - 'HoleCount':  the number of cycles in the skeleton\n"
            "                  (i.e. the number of cavities in the region)\n\n");
}

} // namespace vigra
