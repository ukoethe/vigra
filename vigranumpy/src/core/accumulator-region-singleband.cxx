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

// Helper functions for stringification of macro string constants
#define STR(s) #s
#define XSTR(s) s

template <unsigned int N, class T>
python::object
extractConvexHullFeatures(NumpyArray<N, Singleband<T> > const & labels,
			  python::object ignore_label,
                          bool list_features_only=false)
{
    using namespace vigra::acc;

    #define VIGRA_CONVEX_HULL_FEATURE_INPUT_COUNT "Input Count"
    #define VIGRA_CONVEX_HULL_FEATURE_INPUT_PERIMETER "Input Perimeter"
    #define VIGRA_CONVEX_HULL_FEATURE_INPUT_AREA "Input Area"
    #define VIGRA_CONVEX_HULL_FEATURE_AREA "Area"
    #define VIGRA_CONVEX_HULL_FEATURE_PERIMETER "Perimeter"
    #define VIGRA_CONVEX_HULL_FEATURE_RUGOSITY "Rugosity"
    #define VIGRA_CONVEX_HULL_FEATURE_CONVEXITY "Convexity"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_COUNT "Defect Count"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_MEAN_DISPLACEMENT "Defect Mean Displacement"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_LIST "Defect Area List"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_MEAN "Defect Area Mean"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_VARIANCE "Defect Area Variance"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_SKEWNESS "Defect Area Skewness"
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_KURTOSIS "Defect Area Kurtosis"
    #define VIGRA_CONVEX_HULL_FEATURE_POLYGON "Polygon"

    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_INPUT_CENTER "Input Center"
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_CENTER "Center"
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT_CENTER "Defect Center"

    if(list_features_only)
    {

        python::list res;
	res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_INPUT_COUNT));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_INPUT_PERIMETER));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_INPUT_AREA));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_AREA));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_PERIMETER));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_RUGOSITY));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_CONVEXITY));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_POLYGON));

        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_COUNT));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_MEAN_DISPLACEMENT));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_LIST));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_MEAN));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_VARIANCE));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_SKEWNESS));
        res.append(XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_KURTOSIS));

        res.append(XSTR(VIGRA_CONVEX_HULL_VECTOR_FEATURE_INPUT_CENTER));
        res.append(XSTR(VIGRA_CONVEX_HULL_VECTOR_FEATURE_CENTER));

        res.append(XSTR(VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT_CENTER));

        return res;
    }

    TinyVector<npy_intp, N> permutation = labels.template permuteLikewise<N>();

    AccumulatorChainArray<CoupledArrays<N, T>, 
                          Select<ConvexHull, DataArg<1>, LabelArg<1> >
                         > acc;

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
    python::dict res;
    {
        NumpyArray<1, npy_uint32> array((Shape1(size)));
        for(int k=0; k<size; ++k)
        {
            if(k == ignored_label)
                continue;
            array(k) = get<Count>(acc, k);
        }
	res[XSTR(VIGRA_CONVEX_HULL_FEATURE_INPUT_COUNT)] = array;
    }
    
    #define VIGRA_CONVEX_HULL_FEATURE(TYPE, NAME, FUNCTION) \
    { \
        NumpyArray<1, TYPE> array((Shape1(size)));	\
        for(int k=0; k<size; ++k) \
        { \
            if(k == ignored_label || get<Count>(acc, k) == 0) \
                continue; \
            array(k) = get<ConvexHull>(acc, k).FUNCTION(); \
        } \
	res[XSTR(NAME)] = array;			\
    }
    
    #define VIGRA_CONVEX_HULL_FEATURE_DEFECT(TYPE, NAME, FUNCTION) \
    { \
        NumpyArray<1, double> array((Shape1(size)));	\
        for(int k=0; k<size; ++k) \
        { \
	    if(k == ignored_label || get<Count>(acc, k) == 0) \
	      continue; \
	    array(k) = get<ConvexHull>(acc, k).meanDefectDisplacement(); \
        } \
	res[XSTR(NAME)] = array;			\
    }

    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_INPUT_PERIMETER, inputPerimeter)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_INPUT_AREA, inputArea)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_PERIMETER, hullPerimeter)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_AREA, hullArea)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_CONVEXITY, convexity)
    VIGRA_CONVEX_HULL_FEATURE(double, VIGRA_CONVEX_HULL_FEATURE_RUGOSITY, rugosity)
    VIGRA_CONVEX_HULL_FEATURE(npy_uint32, VIGRA_CONVEX_HULL_FEATURE_DEFECT_COUNT, convexityDefectCount)

    VIGRA_CONVEX_HULL_FEATURE_DEFECT(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_MEAN, convexityDefectAreaMean)
    VIGRA_CONVEX_HULL_FEATURE_DEFECT(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_MEAN_DISPLACEMENT, meanDefectDisplacement)
    VIGRA_CONVEX_HULL_FEATURE_DEFECT(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_VARIANCE, convexityDefectAreaVariance)
    VIGRA_CONVEX_HULL_FEATURE_DEFECT(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_SKEWNESS, convexityDefectAreaSkewness)
    VIGRA_CONVEX_HULL_FEATURE_DEFECT(double, VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_KURTOSIS, convexityDefectAreaKurtosis)
    
    #undef VIGRA_CONVEX_HULL_FEATURE
    
    {
        python::list hulls;
        for(int k=0; k<size; ++k)
        {
            if(k == ignored_label || get<Count>(acc, k) == 0)
            {
                hulls.append(python::object());
                continue;
            }
            int hull_size = get<ConvexHull>(acc, k).hull().size();
            NumpyArray<2, double> array(Shape2(hull_size, N));
            Polygon<TinyVector<double, 2> > poly = (permutation == Shape2(0,1))
                                                       ? get<ConvexHull>(acc, k).hull()
                                                       : reverse(transpose(get<ConvexHull>(acc, k).hull()));
            for(int p=0; p<hull_size; ++p)
            {
                for(int j=0; j<N; ++j)
                    array(p, j) = poly[p][j];
            }
            hulls.append(array);
        }
	res[XSTR(VIGRA_CONVEX_HULL_FEATURE_POLYGON)] = hulls;
    }
    
    {
        NumpyArray<2, double> array(Shape2(size, 3));
        for(int k=0; k<size; ++k)
        {
            if(k == ignored_label || get<Count>(acc, k) == 0)
                continue;
            int defects = min<int>(3, get<ConvexHull>(acc, k).convexityDefectCount());
            for(int j=0; j<defects; ++j)
                array(k, j) = get<ConvexHull>(acc, k).defectAreaList()[j];
        }
        res[XSTR(VIGRA_CONVEX_HULL_FEATURE_DEFECT_AREA_LIST)] = array;
    }
    
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE(NAME, FUNCTION) \
    { \
        NumpyArray<2, double> array(Shape2(size, N)); \
        for(int k=0; k<size; ++k) \
        { \
	    if(k == ignored_label || get<Count>(acc, k) == 0)	\
	        continue; \
            for(int j=0; j<N; ++j) \
		array(k, permutation[j]) = get<ConvexHull>(acc, k).FUNCTION()[j]; \
        } \
	res[XSTR(NAME)] = array; \
    }
    
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT(NAME, FUNCTION) \
    { \
        NumpyArray<2, double> array(Shape2(size, N)); \
        for(int k=0; k<size; ++k) \
        { \
	    if(k == ignored_label || get<Count>(acc, k) == 0)	\
	        continue; \
            for(int j=0; j<N; ++j) \
		array(k, permutation[j]) = get<ConvexHull>(acc, k).FUNCTION()[j]; \
        } \
	res[XSTR(NAME)] = array; \
    }
    
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(VIGRA_CONVEX_HULL_VECTOR_FEATURE_INPUT_CENTER, inputCenter)
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(VIGRA_CONVEX_HULL_VECTOR_FEATURE_CENTER, hullCenter)

    VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT(VIGRA_CONVEX_HULL_VECTOR_FEATURE_DEFECT_CENTER, convexityDefectCenter)

    #undef VIGRA_CONVEX_HULL_VECTOR_FEATURE

    return res;
}

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
    #define VIGRA_SKELETON_VECTOR_FEATURE_CENTER "Center"
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
            for(int j=0; j<N; ++j) \
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
    
    def("extractConvexHullFeatures", 
         registerConverters(&extractConvexHullFeatures<2, npy_uint32>),
          (arg("labels"),
	   arg("ignoreLabel")=python::object(),
           arg("list_features_only")=false),
            "\nExtract convex hull features for each region of a labeled 2D image\n"
            "(with dtype=numpy.uint32) and return a dictionary holding the\n"
            "resulting feature arrays. Argument 'ignoreLabel' can be used to specify\n"
            "an optional background label that is to be skipped. Note that the\n"
            "convex hull itself and its features are computed from the interpixel\n"
            "contour around each region. In the following, 'convexity defects'\n"
            "are defined as the connected components of the set difference\n"
            "between the convex hull and the original region.\n\n"
            "The result dictionary holds the following keys:\n\n"
            "   - 'InputCount':  the number of pixels in the original region\n\n"
            "   - 'InputPerimeter':  the perimeter of the original interpixel contour\n\n"
            "   - 'InputArea':  the areay enclosed by the original interpixel contour\n\n"
            "   - 'InputCenter':  the centroid of the original interpixel contour polygon\n\n"
            "   - 'Perimeter':  the perimeter of the convex hull polygon\n\n"
            "   - 'Area':  the area enclosed by the convex hull polygon\n\n"
            "   - 'Center':  the centroid of the convex hull polygon\n\n"
            "   - 'Rugosity':  ratio between original perimeter and hull perimeter (>= 1)\n\n"
            "   - 'Convexity':  the ratio between hull area and original area (<= 1)\n\n"
            "   - 'DefectCount':  the number of convexity defects\n\n"
            "   - 'DefectCenter':  the combined centroid of the defects\n\n"
            "   - 'MeanDefectDisplacement':  mean distance between the centroids of the\n"
            "                                original object and the centroids of the defects,\n"
            "                                weighted by defect area\n\n"
            "   - 'DefectAreaList':  the area of the three largest convexity defects\n\n"
            "   - 'DefectAreaMean':  mean of the convexity defect areas\n\n"
            "   - 'DefectAreaVariance':  variance of the convexity defect areas\n\n"
            "   - 'DefectAreaSkewness':  skewness of the convexity defect areas\n\n"
            "   - 'DefectAreaKurtosis':  kurtosis of the convexity defect areas\n\n"
            "   - 'Polygon':  the convex hull polygon\n\n");
    
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
