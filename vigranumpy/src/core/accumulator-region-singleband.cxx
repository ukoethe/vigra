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

namespace python = boost::python;

namespace vigra
{

template <unsigned int N, class T>
python::dict
extractConvexHullFeatures(NumpyArray<N, Singleband<T> > const & labels, python::object ignore_label)
{
    using namespace vigra::acc;

    TinyVector<npy_intp, N> permutation = labels.template permuteLikewise<N>();

    AccumulatorChainArray<CoupledArrays<N, T>, 
                          Select<ConvexHull, DataArg<1>, LabelArg<1> >
                         > acc;

    if(ignore_label != python::object())
        acc.ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
        
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
            array(k) = get<Count>(acc, k);
        }
        res["InputCount"] = array;
    }
    
    #define VIGRA_CONVEX_HULL_FEATURE(TYPE, NAME, FUNCTION) \
    { \
        NumpyArray<1, TYPE> array((Shape1(size))); \
        for(int k=0; k<size; ++k) \
        { \
            if(get<Count>(acc, k) == 0) \
                continue; \
            array(k) = get<ConvexHull>(acc, k).FUNCTION(); \
        } \
        res[#NAME] = array; \
    }
    
    VIGRA_CONVEX_HULL_FEATURE(double, InputPerimeter, inputPerimeter)
    VIGRA_CONVEX_HULL_FEATURE(double, InputArea, inputArea)
    VIGRA_CONVEX_HULL_FEATURE(double, Perimeter, hullPerimeter)
    VIGRA_CONVEX_HULL_FEATURE(double, Area, hullArea)
    VIGRA_CONVEX_HULL_FEATURE(double, Convexity, convexity)
    VIGRA_CONVEX_HULL_FEATURE(double, Rugosity, rugosity)
    VIGRA_CONVEX_HULL_FEATURE(npy_uint32, DefectCount, convexityDefectCount)
    VIGRA_CONVEX_HULL_FEATURE(double, DefectAreaMean, convexityDefectAreaMean)
    VIGRA_CONVEX_HULL_FEATURE(double, DefectAreaVariance, convexityDefectAreaVariance)
    VIGRA_CONVEX_HULL_FEATURE(double, DefectAreaSkewness, convexityDefectAreaSkewness)
    VIGRA_CONVEX_HULL_FEATURE(double, DefectAreaKurtosis, convexityDefectAreaKurtosis)
    VIGRA_CONVEX_HULL_FEATURE(double, MeanDefectDisplacement, meanDefectDisplacement)
    
    #undef VIGRA_CONVEX_HULL_FEATURE
    
    {
        python::list hulls;
        for(int k=0; k<size; ++k)
        {
            if(get<Count>(acc, k) == 0)
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
        res["Polygon"] = hulls;
    }
    
    {
        NumpyArray<2, double> array(Shape2(size, 3));
        for(int k=0; k<size; ++k)
        {
            if(get<Count>(acc, k) == 0)
                continue;
            int defects = min<int>(3, get<ConvexHull>(acc, k).convexityDefectCount());
            for(int j=0; j<defects; ++j)
                array(k, j) = get<ConvexHull>(acc, k).defectAreaList()[j];
        }
        res["DefectAreaList"] = array;
    }
    
    #define VIGRA_CONVEX_HULL_VECTOR_FEATURE(NAME, FUNCTION) \
    { \
        NumpyArray<2, double> array(Shape2(size, N)); \
        for(int k=0; k<size; ++k) \
        { \
            if(get<Count>(acc, k) == 0) \
                continue; \
            for(int j=0; j<N; ++j) \
                array(k, permutation[j]) = get<ConvexHull>(acc, k).FUNCTION()[j]; \
        } \
        res[#NAME] = array; \
    }
    
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(InputCenter, inputCenter)
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(Center, hullCenter)
    VIGRA_CONVEX_HULL_VECTOR_FEATURE(DefectCenter, convexityDefectCenter)
    
    #undef VIGRA_CONVEX_HULL_VECTOR_FEATURE

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
          (arg("labels"), arg("ignoreLabel")=python::object()),
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
}

} // namespace vigra
