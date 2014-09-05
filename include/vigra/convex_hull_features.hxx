/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2014 by                                 */
/*               Ullrich Koethe,                                        */
/*               Esteban Pardo                                          */
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

#ifndef VIGRA_CONVEX_HULL_FEATURES_HXX
#define VIGRA_CONVEX_HULL_FEATURES_HXX

#include "multi_array.hxx"
#include "accumulator.hxx"
#include "polygon.hxx"
#include "labelimage.hxx"
// #include "pointordering.hxx"

namespace vigra {

namespace detail {

template<typename T>
int countNonZero(MultiArray<2, T> const &array) {

    using namespace acc;
    MultiArray<2, double> data(array.shape());
    AccumulatorChainArray<CoupledArrays<2, double, T>,
            Select<LabelArg<2>, Count> > a;
    extractFeatures(data, array, a);

    int nonZero = array.size() - get<Count>(a, 0);
    return nonZero;
}

/*
 * Returns first pixel found on a border
 */
template<class IMAGE_TYPE>
bool findAnchorPoint(IMAGE_TYPE const &input_image,
        TinyVector<int, 2> &anchor_point) {

    TinyVector<int, 2> neighbor_coords[4] = { TinyVector<int, 2>(0, -1),
            TinyVector<int, 2>(1, 0), TinyVector<int, 2>(0, 1), TinyVector<int,
                    2>(-1, 0) };

    for (unsigned y = 1; y < input_image.height() - 1; ++y) {
        for (unsigned x = 1; x < input_image.width() - 1; ++x) {
            if (input_image(y, x)) {
                for (unsigned i = 0; i < 4; ++i) {
                    unsigned neighbor_y = y + neighbor_coords[i][1];
                    unsigned neighbor_x = x + neighbor_coords[i][0];
                    // It must have at least a background pixel
                    if (!input_image(neighbor_x, neighbor_y)) {
                        anchor_point = TinyVector<int, 2>(x, y);
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

} // namespace detail

/*
 * Calculates some features from the relation between a polygon, and its convex hull.
 */
template<class IMAGE_TYPE>
class ConvexHullFeatures {
public:

    ConvexHullFeatures(IMAGE_TYPE const &input_image) :
            featuresCalculated(false) {
        calculateFeatures(input_image);

    }

    void calculateFeatures(IMAGE_TYPE const &input_image) {

        using namespace vigra::acc;

        TinyVector<int, 2> anchor_point;
        detail::findAnchorPoint(input_image, anchor_point);

        std::vector<TinyVector<float, 2> > contour_points;

        extractContour(input_image, anchor_point, contour_points);

        // FIXME: this is not the true perimeter of the polygon
        inputPerimeter = contour_points.size();
        if(contour_points.front() == contour_points.back())
            --inputPerimeter;

        vigra::Polygon<TinyVector<float, 2> > convex_hull_points;
        convexHull(contour_points, convex_hull_points);
        std::cerr << "convex hull length: " << convex_hull_points.length() << " " <<
                      (10.0 + 2.0*M_SQRT2) << "\n";

        // for(auto i = convex_hull_points.begin(); i != convex_hull_points.end(); ++i)
        // {
            // std::cerr << *i << "\n";
        // }

        MultiArray<2, int> convex_hull_image(input_image.shape());
        
        fillPolygon(convex_hull_points, convex_hull_image, 1);
        
        // for(auto i = convex_hull_image.begin(); i != convex_hull_image.end(); ++i)
        // {
            // std::cerr << (*i ? "*" : " ");
            // if(i.point()[0] == convex_hull_image.shape(0)-1)
                // std::cerr << "\n";
        // }

        std::vector<TinyVector<float, 2> > convex_hull_contour_points;

        detail::findAnchorPoint(convex_hull_image, anchor_point);
        extractContour(convex_hull_image, anchor_point,
                       convex_hull_contour_points);

        // FIXME: this is not the true perimeter of the polygon
        convexHullPerimeter = convex_hull_contour_points.size() - 1;

        rugosity = (double) inputPerimeter / (double) convexHullPerimeter;

        MultiArray<2, double> diff_image(
                vigra::Shape2(input_image.height(), input_image.width()));

        combineTwoImages(input_image, convex_hull_image, diff_image,
                std::not_equal_to<int>());

        MultiArray<2, int> labels(
                vigra::Shape2(input_image.height(), input_image.width()));

        convexityDefectCount = labelImageWithBackground(diff_image, labels,
                false, 0);

        convexHullArea = detail::countNonZero(convex_hull_image);


        inputArea = detail::countNonZero(input_image);
        convexity = (double) inputArea / (double) convexHullArea;

        AccumulatorChainArray<CoupledArrays<2, double, int>,
                Select<DataArg<1>, LabelArg<2>, Count> > labelAccumulator;

        labelAccumulator.ignoreLabel(0);

        extractFeatures(diff_image, labels, labelAccumulator);

        MultiArray<1, int> convexityDefectAreas(convexityDefectCount);

        for (int i = 0; i < convexityDefectCount; ++i) {
            convexityDefectAreas(i) = get<Count>(labelAccumulator, i + 1); // Interesting labels start at 1 (0 means background)
        }

        AccumulatorChain<int, Select<Mean, Variance, Skewness, Kurtosis> > convexityAccumulator;
        extractFeatures(convexityDefectAreas.begin(),
                convexityDefectAreas.end(), convexityAccumulator);

        convexityDefectAreaMean = get<Mean>(convexityAccumulator);
        convexityDefectAreaVariance = get<Variance>(convexityAccumulator);
        convexityDefectAreaSkewness = get<Skewness>(convexityAccumulator);
        convexityDefectAreaKurtosis = get<Kurtosis>(convexityAccumulator);

        featuresCalculated = true;
    }

    /*
     * Returns the count of pixels contained in the convex hull
     */
    int getConvexHullArea() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexHullArea;
    }

    /*
     * Returns the count of pixels contained in the input polygon
     */
    int getInputArea() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return inputArea;
    }

    /*
     * Returns the ratio between the input area and the convex hull area
     * The closer to 1 the more convex the input polygon is
     */
    double getConvexity() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexity;
    }

    /*
     * Returns the number of convexity defects
     */
    int getConvexityDefectCount() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexityDefectCount;
    }

    /*
     * Returns the mean area of the convexity defects
     */
    double getConvexityDefectAreaMean() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexityDefectAreaMean;
    }

    /*
     * Returns the variance of the convexity defect areas
     */
    double getConvexityDefectAreaVariance() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexityDefectAreaVariance;
    }

    /*
     * Returns the skewness of the convexity defect areas
     */
    double getConvexityDefectAreaSkewness() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexityDefectAreaSkewness;
    }

    /*
     * Returns the kurtosis of the convexity defect areas
     */
    double getConvexityDefectAreaKurtosis() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexityDefectAreaKurtosis;
    }

    /*
     * Returns the perimeter of the input polygon
     * The perimeter is calculated as the number of pixel faces
     * surrounding the polygon.
     */
    int getInputPerimeter() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return inputPerimeter;
    }

    /*
     * Returns the perimeter of the convex hull
     * The perimeter is calculated as the number of pixel faces
     * surrounding the polygon.
     */
    int getConvexHullPerimeter() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return convexHullPerimeter;
    }

    /*
     * Returns the ration between the input perimeter and the convex perimeter
     * The higher the value is, the less convex the input polygon is
     */
    double getRugosity() {
        vigra_precondition(featuresCalculated,
                "Features must be calculated first.");
        return rugosity;
    }

private:

    bool featuresCalculated;

    int inputArea;
    int convexHullArea;
    double convexity;

    int convexityDefectCount;
    double convexityDefectAreaMean;
    double convexityDefectAreaVariance;
    double convexityDefectAreaSkewness;
    double convexityDefectAreaKurtosis;

    int inputPerimeter;
    int convexHullPerimeter;
    double rugosity;

};

namespace acc {

// template <class ITERATOR, class ACCUMULATOR>
// void extractConvexHullFeatures(ITERATOR start, ITERATOR end, ACCUMULATOR & a)
// {
    // for(unsigned int k=1; k <= a.passesRequired(); ++k)
        // for(ITERATOR i=start; i < end; ++i)
            // a.updatePassN(*i, k);
// }


template <class T1, class S1,
          class T2, class S2,
          class ACCUMULATOR,
          class CONVEX_HULL_ACCUMULATOR>
void extractConvexHullFeatures(MultiArrayView<2, T1, S1> const & data, 
                               MultiArrayView<2, T2, S2> const & labels, 
                               ACCUMULATOR const & regions,
                               CONVEX_HULL_ACCUMULATOR & a)
{
    typedef typename LookupTag<RegionContour, ACCUMULATOR>::value_type Poly;
    typedef typename Poly::value_type Point;
    
    MultiArrayIndex ignored_label = a.ignoredLabel();
    MultiArrayIndex max_label = regions.maxRegionLabel();
    
    a.setMaxRegionLabel(max_label);
    a.ignoreLabel(max_label+1);
    for(MultiArrayIndex k=0; k <= max_label; ++k)
    {
        if(get<Count>(regions, k) == 0 || k == ignored_label)
            continue;
            
        Shape2 start = get<Coord<Minimum> >(regions, k),
               stop  = get<Coord<Maximum> >(regions, k) + Shape2(1);
        Point offset(start);
        a.setCoordinateOffset(k, start);
        
        Poly convex_hull;
        convexHull(get<RegionContour>(regions, k), convex_hull);
        convex_hull -= offset;
        
        MultiArray<2, T2> convex_hull_image(stop-start, max_label+1);
        fillPolygon(convex_hull, convex_hull_image, k);
        
        extractFeatures(data.subimage(start, stop), convex_hull_image, a);
    }
    a.ignoreLabel(ignored_label);
}

template <class T1, class S1,
          class T2, class S2,
          class CONVEX_HULL_ACCUMULATOR>
void extractConvexHullFeatures(MultiArrayView<2, T1, S1> const & data, 
                               MultiArrayView<2, T2, S2> const & labels, 
                               CONVEX_HULL_ACCUMULATOR & a)
{
    AccumulatorChainArray<CoupledArrays<2, T1, T2>,
                          Select<DataArg<1>, LabelArg<2>,
                                 Count, BoundingBox, RegionPerimeter> > regions;
    regions.ignoreLabel(a.ignoredLabel());
    extractFeatures(data, labels, regions);
    extractConvexHullFeatures(data, labels, regions, a);
}

} // namespace acc

} // namespace vigra

#endif // VIGRA_CONVEX_HULL_FEATURES_HXX
