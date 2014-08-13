
#ifndef VIGRA_CONVEXHULL_HXX
#define VIGRA_CONVEXHULL_HXX

#include "multi_array.hxx"
#include "accumulator.hxx"
#include "polygon.hxx"
#include "labelimage.hxx"

namespace vigra {

namespace detail {

static const bool ASCENDING = false;
static const bool DESCENDING = true;

/*
 * functor to order Points by x-coordinate
 */
template<typename T, bool DESCENDING>
class PointXOrdering {
public:
    inline
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        return (DESCENDING) ? !(a[0] < b[0]) : (a[0] < b[0]);
    }
};

/*
 * functor to order Points by y-coordinate
 */
template<typename T, bool DESCENDING>
class PointYOrdering {
public:
    inline
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        return (DESCENDING) ? !(a[1] < b[1]) : (a[1] < b[1]);
    }
};

/*
 * functor to order Points by yx-coordinates (y first, x after)
 */
template<typename T, bool DESCENDING>
class PointYXOrdering {
public:
    inline
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        if (a[1] == b[1])
            return PointXOrdering<T, DESCENDING>()(a, b);
        else
            return PointYOrdering<T, DESCENDING>()(a, b);
    }
};

/*
 * functor to order Points by xy-coordinates (x first, y after)
 */
template<typename T, bool DESCENDING>
class PointXYOrdering {
public:
    inline
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        if (a[0] == b[0])
            return PointYOrdering<T, DESCENDING>()(a, b);
        else
            return PointXOrdering<T, DESCENDING>()(a, b);
    }
};

/*
 * This is used to get the area of an image
 */
template<typename T>
int countNonZero(MultiArray<2, T> const &array) {

    using namespace vigra::acc;

    vigra::MultiArray<2, double> data(array.width(), array.height());
    AccumulatorChainArray<CoupledArrays<2, double, T>,
            Select<LabelArg<2>, Count> > a;
    extractFeatures(data, array, a);

    int nonZero = array.width() * array.height() -  get<Count>(a, 0);
    return nonZero;
}

/*
 * Puts all the points of the line between p1 and p2
 * to the point-vector result.
 * Works either when working with discrete pixel coordinates
 * and when working with floating point coordinates.
 */
template<typename T>
void pushLinePoints(TinyVector<T, 2> const &p1, TinyVector<T, 2> const &p2,
        std::vector<TinyVector<T, 2> > &result) {

    TinyVector<T, 2> ps, diff = p2 - p1;
    bool yLonger = false;
    int incVal, endVal;
    float shortLen = diff[1];
    float longLen = diff[0];

    if (abs(shortLen) > abs(longLen)) {
        std::swap(shortLen, longLen);
        yLonger = true;
    }
    endVal = longLen;
    if (longLen < 0) {
        incVal = -1;
        longLen = -longLen;
    } else
        incVal = 1;
    float decInc;
    if (longLen == 0)
        decInc = 0;
    else
        decInc = shortLen / longLen;
    if (yLonger) {
        int i = incVal;
        float j = decInc;
        for (; abs(i) < abs(endVal); i += incVal, j += decInc)
            result.push_back(TinyVector<T, 2>(p1[0] + j, p1[1] + i));
    } else {
        int i = incVal;
        float j = decInc;
        for (; abs(i) < abs(endVal); i += incVal, j += decInc)
            result.push_back(TinyVector<T, 2>(p1[0] + i, p1[1] + j));
    }

}

/*
 * Scan line filling, works on both convex and non convex polygons
 * The contour points are outside the region to draw
 * Adjacent polygon points must be next to each other in p too.
 */
template<class T>
void fillPolygon(std::vector<TinyVector<float, 2> > const &p,
        MultiArray<2, T> &output_image, T value) {

    std::vector<TinyVector<float, 2> > contour_points;

    std::vector<TinyVector<float, 2> >::const_iterator ip = p.begin();

    for (; ip < p.end() - 1; ++ip) {
        contour_points.push_back(*ip);
        pushLinePoints(*ip, *(ip + 1), contour_points);
    }
    contour_points.push_back(*(p.end() - 1));
    pushLinePoints(p.back(), p.front(), contour_points);

    sort(contour_points.begin(), contour_points.end(),
            PointYXOrdering<float, ASCENDING>());

    std::vector<TinyVector<float, 2> >::iterator points_iterator =
            contour_points.begin();

    float min_y = (contour_points.front())[1];
    float max_y = (contour_points.back())[1];

    while (points_iterator != contour_points.end()) {
        float y = (*points_iterator)[1];
        if ((y > min_y) && (y < max_y)) { // We are inside the polygon
            TinyVector<float, 2> current = *points_iterator;
            float min_x = current[0];
            current[0] = ceil(current[0]);
            current[1] = ceil(current[1]);

            bool drawing = true;

            while ((*points_iterator)[1] == y) {
                ++points_iterator;

                TinyVector<float, 2> endPoint = *points_iterator;
                float max_x = endPoint[0];

                // Draw the scan line

                for (; current[0] < endPoint[0]; ++current[0]) {
                    if ((current[0] > min_x) && (current[0] < max_x)) {
                        if (drawing) {
                            output_image(current[0], current[1]) = value;
                        }
                    }
                }

                drawing = !drawing;
            }

        } else { // We are on an edge, just continue
            ++points_iterator;
        }
    }

}

/*
 * Left hand on the wall contour extraction
 * label_image 1 is object, 0 is background
 * anchor_point is a point inside the object
 * contour_points contains a point for each section of the wall
 */
template<class IMAGE_TYPE>
void extractContour(IMAGE_TYPE const &label_image,
        TinyVector<int, 2> anchor_point,
        std::vector<TinyVector<float, 2> > &contour_points) {

    TinyVector<int, 2> neightbour_coords[4] = { TinyVector<int, 2>(0, -1),
            TinyVector<int, 2>(1, 0), TinyVector<int, 2>(0, 1), TinyVector<int,
                    2>(-1, 0) };

    TinyVector<int, 2> directions[4] = { TinyVector<int, 2>(-1, 0), TinyVector<
            int, 2>(0, -1), TinyVector<int, 2>(1, 0), TinyVector<int, 2>(0, 1) };

    int direction = 0;
    TinyVector<int, 2> position; // Position outside the object from which we place the hand
    TinyVector<int, 2> initial_position;

    // Place the left hand on the wall
    for (int i = 0; i < 4; ++i) {
        TinyVector<int, 2> neighbour = anchor_point + neightbour_coords[i];

        bool neighbour_outside = label_image(neighbour[0], neighbour[1]) == 0;

        if (neighbour_outside) {
            direction = i;
            position = neighbour;
            initial_position = neighbour;
            break;
        }
    }

    // Go around the object
    do {

        // Go forward until we don't touch the wall or there is a right turn
        do {

            TinyVector<int, 2> neighbour = position
                    - neightbour_coords[direction];

            contour_points.push_back((position + neighbour) / 2.0);

            position += directions[direction];

            // We have bumped into a wall
            if (label_image(position[0], position[1])) {
                // Go back and...
                position -= directions[direction];
                // ...Turn right, we should touch the wall again
                direction++;
                if (direction > 3) {
                    direction = 0;
                }
                break;
            }

            neighbour = position - neightbour_coords[direction];

            bool neighbour_outside = label_image(neighbour[0], neighbour[1])
                    == 0;

            // We have lost the wall
            if (neighbour_outside) {
                // Turn left and move forward, we should touch the wall again
                direction--;
                if (direction < 0) {
                    direction = 3;
                }
                position += directions[direction];
                break;
            }

        } while (position != initial_position);

    } while (position != initial_position);

}

/*
 * Returns first pixel found on a border
 */
template<class IMAGE_TYPE>
bool findAnchorPoint(IMAGE_TYPE const &input_image,
        TinyVector<int, 2> &anchor_point) {

    TinyVector<int, 2> neightbour_coords[4] = { TinyVector<int, 2>(0, -1),
            TinyVector<int, 2>(1, 0), TinyVector<int, 2>(0, 1), TinyVector<int,
                    2>(-1, 0) };

    for (unsigned y = 1; y < input_image.height() - 1; ++y) {
        for (unsigned x = 1; x < input_image.width() - 1; ++x) {
            if (input_image(y, x)) {
                for (unsigned i = 0; i < 4; ++i) {
                    unsigned neighbour_y = y + neightbour_coords[i][1];
                    unsigned neighbour_x = x + neightbour_coords[i][0];
                    // It must have at least a background pixel
                    if (!input_image(neighbour_x, neighbour_y)) {
                        anchor_point = TinyVector<int, 2>(x, y);
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

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
        findAnchorPoint(input_image, anchor_point);

        std::vector<TinyVector<float, 2> > contour_points;

        extractContour(input_image, anchor_point, contour_points);

        inputPerimeter = contour_points.size();

        std::vector<TinyVector<float, 2> > convex_hull_points;
        convexHull(contour_points, convex_hull_points);

        convex_hull_points.pop_back(); // Last point is the same as the first one

        MultiArray<2, int> convex_hull_image(
                vigra::Shape2(input_image.height(), input_image.width()));

        fillPolygon(convex_hull_points, convex_hull_image, 1);

        std::vector<TinyVector<float, 2> > convex_hull_contour_points;

        findAnchorPoint(convex_hull_image, anchor_point);
        extractContour(convex_hull_image, anchor_point,
                convex_hull_contour_points);

        convexHullPerimeter = convex_hull_contour_points.size();

        rugosity = (double) inputPerimeter / (double) convexHullPerimeter;

        MultiArray<2, double> diff_image(
                vigra::Shape2(input_image.height(), input_image.width()));

        combineTwoImages(input_image, convex_hull_image, diff_image,
                std::not_equal_to<int>());

        MultiArray<2, int> labels(
                vigra::Shape2(input_image.height(), input_image.width()));

        convexityDefectCount = labelImageWithBackground(diff_image, labels,
                false, 0);

        convexHullArea = countNonZero(convex_hull_image);


        inputArea = countNonZero(input_image);
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

}

}

#endif
