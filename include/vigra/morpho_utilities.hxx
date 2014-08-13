// Author(s): Thomas Walter, Esteban Pardo
// $Date$
// $Rev$
// $URL$

#ifndef MORPHO_UTILITIES_HXX_
#define MORPHO_UTILITIES_HXX_

#include <vector>
#include <iostream>
#include <algorithm>

#include "diff2d.hxx"
#include "inspectimage.hxx"

namespace vigra {
namespace morpho {

/*
 * Morphological functors
 */

template<class VALUETYPE>
class Max {
public:

    /** the functor's argument type
     */
    typedef VALUETYPE argument_type;

    /** the functor's result type
     */
    typedef VALUETYPE result_type;

    /** \deprecated use argument_type
     */
    typedef VALUETYPE value_type;

    void reset() {
        findMinMax.reset();
    }

    VALUETYPE operator()() {
        return findMinMax.max;
    }

    /** update min and max
     */
    void operator()(argument_type const & v) {
        findMinMax(v);
    }

    /** update min and max with components of RGBValue<VALUETYPE>
     */
    void operator()(RGBValue<VALUETYPE> const & v) {
        findMinMax(v);
    }

    /** merge two statistics
     */
    void operator()(Max const & v) {
        this->findMinMax = v.findMinMax;
    }

    FindMinMax<VALUETYPE> findMinMax;

};

template<class VALUETYPE>
class Min {
public:

    /** the functor's argument type
     */
    typedef VALUETYPE argument_type;

    /** the functor's result type
     */
    typedef VALUETYPE result_type;

    /** \deprecated use argument_type
     */
    typedef VALUETYPE value_type;

    void reset() {
        findMinMax.reset();
    }

    VALUETYPE operator()() {
        return findMinMax.min;
    }

    /** update min and max
     */
    void operator()(argument_type const & v) {
        findMinMax(v);
    }

    /** update min and max with components of RGBValue<VALUETYPE>
     */
    void operator()(RGBValue<VALUETYPE> const & v) {
        findMinMax(v);
    }

    /** merge two statistics
     */
    void operator()(Min const & v) {
        this->findMinMax = v.findMinMax;
    }

    FindMinMax<VALUETYPE> findMinMax;

};

template<class VALUETYPE>
class Rank {
public:

    /** the functor's argument type
     */
    typedef VALUETYPE argument_type;

    /** the functor's result type
     */
    typedef VALUETYPE result_type;

    /** \deprecated use argument_type
     */
    typedef VALUETYPE value_type;

    Rank(double r) {
        rank = r;
        pixelCount = 0;
    }

    void reset() {
        std::fill(histogram, histogram + 256, 0);
        pixelCount = 0;
    }

    VALUETYPE operator()() {
        long leftsum = 0;
        unsigned rankPosition = 0;

        if (rank == 0.0) {
            for (unsigned i = 0; i < 256; i++) {
                if (histogram[i]) {
                    rankPosition = i;
                    break;
                }
            }
        } else {
            for (unsigned i = 0; i < 256; i++) {
                if ((float) (histogram[i] + leftsum) / pixelCount >= rank) {
                    rankPosition = i;
                    break;
                }
                leftsum += histogram[i];
            }
        }
        return rankPosition;
    }

    /** update histogram
     */
    void operator()(argument_type const & v) {
        histogram[v]++;
        pixelCount++;
    }

    /** merge two statistics
     */
    void operator()(Rank const & v) {
        rank = v.rank;
        pixelCount = v.pixelCount;
        std::copy(v.histogram, v.histogram + 256, histogram);
    }

    double rank;
    long histogram[256];
    long pixelCount;

};

// HELP FUNCTORS
template<class T, class S>
struct IsGreaterEqual {
    IsGreaterEqual() {
    }
    bool operator()(T const &a, S const &b) {
        return ((a >= b));
    }
};

template<class T, class S>
struct IsSmallerEqual {
    IsSmallerEqual() {
    }
    bool operator()(T const &a, S const &b) {
        return ((a <= b));
    }
};

template<class T, class S>
struct IsGreater {
    IsGreater() {
    }

    bool operator()(T const &a, S const &b) {
        return ((a > b));
    }
};

template<class T, class S>
struct IsSmaller {
    IsSmaller() {
    }
    bool operator()(T const &a, S const &b) {
        return ((a < b));
    }
};

template<class T, class S>
struct IsEqual {
    bool operator()(T const &a, S const &b) {
        return (a == b);
    }
};

template<class T, class S>
struct IsUnequal {
    bool operator()(T const &a, S const &b) {
        return (a != b);
    }
};

// MinFunctor returns the minimum of two numbers.
template<class T>
struct MinFunctor {
    const T neutralValue;
    MinFunctor() :
            neutralValue(vigra::NumericTraits<T>::maxConst) {
    }
    T operator()(T const &a, T const &b) {
        return ((a < b) ? a : b);
    }
};

// MaxFunctor returns the maximum of two numbers.
template<class T>
struct MaxFunctor {
    const T neutralValue;
    MaxFunctor() :
            neutralValue(vigra::NumericTraits<T>::minConst) {
    }
    T operator()(T const &a, T const &b) {
        return ((a > b) ? a : b);
    }
};

// minus and plus with clipping.
template<class T>
struct minusConstantClipp {
    typedef typename vigra::NumericTraits<T> NumTraits;
    typedef typename vigra::NumericTraits<T>::Promote SUMTYPE;

    T c;

    minusConstantClipp(T const & val) :
            c(val) {
    }

    T operator()(T const & a) const {
        SUMTYPE sum = NumTraits::toPromote(a) - NumTraits::toPromote(c);
        return (NumTraits::fromPromote(sum));
    }
};

template<class T>
struct plusConstantClipp {
    typedef typename vigra::NumericTraits<T> NumTraits;
    typedef typename vigra::NumericTraits<T>::Promote SUMTYPE;

    T c;

    plusConstantClipp(T const & val) :
            c(val) {
    }

    T operator()(T const & a) const {
        SUMTYPE sum = NumTraits::toPromote(c) + NumTraits::toPromote(a);
        return (NumTraits::fromPromote(sum));
    }
};

// The Pixel2D structure is used for priority queues (reconstruction, watershed & co)
template<class ValueType>
class Pixel2D {
public:
    Pixel2D() :
            value(0), offset(vigra::Diff2D(0, 0)), insertionOrder(0) {
    }

    Pixel2D(const ValueType & val, const vigra::Diff2D & loc,
            unsigned long insOrder = 0) :
            value(val), offset(loc), insertionOrder(insOrder) {
    }

    bool operator>(const Pixel2D & p1) const {
        bool res = (value == p1.value);
        if (res)
            return (insertionOrder > p1.insertionOrder);
        else
            return (value > p1.value);
    }

    bool operator<(const Pixel2D & p1) const {
        bool res = (value == p1.value);
        if (res)
            return (insertionOrder < p1.insertionOrder);
        else
            return (value < p1.value);
    }

    ValueType value;
    vigra::Diff2D offset;
    unsigned long insertionOrder;
};

// Priority for algorithms starting with local maxima.
template<class T>
struct PriorityTopDown {
    inline bool operator()(const Pixel2D<T> & p1, const Pixel2D<T> &p2) {
        return p1 < p2;
    }
};

// Priority for algorithms starting with local minima.
template<class T>
struct PriorityBottomUp {
    inline bool operator()(const Pixel2D<T> & p1, const Pixel2D<T> &p2) {
        return p1 > p2;
    }
};

using vigra::Diff2D;

// The following definitions are used by the class neighborhood2D
const Diff2D NB8_WC[9] = { Diff2D(0, 0), Diff2D(1, 0), Diff2D(1, 1), Diff2D(0,
        1), Diff2D(-1, 1), Diff2D(-1, 0), Diff2D(-1, -1), Diff2D(0, -1), Diff2D(
        1, -1) };

const Diff2D SEG_Y[3] = { Diff2D(0, -1), Diff2D(0, 0), Diff2D(0, 1) };

const Diff2D SEG_X[3] = { Diff2D(-1, 0), Diff2D(0, 0), Diff2D(1, 0) };

const Diff2D NB4_WC[5] = { Diff2D(0, 0), Diff2D(1, 0), Diff2D(0, 1), Diff2D(-1,
        0), Diff2D(0, -1) };

const Diff2D NB8[8] = { Diff2D(1, 0), Diff2D(1, 1), Diff2D(0, 1), Diff2D(-1, 1),
        Diff2D(-1, 0), Diff2D(-1, -1), Diff2D(0, -1), Diff2D(1, -1) };

const Diff2D NB4[4] =
        { Diff2D(1, 0), Diff2D(0, 1), Diff2D(-1, 0), Diff2D(0, -1) };

class NeighborDefinition {
public:
    const Diff2D *nbList;
    const unsigned int nbPixels;
    NeighborDefinition(const Diff2D *l, unsigned int n) :
            nbList(l), nbPixels(n) {
    }
};

// 2D neighborhoods for 4 and 8 connection
const NeighborDefinition WITHCENTER8(NB8_WC, 9);
const NeighborDefinition WITHCENTER4(NB4_WC, 5);
const NeighborDefinition WITHOUTCENTER8(NB8, 8);
const NeighborDefinition WITHOUTCENTER4(NB4, 4);

// neighbors in horizontal and vertical direction
const NeighborDefinition XSEGMENT(SEG_X, 3);
const NeighborDefinition YSEGMENT(SEG_Y, 3);

// neighborPixels handles the list of 2D neighbors.
// it also calculates the maximal extension of the neighborhood in x- and y-direction
// and
class neighborPixels {
protected:
    std::vector<Diff2D> support;
    unsigned long nbPixels_;

    // for border treatment
    Diff2D minOffset_;
    Diff2D maxOffset_;

public:
    typedef std::vector<Diff2D>::iterator ITERATORTYPE;
    typedef std::vector<Diff2D>::size_type SIZETYPE;

    // Constructors:
    neighborPixels() :
            nbPixels_(0) {
    }

    neighborPixels(std::vector<Diff2D> supportParam) :
            support(supportParam), nbPixels_(supportParam.size()) {
        CalculateExtension();
    }

    neighborPixels(const Diff2D *beg, const Diff2D *end) :
            support(beg, end), nbPixels_(end - beg) {
        CalculateExtension();
    }

    neighborPixels(const Diff2D *beg, int nbPixels) :
            support(beg, beg + nbPixels), nbPixels_(nbPixels) {
        CalculateExtension();
    }

    neighborPixels(const NeighborDefinition &nd) :
            support(nd.nbList, nd.nbList + nd.nbPixels), nbPixels_(nd.nbPixels) {
        CalculateExtension();
    }

    Diff2D minOffset() {
        return (minOffset_);
    }

    Diff2D maxOffset() {
        return (maxOffset_);
    }

    // Calculate the maximal extensions of the structuring element.
    void CalculateExtension() {
        if (!support.empty()) {
            minOffset_ = *(support.begin());
            maxOffset_ = *(support.begin());

            for (ITERATORTYPE iter = support.begin(); iter != support.end();
                    ++iter) {
                minOffset_.x = std::min(minOffset_.x, (*iter).x);
                minOffset_.y = std::min(minOffset_.y, (*iter).y);
                maxOffset_.x = std::max(maxOffset_.x, (*iter).x);
                maxOffset_.y = std::max(maxOffset_.y, (*iter).y);
            }
        }
    }

    ITERATORTYPE begin() {
        return (support.begin());
    }

    ITERATORTYPE end() {
        return (support.end());
    }

    void output() {

        std::cout << "Coordinates of the structuring elements: " << std::endl;
        for (std::vector<Diff2D>::iterator iter = support.begin();
                iter != support.end(); ++iter) {
            std::cout << "(" << (*iter).x << ", " << (*iter).y << ")  ";
        }
        std::cout << std::endl;
        std::cout << "Maximal Extensions of the SE: " << std::endl;
        std::cout << "maxOffset: " << maxOffset_.x << " " << maxOffset_.y
                << std::endl;
        std::cout << "minOffset:  " << minOffset_.x << " " << minOffset_.y
                << std::endl;
        return;
    }

    unsigned long numberOfPixels() {
        return (nbPixels_);
    }

    const Diff2D &operator[](int i) {
        return (support[i]);
    }

};

// The class structuringElement2D inherits from neighborPixels
// and contains an additional size parameter. This parameter is used in erosions and
// dilations for iteration.
// The point (0,0) is assumed to be the center of the structuring element.
// The center of the structuring element does not have to belong to it.
class structuringElement2D: public neighborPixels {

public:

    // default constructor
    structuringElement2D() :
            size(0) {
    }

//    structuringElement2D(std::vector<int[2]> supportParam, int sizeParam = 1):
//                         size(sizeParam)
//    {
//        for(std::vector<int[2]>::iterator iter = supportParam.begin();
//            iter != supportParam.end();
//            ++iter)
//        {
//            support.push_back(Diff2D((*iter)[0], (*iter)[1]));
//        }
//    }

//    structuringElement2D(std::vector<float[2]> supportParam, int sizeParam = 1):
//                         size(sizeParam)
//    {
//        for(std::vector<int[2]>::iterator iter = supportParam.begin();
//            iter != supportParam.end();
//            ++iter)
//        {
//            support.push_back(Diff2D((int)(*iter)[0], (int)(*iter)[1]));
//        }
//    }

    structuringElement2D(std::vector<Diff2D> supportParam, int sizeParam = 1) :
            neighborPixels(supportParam), size(sizeParam) {
    }

    structuringElement2D(const Diff2D *beg, const Diff2D *end,
            int sizeParam = 1) :
            neighborPixels(beg, end), size(sizeParam) {
    }

    structuringElement2D(const Diff2D *beg, int nbPixels, int sizeParam = 1) :
            neighborPixels(beg, nbPixels), size(sizeParam) {
    }

    structuringElement2D(const NeighborDefinition &nd, int sizeParam = 1) :
            neighborPixels(nd), size(sizeParam) {
    }

    void transpose() {
        for (std::vector<Diff2D>::iterator iter = support.begin();
                iter != support.end(); ++iter) {
            (*iter).x = (-1) * (*iter).x;
            (*iter).y = (-1) * (*iter).y;
        }
        CalculateExtension();
    }

    int size;
};

// The class neighborhood2D inherits from class neighborPixels
// and has also its border treatment.
class neighborhood2D: public neighborPixels {
public:

    neighborhood2D(std::vector<Diff2D> supportParam, Diff2D imageSize) :
            neighborPixels(supportParam), imageSize_(imageSize) {
    }

    neighborhood2D(const Diff2D *beg, const Diff2D *end,
            const Diff2D &imageSize) :
            neighborPixels(beg, end), imageSize_(imageSize) {
    }

    neighborhood2D(const Diff2D *beg, int nbPixels, const Diff2D &imageSize) :
            neighborPixels(beg, nbPixels), imageSize_(imageSize) {
    }

    neighborhood2D(const NeighborDefinition &nd, const Diff2D &imageSize) :
            neighborPixels(nd), imageSize_(imageSize) {
    }

    neighborhood2D(const NeighborDefinition &nd, const Size2D &imageSize) :
            neighborPixels(nd), imageSize_(imageSize) {
    }

    neighborhood2D(const NeighborDefinition &nd, Size2D &imageSize) :
            neighborPixels(nd), imageSize_(imageSize) {
    }

    bool isBorderPixel(const Diff2D &pixel) {
        return ((pixel.x == 0) || (pixel.x == imageSize_.x - 1)
                || (pixel.y == 0) || (pixel.y == imageSize_.y - 1));
    }

    bool isOutsidePixel(const Diff2D &pixel, const Diff2D &offset) {
        return (((pixel + offset).x < 0)
                || ((pixel + offset).x > imageSize_.x - 1)
                || ((pixel + offset).y < 0)
                || ((pixel + offset).y > imageSize_.y - 1));
    }

    bool isOutsidePixel(const Diff2D &pixel) {
        return ((pixel.x < 0) || (pixel.x > imageSize_.x - 1) || (pixel.y < 0)
                || (pixel.y > imageSize_.y - 1));
    }
private:
    //Diff2D imageSize_;
    Size2D imageSize_;
};

/*
 * Creates a disc structuring element
 */
void generateDiscSE(int radius, structuringElement2D &SE) {

    std::vector<Diff2D> strElCoordinates;

    double r2 = (double) radius * radius;

    // Go through each of the disk sections
    for (int y = -radius; y <= radius; ++y) {
        double sectionWidth = VIGRA_CSTD::sqrt(r2 - y * y);
        for (int x = -sectionWidth; x <= sectionWidth; ++x) {
            strElCoordinates.push_back(Diff2D(x, y));
        }
    }

    SE = structuringElement2D(strElCoordinates);
}

}
;
}
;

#endif /*MORPHO_UTILITIES_HXX_*/

