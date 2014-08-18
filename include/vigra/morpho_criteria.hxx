/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2014 by                                 */
/*               Ullrich Koethe,                                        */
/*               Thomas Walter,                                         */
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

#ifndef MORPHO_CRITERIA_HXX_
#define MORPHO_CRITERIA_HXX_

#include "morpho_utilities.hxx"
#include "morpho_geodesy.hxx"

namespace vigra {
namespace morpho {

const int LAB_QUEUED = -1;
const int LAB_NOT_PROCESSED = 0;

template<class PointType, class ValType>
class LakeArea {
public:
    LakeArea() :
            area_(0) {
    }

    LakeArea(const int& areaParam) :
            area_(areaParam) {
    }

    LakeArea(const LakeArea& areaParam) :
            area_(areaParam.area()) {
    }

    int area() const {
        return (area_);
    }

    void actualize(const PointType & point, const ValType & val) {
        ++area_;
    }

    LakeArea & operator=(const int & a) {
        area_ = a;
        return *this;
    }

    LakeArea & operator=(const LakeArea &a) {
        area_ = a.area();
        return *this;
    }

    bool operator==(const LakeArea & a) {
        return (area_ == a.area());
    }

    bool operator!=(const LakeArea & a) {
        return (area_ != a.area());
    }

    bool operator>(const LakeArea & a) {
        return (area_ > a.area());
    }

    bool operator>=(const LakeArea & a) {
        return (area_ >= a.area());
    }

    bool operator<(const LakeArea & a) {
        return (area_ < a.area());
    }

    bool operator<=(const LakeArea & a) {
        return (area_ <= a.area());
    }

    //
    bool operator==(const int & a) {
        return (area_ == a);
    }

    bool operator!=(const int & a) {
        return (area_ != a);
    }

    bool operator>(const int & a) {
        return (area_ > a);
    }

    bool operator>=(const int & a) {
        return (area_ >= a);
    }

    bool operator<(const int & a) {
        return (area_ < a);
    }

    bool operator<=(const int & a) {
        return (area_ <= a);
    }

    //
    LakeArea& operator++() {
        ++area_;
        return *this;
    }

    LakeArea operator++(int) {
        LakeArea ret(*this);
        ++*this;
        return (ret);
    }

    LakeArea & operator+=(const LakeArea &a) {
        area_ += a.area();
        return *this;
    }

    LakeArea & operator+=(const int &a) {
        area_ += a;
        return *this;
    }

private:
    int area_;

};

/////////////////////////////////
// morphoDiameterOpening/morphoDiameterClosing
template<class PointType, class ValType>
class LakeDiameter {
public:
    LakeDiameter() :
            minCoordinates_(vigra::Diff2D(2147483647, 2147483647)), maxCoordinates_(
                    vigra::Diff2D(0, 0)) {
    }

    LakeDiameter(vigra::Diff2D minParam, vigra::Diff2D maxParam) :
            minCoordinates_(minParam), maxCoordinates_(maxParam) {
    }

    vigra::Diff2D MinCoord() const {
        return (minCoordinates_);
    }
    vigra::Diff2D MaxCoord() const {
        return (maxCoordinates_);
    }

    void actualize(const PointType & point, const ValType & val) {
        minCoordinates_.x = std::min(minCoordinates_.x, point.x);
        minCoordinates_.y = std::min(minCoordinates_.y, point.y);
        maxCoordinates_.x = std::max(maxCoordinates_.x, point.x);
        maxCoordinates_.y = std::max(maxCoordinates_.y, point.y);
    }

    int Extension() const {
        return (std::max(maxCoordinates_.x - minCoordinates_.x,
                maxCoordinates_.y - minCoordinates_.y) + 1);
    }

    LakeDiameter & operator+=(const LakeDiameter &a) {
        vigra::Diff2D minC(a.MinCoord());
        vigra::Diff2D maxC(a.MaxCoord());

        minCoordinates_.x = std::min(minCoordinates_.x, minC.x);
        minCoordinates_.y = std::min(minCoordinates_.y, minC.y);
        maxCoordinates_.x = std::max(maxCoordinates_.x, maxC.x);
        maxCoordinates_.y = std::max(maxCoordinates_.y, maxC.y);

        return *this;
    }

    LakeDiameter & operator+=(const vigra::Diff2D &a) {
        minCoordinates_.x = std::min(minCoordinates_.x, a.x);
        minCoordinates_.y = std::min(minCoordinates_.y, a.y);
        maxCoordinates_.x = std::max(maxCoordinates_.x, a.x);
        maxCoordinates_.y = std::max(maxCoordinates_.y, a.y);

        return *this;
    }

    bool operator>(const int & a) const {
        return (Extension() > a);
    }

    bool operator<(const int & a) const {
        return (Extension() < a);
    }

    bool operator>=(const int & a) const {
        return (Extension() >= a);
    }

    bool operator<=(const int & a) const {
        return (Extension() <= a);
    }

    bool operator==(const int & a) const {
        return (Extension() == a);
    }

    bool operator!=(const int & a) const {
        return (Extension() != a);
    }

    bool operator>(const LakeDiameter & a) const {
        return (Extension() > a.Extension());
    }

    bool operator<(const LakeDiameter & a) const {
        return (Extension() < a.Extension());
    }

    bool operator>=(const LakeDiameter & a) const {
        return (Extension() >= a.Extension());
    }

    bool operator<=(const LakeDiameter & a) const {
        return (Extension() <= a.Extension());
    }

    bool operator==(const LakeDiameter & a) const {
        return (Extension() == a.Extension());
    }

    bool operator!=(const LakeDiameter & a) const {
        return (Extension() != a.Extension());
    }

private:
    vigra::Diff2D minCoordinates_;
    vigra::Diff2D maxCoordinates_;

};

/*
 * Floods an intensity image using the pixel intensity as height.
 * This algorithm is one step of the watershed algorithm
 */
template<class AttributeStructure, class Iterator1, class Accessor1,
        class Iterator2, class Accessor2, class NBTYPE, class MinmaxFunctor,
        class PriorityFunctor, class StopLevelInitFunctor,
        class LevelInitFunctor>
void morphoMakeFlooding(Iterator1 srcUpperLeft, Iterator1 srcLowerRight,
        Accessor1 srca, Iterator2 destUpperLeft, Accessor2 desta, int areaMax,
        NBTYPE & nbOffset, MinmaxFunctor minmax, PriorityFunctor priority,
        StopLevelInitFunctor stopLevelInit, LevelInitFunctor levelInit) {

    typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
    typedef typename NBTYPE::SIZETYPE SIZETYPE;
    typedef typename Accessor1::value_type VALUETYPE;

    typedef Pixel2D<VALUETYPE> PIX;
    typedef int LABTYPE;

    // Settings for areaclose
    // 1.) priority is PriorityBottomUp
    // 2.) minmax is IsSmaller (finding the minima in utilityMinMaxLabel)
    // 3.) stopLevel (initVal) is calculated as max (stop levels are set
    //      to the image maximum.
    // 4.) level_before is min (the starting level is set to the image
    //      minimum.

//    std::priority_queue<PIX, std::vector<PIX>, PriorFunctor> PQ(priority);
//    typedef PriorityBottomUp<typename Accessor1::value_type> PriorFunctor;
//    PriorFunctor priority;

    int width = srcLowerRight.x - srcUpperLeft.x;
    int height = srcLowerRight.y - srcUpperLeft.y;

    vigra::BasicImage<int> labelImage(width, height), labelImageFinal(width,
            height);
    vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft(),
            labFinalUpperLeft = labelImageFinal.upperLeft();
    vigra::BasicImage<int>::Accessor lab;

    int numberOfMinima = utilityMinMaxLabel(srcUpperLeft, srcLowerRight, srca,
            labUpperLeft, lab, minmax, nbOffset);
    std::vector<AttributeStructure> area(numberOfMinima + 1);

    // equivalence takes the label of the lake with which it has been fused.
    // at the moment, this is simply i.
    std::vector<int> equivalence(numberOfMinima + 1);
    for (std::vector<int>::size_type i = 0; i < equivalence.size(); ++i)
        equivalence[i] = i;

    // the stop levels
    std::vector<VALUETYPE> stopLevel(numberOfMinima + 1);

    // to take the values of the minma
    std::vector<VALUETYPE> valOfMin(numberOfMinima + 1);

    std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);

    Diff2D o0(0, 0);

    // levelBefore and initVal are set to the value of the first pixel encountered.
    VALUETYPE levelBefore = srca(srcUpperLeft, o0);
    VALUETYPE initVal = srca(srcUpperLeft, o0);

    unsigned long insertionOrder = 0;

    // initialization of the hierarchical queue
    for (o0.y = 0; o0.y < height; ++o0.y) {
        for (o0.x = 0; o0.x < width; ++o0.x) {
            VALUETYPE val = srca(srcUpperLeft, o0);
            LABTYPE label = lab(labUpperLeft, o0);

            // labelFinal is just a copy of label in the beginning.
            lab.set(label, labFinalUpperLeft, o0);
            if (label > LAB_NOT_PROCESSED) {
                // the area of the minimum is incremented.
                area[label].actualize(o0, val);
                valOfMin[label] = val;

                // we determine levelBefore as the absolute minimum of the input image.
                levelBefore = levelInit(val, levelBefore);

                // the output image is set to 1
                desta.set(1, destUpperLeft, o0);

                // look to the neighborhood.
                for (ITERATORTYPE iter = nbOffset.begin();
                        iter != nbOffset.end(); ++iter) {
                    Diff2D o1 = o0 + (*iter);
                    // if the neighbor is not outside the image
                    // and if it has no label and if it is not in the queue
                    if ((!nbOffset.isOutsidePixel(o1))
                            && (lab(labUpperLeft, o1) == LAB_NOT_PROCESSED)) {
                        PQ.push(
                                PIX(srca(srcUpperLeft, o1), o1,
                                        insertionOrder++));
                        lab.set(LAB_QUEUED, labUpperLeft, o1);
                    }
                } // end for neighborhood
            } // end if label
            else {
                initVal = stopLevelInit(initVal, val);
            } // end else (not label)
        } // end x-loop
    } // end y-loop

    for (std::vector<int>::size_type i = 1; i < area.size(); ++i) {
        // stop_levels must be set if the area of the lake
        // already exceeds areaMax.
        if (area[i] >= areaMax)
            stopLevel[i] = valOfMin[i];
        else
            stopLevel[i] = initVal;
    }

    while (!PQ.empty()) {
        PIX px = PQ.top();
        PQ.pop();
        VALUETYPE level = px.value;
        Diff2D o0 = px.offset;

        if (levelBefore != level) {
            // in this case, the stop_levels have to be determined
            for (std::vector<int>::size_type i = 1; i < area.size(); ++i) {
                // first we find the dominating lake, i.e. the label
                // that has been propagated when two lakes were fused.
                int lab = i;
                while (lab != equivalence[lab])
                    lab = equivalence[lab];

                // stop_levels must be set if the area of the lake
                // already exceeds areaMax.
                if ((area[lab] >= areaMax) && (stopLevel[i] == initVal))
                    stopLevel[i] = (VALUETYPE) levelBefore;

            }
            levelBefore = level;
        } // end of level change

        // normal flooding procedure
        int label1 = 0;
        int label2 = 0;

        // look to the neighborhood to determine the label of pixel o0.
        for (ITERATORTYPE iter = nbOffset.begin(); iter != nbOffset.end();
                ++iter) {
            Diff2D o1 = o0 + *iter;
            if (!nbOffset.isOutsidePixel(o1)) {
                LABTYPE label_o1 = lab(labUpperLeft, o1);

                // first case: pixel has not been processed.
                if (label_o1 == LAB_NOT_PROCESSED) {
                    PQ.push(PIX(srca(srcUpperLeft, o1), o1, insertionOrder++));
                    lab.set(LAB_QUEUED, labUpperLeft, o1);
                }

                // second case: neighbor pixel is already in the queue:
                // nothing is to be done, then.

                // third case: the neighbor has a label
                if (label_o1 > LAB_NOT_PROCESSED) {
                    label2 = label_o1;
                    while (label2 != equivalence[label2])
                        label2 = equivalence[label2];

                    if (label1 == 0) {
                        // in this case, the label is the first
                        // which has been found in the neighborhood.
                        label1 = label2;

                        lab.set(label1, labUpperLeft, o0);
                        if (stopLevel[label1] == initVal)
                            lab.set(label1, labFinalUpperLeft, o0);

                        area[label1].actualize(o0, srca(srcUpperLeft, o0));
                    } else {
                        // in this case, a label has already been assigned to o0.
                        if (label1 != label2) {
                            // in this case, we have a meeting point of two lakes.
                            // we therefore have to fuse the two lakes.
                            if (area[label1] > area[label2]) {
                                area[label1] += area[label2];
                                equivalence[label2] = label1;
                            } else {
                                area[label2] += area[label1];
                                equivalence[label1] = label2;
                                label1 = label2;
                            }

                        }
                    }

                }
            }
        } // end for neighborhood

    } // end of PRIORITY QUEUE

    // final assignment
    for (o0.y = 0; o0.y < height; ++o0.y) {
        for (o0.x = 0; o0.x < width; ++o0.x) {
            LABTYPE label = lab(labFinalUpperLeft, o0);
            if (label > LAB_NOT_PROCESSED) {
                desta.set(stopLevel[label], destUpperLeft, o0);
            } else
                desta.set(srca(srcUpperLeft, o0), destUpperLeft, o0);
        } // end of x-loop
    } // end of y-loop

} // end of function

/////////////////
// Area closing
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoAreaClosing(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::pair<Iterator2, Accessor2> dest, int areaMax,
        NBTYPE & neighborOffset) {

    typedef typename Accessor1::value_type val_type;

    morphoMakeFlooding<LakeArea<vigra::Diff2D, val_type> >(src.first,
            src.second, src.third, dest.first, dest.second, areaMax,
            neighborOffset,
            IsSmaller<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityBottomUp<val_type>(),
            MaxFunctor<typename Accessor1::value_type>(),
            MinFunctor<typename Accessor1::value_type>());

}

template<class Image1, class Image2, class NBTYPE>
void morphoAreaClosing(const Image1 & imin, Image2 & imout, int areaMax,
        NBTYPE & nbOffset) {
    morphoAreaClosing(srcImageRange(imin), destImage(imout), areaMax, nbOffset);
}

/////////////////
// Area opening
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoAreaOpening(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::pair<Iterator2, Accessor2> dest, int areaMax,
        NBTYPE & neighborOffset) {
    typedef typename Accessor1::value_type val_type;
    morphoMakeFlooding<LakeArea<vigra::Diff2D, val_type> >(src.first,
            src.second, src.third, dest.first, dest.second, areaMax,
            neighborOffset,
            IsGreater<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityTopDown<val_type>(),
            MinFunctor<typename Accessor1::value_type>(),
            MaxFunctor<typename Accessor1::value_type>());
}

template<class Image1, class Image2, class NBTYPE>
void morphoAreaOpening(const Image1 & imin, Image2 & imout, int areaMax,
        NBTYPE & nbOffset) {
    morphoAreaOpening(srcImageRange(imin), destImage(imout), areaMax, nbOffset);
}

///////////////////
// Diameter closing
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoDiameterClosing(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::pair<Iterator2, Accessor2> dest, int areaMax,
        NBTYPE & neighborOffset) {
    typedef typename Accessor1::value_type val_type;

    morphoMakeFlooding<LakeDiameter<vigra::Diff2D, val_type> >(src.first,
            src.second, src.third, dest.first, dest.second, areaMax,
            neighborOffset,
            IsSmaller<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityBottomUp<typename Accessor1::value_type>(),
            MaxFunctor<typename Accessor1::value_type>(),
            MinFunctor<typename Accessor1::value_type>());
}

template<class Image1, class Image2, class NBTYPE>
void morphoDiameterClosing(const Image1 & imin, Image2 & imout, int areaMax,
        NBTYPE & nbOffset) {
    morphoDiameterClosing(srcImageRange(imin), destImage(imout), areaMax,
            nbOffset);
}

///////////////////
// Diameter opening
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoDiameterOpening(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::pair<Iterator2, Accessor2> dest, int areaMax,
        NBTYPE & neighborOffset) {
    typedef typename Accessor1::value_type val_type;
    morphoMakeFlooding<LakeDiameter<vigra::Diff2D, val_type> >(src.first,
            src.second, src.third, dest.first, dest.second, areaMax,
            neighborOffset,
            IsGreater<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityTopDown<val_type>(),
            MinFunctor<typename Accessor1::value_type>(),
            MaxFunctor<typename Accessor1::value_type>());
}

template<class Image1, class Image2, class NBTYPE>
void morphoDiameterOpening(const Image1 & imin, Image2 & imout, int areaMax,
        NBTYPE & nbOffset) {
    morphoDiameterOpening(srcImageRange(imin), destImage(imout), areaMax,
            nbOffset);
}

}
;
}
;

#endif /* MORPHO_CRITERIA_HXX_ */


