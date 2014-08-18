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

#ifndef MORPHO_DYNAMIC_HXX_
#define MORPHO_DYNAMIC_HXX_

#include "morpho_utilities.hxx"
#include "morpho_geodesy.hxx"

namespace vigra {
namespace morpho {

///////////////////
// ExtinctionValues
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE, class MinmaxFunctor, class PriorityFunctor>
void morphoExtinctionValues(Iterator1 srcUpperLeft, Iterator1 srcLowerRight,
        Accessor1 srca, Iterator2 destUpperLeft, Accessor2 desta,
        std::vector<typename Accessor1::value_type> & dynamics,
        NBTYPE & nbOffset, MinmaxFunctor minmax, PriorityFunctor priority) {

    typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
    typedef typename NBTYPE::SIZETYPE SIZETYPE;
    typedef typename Accessor1::value_type VALUETYPE;

    typedef Pixel2D<VALUETYPE> PIX;
    typedef int LABTYPE;

    // Settings for morphoDynMinima
    // 1.) priority is PriorityBottomUp
    // 2.) minmax is IsSmaller (finding the minima in morphoMinMaxLabel)

    int width = srcLowerRight.x - srcUpperLeft.x;
    int height = srcLowerRight.y - srcUpperLeft.y;

    vigra::BasicImage<int> labelImage(width, height);
    vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft();
    vigra::BasicImage<int>::Accessor lab;

    int numberOfMinima = utilityMinMaxLabel(srcUpperLeft, srcLowerRight, srca,
            labUpperLeft, lab, minmax, nbOffset);

    // equivalence takes the label of the lake with which it has been fused.
    // at the moment, this is simply i.
    std::vector<int> equivalence(numberOfMinima + 1);
    for (std::vector<int>::size_type i = 0; i < equivalence.size(); ++i) {
        equivalence[i] = i;
    }

    // the vector containing the dynamics
    //std::vector<VALUETYPE> dynamics(numberOfMinima + 1);
    for (int i = 0; i < numberOfMinima + 1; ++i) {
        dynamics.push_back(0);
    }

    VALUETYPE maxDyn = 0;
    LABTYPE labelOfMaxDyn = -1;

    // as second criterion, we take the area
    std::vector<unsigned> area(numberOfMinima + 1);

    // to take the values of the minma
    std::vector<VALUETYPE> valOfMin(numberOfMinima + 1);

    // Priority queue
    std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);

    Diff2D o0(0, 0);

    unsigned long insertionOrder = 0;

    // initialization of the hierarchical queue
    for (o0.y = 0; o0.y < height; ++o0.y) {
        for (o0.x = 0; o0.x < width; ++o0.x) {
            VALUETYPE val = srca(srcUpperLeft, o0);
            LABTYPE label = lab(labUpperLeft, o0);

            if (label > LAB_NOT_PROCESSED) {
                desta.set(1, destUpperLeft, o0);

                // the area of the minimum is incremented.
                area[label]++;

                // the value of the minimum is assigned.
                valOfMin[label] = val;

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
                desta.set(0, destUpperLeft, o0);
            }
        } // end x-loop
    } // end y-loop

    // Entering the priority queue ...
    while (!PQ.empty()) {
        PIX px = PQ.top();
        PQ.pop();
        VALUETYPE level = px.value;
        Diff2D o0 = px.offset;

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
                        area[label1] = area[label1] + 1;
                    } else {
                        // in this case, a label has already been assigned to o0.
                        if (label1 != label2) {
                            // in this case, we have a meeting point of two lakes.
                            // we can therefore determine the dynamic
                            if (minmax(valOfMin[label1], valOfMin[label2])
                                    || ((valOfMin[label1] == valOfMin[label2])
                                            && (area[label1] > area[label2]))) {
                                dynamics[label2] = level - valOfMin[label2];
                                if (dynamics[label2] > maxDyn) {
                                    maxDyn = dynamics[label2];
                                    labelOfMaxDyn = label2;
                                }
                                area[label1] += area[label2];
                                equivalence[label2] = label1;
                            } else {
                                dynamics[label1] = level - valOfMin[label1];
                                if (dynamics[label1] > maxDyn) {
                                    maxDyn = dynamics[label1];
                                    labelOfMaxDyn = label1;
                                }
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

    for (int i = 1; i < numberOfMinima + 1; ++i) {
        if (dynamics[i] == 0) {
            dynamics[i] = maxDyn + valOfMin[labelOfMaxDyn] - valOfMin[i];
            break;
        }
    }

    // initialization of the hierarchical queue
    for (o0.y = 0; o0.y < height; ++o0.y) {
        for (o0.x = 0; o0.x < width; ++o0.x) {
            VALUETYPE val = desta(destUpperLeft, o0);
            if (val > 0) {
                desta.set(dynamics[lab(labUpperLeft, o0)], destUpperLeft, o0);
            }
        }
    }

} // end of function

/////////////////////////////
// ExtinctionValues with mask
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE, class MinmaxFunctor, class PriorityFunctor>
void morphoMaskExtinctionValues(Iterator1 srcUpperLeft, Iterator1 srcLowerRight,
        Accessor1 srca, Iterator2 maskUpperLeft, Accessor2 maska,
        std::vector<typename Accessor1::value_type> & dynamics,
        std::vector<typename Accessor1::value_type> & valOfMin,
        NBTYPE & nbOffset, MinmaxFunctor minmax, PriorityFunctor priority,
        typename Accessor2::value_type maskLabel) {

    typedef typename NBTYPE::ITERATORTYPE ITERATORTYPE;
    typedef typename NBTYPE::SIZETYPE SIZETYPE;
    typedef typename Accessor1::value_type VALUETYPE;

    typedef Pixel2D<VALUETYPE> PIX;
    typedef int LABTYPE;

    // Settings for morphoDynMinima
    // 1.) priority is PriorityBottomUp
    // 2.) minmax is IsSmaller (finding the minima in morphoMinMaxLabel)

    int width = srcLowerRight.x - srcUpperLeft.x;
    int height = srcLowerRight.y - srcUpperLeft.y;

    vigra::BasicImage<int> labelImage(width, height);
    vigra::BasicImage<int>::Iterator labUpperLeft = labelImage.upperLeft();
    vigra::BasicImage<int>::Accessor lab;

    int numberOfMinima = utilityMinMaxLabel(srcUpperLeft, srcLowerRight, srca,
            labUpperLeft, lab, maskUpperLeft, maska, minmax, nbOffset,
            maskLabel);

    // equivalence takes the label of the lake with which it has been fused.
    // at the moment, this is simply i.
    std::vector<int> equivalence(numberOfMinima + 1);
    for (std::vector<int>::size_type i = 0; i < equivalence.size(); ++i) {
        equivalence[i] = i;
    }

    // the vector containing the dynamics
    //std::vector<VALUETYPE> dynamics(numberOfMinima + 1);
    for (int i = 0; i < numberOfMinima + 1; ++i) {
        dynamics.push_back(0);
    }
    VALUETYPE maxVal = 0;
    VALUETYPE minVal = 0;

    // as second criterion, we take the area
    std::vector<unsigned> area(numberOfMinima + 1);

    // to take the values of the minma
    //std::vector<VALUETYPE> valOfMin(numberOfMinima + 1);
    for (int i = 0; i < numberOfMinima + 1; ++i) {
        valOfMin.push_back(0);
    }

    // Priority queue
    std::priority_queue<PIX, std::vector<PIX>, PriorityFunctor> PQ(priority);

    Diff2D o0(0, 0);

    unsigned long insertionOrder = 0;

    // initialization of the hierarchical queue
    for (o0.y = 0; o0.y < height; ++o0.y) {
        for (o0.x = 0; o0.x < width; ++o0.x) {
            VALUETYPE val = srca(srcUpperLeft, o0);
            LABTYPE label = lab(labUpperLeft, o0);

            if (val > maxVal)
                maxVal = val;

            if (label > LAB_NOT_PROCESSED) {
                // the area of the minimum is incremented.
                area[label]++;

                // the value of the minimum is assigned.
                valOfMin[label] = val;

                // look to the neighborhood.
                for (ITERATORTYPE iter = nbOffset.begin();
                        iter != nbOffset.end(); ++iter) {
                    Diff2D o1 = o0 + (*iter);
                    // if the neighbor is not outside the image
                    // and if it has no label and if it is not in the queue
                    if ((!nbOffset.isOutsidePixel(o1))
                            && (maska(maskUpperLeft, o1) == maskLabel)
                            && (lab(labUpperLeft, o1) == LAB_NOT_PROCESSED)) {
                        PQ.push(
                                PIX(srca(srcUpperLeft, o1), o1,
                                        insertionOrder++));
                        lab.set(LAB_QUEUED, labUpperLeft, o1);
                    }
                } // end for neighborhood
            } // end if label
        } // end x-loop
    } // end y-loop

    // Get the absolute minimum:
    minVal = valOfMin[0];
    for (int i = 1; i < valOfMin.size(); ++i) {
        if (minVal > valOfMin[i])
            minVal = valOfMin[i];
    }

    // Entering the priority queue ...
    while (!PQ.empty()) {
        PIX px = PQ.top();
        PQ.pop();
        VALUETYPE level = px.value;
        Diff2D o0 = px.offset;

        // normal flooding procedure
        int label1 = 0;
        int label2 = 0;

        // look to the neighborhood to determine the label of pixel o0.
        for (ITERATORTYPE iter = nbOffset.begin(); iter != nbOffset.end();
                ++iter) {
            Diff2D o1 = o0 + *iter;
            if ((!nbOffset.isOutsidePixel(o1))
                    && (maska(maskUpperLeft, o1) == maskLabel)) {
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
                        area[label1] = area[label1] + 1;
                    } else {
                        // in this case, a label has already been assigned to o0.
                        if (label1 != label2) {
                            // in this case, we have a meeting point of two lakes.
                            // we therefore have to fuse the two lakes.
                            if (minmax(valOfMin[label1], valOfMin[label2])
                                    || ((valOfMin[label1] == valOfMin[label2])
                                            && (area[label1] > area[label2]))) {
                                dynamics[label2] = level - valOfMin[label2];
                                area[label1] += area[label2];
                                equivalence[label2] = label1;
                            } else {
                                dynamics[label1] = level - valOfMin[label1];
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

    for (int i = 1; i < numberOfMinima + 1; ++i) {
        if (dynamics[i] == 0) {
            dynamics[i] = maxVal - minVal;
            break;
        }
    }

} // end of function

/////////////////
// DynMinima
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoDynMinima(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::pair<Iterator2, Accessor2> dest,
        std::vector<typename Accessor1::value_type> & dynamics,
        //std::vector<VECTYPE> & dynamics,
        NBTYPE & neighborOffset) {

    typedef typename Accessor1::value_type val_type;
    //std::vector<val_type> dynamics;

    morphoExtinctionValues(src.first, src.second, src.third, dest.first,
            dest.second, dynamics, neighborOffset,
            IsSmaller<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityBottomUp<val_type>());

}

//////////////////
// DynMinima
template<class Image1, class Image2, class NBTYPE>
void morphoDynMinima(Image1 & src, Image2 & dest, NBTYPE & neighborOffset) {
    //vigra::BasicImage::PixelType
    std::vector<typename Image1::value_type> dynamics;
    morphoDynMinima(srcImageRange(src), destImage(dest), dynamics,
            neighborOffset);
}

//////////////////////
// DynMinima with mask
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoMaskDynMinima(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::pair<Iterator2, Accessor2> mask,
        std::vector<typename Accessor1::value_type> & dynamics,
        std::vector<typename Accessor1::value_type> & values,
        NBTYPE & neighborOffset,
        typename Accessor2::value_type maskLabel = 255) {

    typedef typename Accessor1::value_type val_type;

    morphoMaskExtinctionValues(src.first, src.second, src.third, mask.first,
            mask.second, dynamics, values, neighborOffset,
            IsSmaller<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityBottomUp<val_type>(), maskLabel);
}

//////////////////////
// DynMinima with mask
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class NBTYPE>
void morphoMaskDynMinima(Iterator1 srcUpperLeft, Iterator1 srcLowerRight,
        Accessor1 srca, Iterator2 maskUpperLeft, Accessor2 maska,
        std::vector<typename Accessor1::value_type> & dynamics,
        std::vector<typename Accessor1::value_type> & values,
        NBTYPE & neighborOffset,
        typename Accessor2::value_type maskLabel = 255) {

    typedef typename Accessor1::value_type val_type;

    morphoMaskExtinctionValues(srcUpperLeft, srcLowerRight, srca, maskUpperLeft,
            maska, dynamics, values, neighborOffset,
            IsSmaller<typename Accessor1::value_type,
                    typename Accessor1::value_type>(),
            PriorityBottomUp<val_type>(), maskLabel);

}

}
;
}
;
#endif /*MORPHO_DYNAMIC_HXX_*/

