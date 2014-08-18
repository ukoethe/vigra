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

#ifndef MORPHO_BASIC_HXX_
#define MORPHO_BASIC_HXX_

#include "morpho_utilities.hxx"
#include "copyimage.hxx"
#include "combineimages.hxx"
#include "pointOrdering.hxx"

namespace vigra {
namespace morpho {

/*
 * Structuring element erosion.
 */
template<class SElement>
void erodeSE(SElement inputSE, SElement SE, SElement &outputSE) {

    using namespace std;
    using namespace vigra::detail;

    std::vector<Diff2D> outputSECoordinates;
    int xMin = inputSE.minOffset().x;
    int xMax = inputSE.maxOffset().x;
    int yMin = inputSE.minOffset().y;
    int yMax = inputSE.maxOffset().y;

    std::vector<Diff2D> inputSECoordinates(
            distance(inputSE.begin(), inputSE.end()));
    std::copy(inputSE.begin(), inputSE.end(), inputSECoordinates.begin());
    sort(inputSECoordinates.begin(), inputSECoordinates.end(),
            PointYXOrdering<ASCENDING>());
    inputSECoordinates.erase(
            unique(inputSECoordinates.begin(), inputSECoordinates.end()),
            inputSECoordinates.end());

    std::vector<Diff2D> SECoordinates(distance(SE.begin(), SE.end()));
    std::copy(SE.begin(), SE.end(), SECoordinates.begin());
    sort(SECoordinates.begin(), SECoordinates.end(),
            PointYXOrdering<ASCENDING>());
    SECoordinates.erase(unique(SECoordinates.begin(), SECoordinates.end()),
            SECoordinates.end());

    for (int y = yMin; y <= yMax; ++y) {
        for (int x = xMin; x <= xMax; ++x) {
            std::vector<Diff2D>::const_iterator itSE = SECoordinates.begin();
            std::vector<Diff2D>::const_iterator itInputSE;
            Diff2D neighbour(x + itSE->x, y + itSE->y);
            itInputSE = find(inputSECoordinates.begin(),
                    inputSECoordinates.end(), neighbour);

            while (itSE != SECoordinates.end()
                    && itInputSE != inputSECoordinates.end()) {
                neighbour = Diff2D(x + itSE->x, y + itSE->y);
                while (neighbour != *itInputSE
                        && itInputSE != inputSECoordinates.end()) {
                    itInputSE++;
                }
                if (neighbour == *itInputSE) {
                    itSE++;
                }
            }

            if (itSE == SECoordinates.end()) {
                outputSECoordinates.push_back(Diff2D(x, y));
            }
        }
    }

    outputSE = StructuringElement2D(outputSECoordinates);
}

/*
 * Structuring element dilation.
 */
template<class SElement>
void dilateSE(SElement inputSE, SElement SE, SElement &outputSE) {

    using namespace std;
    using namespace vigra::detail;

    std::vector<Diff2D> outputSECoordinates;
    int xMin = inputSE.minOffset().x + SE.minOffset().x;
    int xMax = inputSE.maxOffset().x + SE.maxOffset().x;
    int yMin = inputSE.minOffset().y + SE.minOffset().y;
    int yMax = inputSE.maxOffset().y + SE.maxOffset().y;

    std::vector<Diff2D> inputSECoordinates(
            distance(inputSE.begin(), inputSE.end()));
    std::copy(inputSE.begin(), inputSE.end(), inputSECoordinates.begin());
    sort(inputSECoordinates.begin(), inputSECoordinates.end(),
            PointYXOrdering<ASCENDING>());
    inputSECoordinates.erase(
            unique(inputSECoordinates.begin(), inputSECoordinates.end()),
            inputSECoordinates.end());

    std::vector<Diff2D> SECoordinates(distance(SE.begin(), SE.end()));
    std::copy(SE.begin(), SE.end(), SECoordinates.begin());
    sort(SECoordinates.begin(), SECoordinates.end(),
            PointYXOrdering<ASCENDING>());
    SECoordinates.erase(unique(SECoordinates.begin(), SECoordinates.end()),
            SECoordinates.end());

    for (int y = yMin; y <= yMax; ++y) {
        for (int x = xMin; x <= xMax; ++x) {
            std::vector<Diff2D>::const_iterator itSE = SECoordinates.begin();

            for (itSE = SECoordinates.begin(); itSE != SECoordinates.end();
                    ++itSE) {
                Diff2D neighbour(x + itSE->x, y + itSE->y);
                std::vector<Diff2D>::const_iterator itInputSE = find(
                        inputSECoordinates.begin(), inputSECoordinates.end(),
                        neighbour);

                if (std::binary_search(inputSECoordinates.begin(),
                        inputSECoordinates.end(), neighbour,
                        PointYXOrdering<ASCENDING>())) {
                    break;
                }
            }

            if (itSE != SECoordinates.end()) {
                outputSECoordinates.push_back(Diff2D(x, y));
            }
        }
    }

    outputSE = StructuringElement2D(outputSECoordinates);
}

/*
 *  Structuring element operation.
 *  By providing a functor you can achieve erosion, dilation,
 *  rank filtering or any other morphological operation you
 *  can express with this kind of functors.
 */
/*
 * This code is less efficient than discRankOrderFilter.
 * That may be because discRankOrderFilter doesn't iterate over the elements
 * of an structuring element.
 * If that's true, the efficiency of this function could be improved processing
 * the structuring element in chunks.
 */
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement, class Functor>
void morphoBasicSEOperation(Iterator1 srcUpperLeft, Iterator1 srcLowerRight,
        Accessor1 srca, Iterator2 destUpperLeft, Accessor2 desta, SElement & se,
        Functor f) {

    // border treatment
    // offsets correspond to the maximal extension of the SE.
    Diff2D minOffset = se.minOffset();
    Diff2D maxOffset = se.maxOffset();

    const Iterator1 upperLeftCorner = srcUpperLeft;
    const Iterator1 lowerRightCorner = srcLowerRight;
    const Iterator1 upperLeftCritical = srcUpperLeft - minOffset;
    const Iterator1 lowerRightCritical = srcLowerRight - maxOffset;

    /*
     * If the structuring element exceeds the image,
     * always perform boundary checking.
     * There will always be part of the structuring element outside the image.
     */
    bool seExceedsImage = (maxOffset.x - minOffset.x)
            > (srcLowerRight.x - srcUpperLeft.x)
            || (maxOffset.y - minOffset.y) > (srcLowerRight.y - srcUpperLeft.y);

    if (seExceedsImage) {

        for (; srcUpperLeft.y < srcLowerRight.y;
                ++srcUpperLeft.y, ++destUpperLeft.y) {
            Iterator1 scurrent = srcUpperLeft;
            Iterator2 dcurrent = destUpperLeft;

            for (; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x) {
                f.reset();
                for (typename SElement::ITERATORTYPE iter = se.begin();
                        iter != se.end(); ++iter) {
                    Iterator1 neighbour = scurrent + *iter;
                    if ((neighbour.y >= upperLeftCorner.y)
                            && (neighbour.y < lowerRightCorner.y)
                            && (neighbour.x >= upperLeftCorner.x)
                            && (neighbour.x < lowerRightCorner.x))
                        f(srca(neighbour));
                }
                desta.set(f(), dcurrent);
            } // end of x loop
        } // end for the first y-loop.

    } else {
        /*
         * If the structuring element is smaller than the image,
         * just perform boundary checking when needed.
         * This method will minimize the boundary checking,
         */

        // Upper section
        for (; srcUpperLeft.y < upperLeftCritical.y;
                ++srcUpperLeft.y, ++destUpperLeft.y) {
            Iterator1 scurrent = srcUpperLeft;
            Iterator2 dcurrent = destUpperLeft;

            for (; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x) {
                f.reset();
                for (typename SElement::CONST_ITERATORTYPE iter = se.begin();
                        iter != se.end(); ++iter) {
                    Iterator1 neighbour = scurrent + *iter;
                    if ((neighbour.y >= upperLeftCorner.y)
                            && (neighbour.x >= upperLeftCorner.x)
                            && (neighbour.x < lowerRightCorner.x))
                        f(srca(neighbour));
                }
                desta.set(f(), dcurrent);
            } // end of x loop
        } // end for the first y-loop.

        for (; srcUpperLeft.y < lowerRightCritical.y;
                ++srcUpperLeft.y, ++destUpperLeft.y) {
            Iterator1 scurrent = srcUpperLeft;
            Iterator2 dcurrent = destUpperLeft;

            // x-loop: the left side
            for (; scurrent.x < upperLeftCritical.x;
                    ++scurrent.x, ++dcurrent.x) {
                f.reset();
                for (typename SElement::CONST_ITERATORTYPE iter = se.begin();
                        iter != se.end(); ++iter) {
                    Iterator1 neighbour = scurrent + *iter;
                    if (neighbour.x >= upperLeftCorner.x)
                        f(srca(neighbour));
                }
                desta.set(f(), dcurrent);
            } // end of x loop (left)

            for (; scurrent.x < lowerRightCritical.x;
                    ++scurrent.x, ++dcurrent.x) {
                f.reset();
                for (typename SElement::CONST_ITERATORTYPE iter = se.begin();
                        iter != se.end(); ++iter) {
                    Iterator1 neighbour = scurrent + *iter;
                    f(srca(neighbour));
                }
                desta.set(f(), dcurrent);
            } // end of the middle x loop

            // the right side
            for (; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x) {
                f.reset();
                for (typename SElement::CONST_ITERATORTYPE iter = se.begin();
                        iter != se.end(); ++iter) {
                    Iterator1 neighbour = scurrent + *iter;
                    if (neighbour.x < lowerRightCorner.x)
                        f(srca(neighbour));
                }
                desta.set(f(), dcurrent);
            } // end of the right x loop
        } // end of y loop (middle)

        // Lower section
        for (; srcUpperLeft.y < srcLowerRight.y;
                ++srcUpperLeft.y, ++destUpperLeft.y) {
            Iterator1 scurrent = srcUpperLeft;
            Iterator2 dcurrent = destUpperLeft;

            for (; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x) {
                f.reset();
                for (typename SElement::CONST_ITERATORTYPE iter = se.begin();
                        iter != se.end(); ++iter) {
                    Iterator1 neighbour = scurrent + *iter;
                    if ((neighbour.y < lowerRightCorner.y)
                            && (neighbour.x < lowerRightCorner.x)
                            && (neighbour.x >= upperLeftCorner.x))
                        f(srca(neighbour));
                }
                desta.set(f(), dcurrent);
            } // end of x loop
        } // end for the lower y-loop.

    }
} // end of morphoBasicSEOperation

/*
 * Perform the morphological operation as many times as
 * the structuring element needs reading sizeMultiplier
 */
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement, class Functor>
void runPasses(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se,
        Functor f) {

    vigra_precondition(se.sizeMultiplier > 0,
            "Structuring element size multiplier must be higher than 0.");

    morphoBasicSEOperation(src.first, src.second, src.third, dest.first,
            dest.third, se, f);

    if (se.sizeMultiplier > 1) {

        vigra::BasicImage<typename Accessor2::value_type> temp(
                dest.second - dest.first);

        // a morphological rank filtering with se of size n
        // corresponds to n morphological filterings with size 1.
        for (int i = 1; i < se.sizeMultiplier; i++) {
            if (i % 2 == 0)
                morphoBasicSEOperation(temp.upperLeft(), temp.lowerRight(),
                        temp.accessor(), dest.first, dest.third, se, f);
            else
                morphoBasicSEOperation(dest.first, dest.second, dest.third,
                        temp.upperLeft(), temp.accessor(), se, f);
        }

        if (se.sizeMultiplier % 2 == 0)
            vigra::copyImage(temp.upperLeft(), temp.lowerRight(),
                    temp.accessor(), dest.first, dest.third);
    }

}

/////////////////////////////////////////////////////////////////////////
// RANK FILTER
/////////////////////////////////////////////////////////////////////////
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoRankFilter(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se,
        double rank) {

    runPasses(src, dest, se, Rank<typename Accessor1::value_type>(rank));

}

/////////////////////////////////////////////////////////////////////////
// EROSION AND DILATION
/////////////////////////////////////////////////////////////////////////

// Morphological dilation
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoDilation(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se) {

    runPasses(src, dest, se, Max<typename Accessor1::value_type>());
} // end of dilation

// Morphological erosion
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoErosion(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se) {

    runPasses(src, dest, se, Min<typename Accessor1::value_type>());

} // end of erosion

/////////////////////////////////////////////////////////////////////////
// OPENING AND CLOSING
/////////////////////////////////////////////////////////////////////////

// Morphological opening
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoOpening(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se) {
    vigra::BasicImage<typename Accessor1::value_type> temp(
            src.second - src.first);
    morphoErosion(src, vigra::destImageRange(temp), se);
    se.transpose();
    morphoDilation(vigra::srcImageRange(temp), dest, se);
    se.transpose();
} // end of opening

// Morphological closing
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoClosing(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se) {
    vigra::BasicImage<typename Accessor1::value_type> temp(
            src.second - src.first);
    morphoDilation(src, vigra::destImageRange(temp), se);
    se.transpose();
    morphoErosion(vigra::srcImageRange(temp), dest, se);
    se.transpose();
} // end of closing

// Morphological gradient
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoGradient(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se,
        typename Accessor2::value_type markerVal = 255) {
    typedef typename Accessor1::value_type INTYPE;
    typedef typename Accessor2::value_type OUTTYPE;

    vigra::BasicImage<INTYPE> dil(src.second - src.first);
    vigra::BasicImage<INTYPE> ero(src.second - src.first);
    morphoDilation(src, vigra::destImageRange(dil), se);
    morphoErosion(src, vigra::destImageRange(ero), se);
    vigra::combineTwoImages(srcImageRange(dil), srcImage(ero),
            destIter(dest.first, dest.third), std::minus<INTYPE>());

} // end of morphogradient

// External gradient
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoExternalGradient(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se,
        typename Accessor2::value_type markerVal = 255) {
    typedef typename Accessor1::value_type INTYPE;
    typedef typename Accessor2::value_type OUTTYPE;

    vigra::BasicImage<INTYPE> dil(src.second - src.first);
    morphoDilation(src, vigra::destImageRange(dil), se);
    vigra::combineTwoImages(srcImageRange(dil), srcIter(src.first, src.third),
            destIter(dest.first, dest.third), std::minus<INTYPE>());

}

// Internal gradient
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,
        class SElement>
void morphoInternalGradient(vigra::triple<Iterator1, Iterator1, Accessor1> src,
        vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se,
        typename Accessor2::value_type markerVal = 255) {
    typedef typename Accessor1::value_type INTYPE;
    typedef typename Accessor2::value_type OUTTYPE;

    vigra::BasicImage<INTYPE> ero(src.second - src.first);
    morphoErosion(src, vigra::destImageRange(ero), se);
    vigra::combineTwoImages(srcIterRange(src.first, src.second, src.third),
            srcImage(ero), destIter(dest.first, dest.third),
            std::minus<INTYPE>());

}

template<class Image1, class Image2, class SElement>
void morphoInternalGradient(const Image1 & imin, Image2 & imout, SElement se) {
    morphoInternalGradient(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoExternalGradient(const Image1 & imin, Image2 & imout, SElement se) {
    morphoExternalGradient(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoGradient(const Image1 & imin, Image2 & imout, SElement se) {
    morphoGradient(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoRankFilter(const Image1 & imin, Image2 & imout, SElement se,
        double rank) {
    morphoRankFilter(srcImageRange(imin), destImageRange(imout), se, rank);
}

template<class Image1, class Image2, class SElement>
void morphoErosion(const Image1 & imin, Image2 & imout, SElement se) {
    morphoErosion(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoDilation(const Image1 & imin, Image2 & imout, SElement se) {
    morphoDilation(srcImageRange(imin), destImageRange(imout), se);
}

// Open and close
template<class Image1, class Image2, class SElement>
void morphoOpening(const Image1 & imin, Image2 & imout, SElement se) {
    morphoOpening(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoClosing(const Image1 & imin, Image2 & imout, SElement se) {
    morphoClosing(srcImageRange(imin), destImageRange(imout), se);
}

}
;
/* end of namespace morpho */
}
;
/* end of namespace vigra */

#endif /*BASIC_MORPHO_HXX_*/

