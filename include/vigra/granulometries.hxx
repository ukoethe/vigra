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

#include "multi_array.hxx"
#include "accumulator.hxx"
#include "mathutil.hxx"
#include "morpho_basic.hxx"

namespace vigra {
namespace morpho {

using namespace vigra::acc;

void granulometricOpening(MultiArray<2, int> const &inputImage,
        std::vector<StructuringElement2D> const &SEs,
        std::vector<double> &areas) {

    if (SEs.size() < 1) {
        return;
    }

    areas.clear();

    MultiArray<2, int> openImage(inputImage.width(), inputImage.height());

    int inputArea = countNonZero(inputImage);

    std::vector<StructuringElement2D>::const_iterator it;
    for (it = SEs.begin(); it != SEs.end(); ++it) {
        morphoOpening(inputImage, openImage, *it);
        int openArea = countNonZero(openImage);
        double areaRatio = (double) (inputArea - openArea) / (double) inputArea;
        areas.push_back(areaRatio);

    }

}

void createSEPyramid(int radius, std::vector<StructuringElement2D> &SEs,
        int count) {

    if (count < 1) {
        return;
    }

    SEs.clear();

    int currentRadius = radius;

    for (int i = 0; i < count; ++i) {
        StructuringElement2D SE;
        generateDiscSE(currentRadius, SE);
        SEs.push_back(SE);
        currentRadius += radius;
    }

}

void createSEPyramid(StructuringElement2D const &SE,
        std::vector<StructuringElement2D> &SEs, int count) {

    if (count < 1) {
        return;
    }

    SEs.clear();

    SEs.push_back(SE);

    StructuringElement2D dilatedSE = SE;

    for (int i = 1; i < count; ++i) {
        StructuringElement2D tmp;
        dilateSE(dilatedSE, SE, tmp);
        SEs.push_back(tmp);
        dilatedSE = tmp;
    }

}

int getSEArea(StructuringElement2D SE) {
    return std::distance(SE.begin(), SE.end());
}

////////////////////////////////////////////////////////////////////////////////////////////
//
//  Granulometries:
//    Provide either a radius or an structuring element to perform the calculation.
//    The output will be a granulometric curve expressed as
//    the normalized areas of the elements covered by the structuring element.
//    And the size of the structuring element at each point of the curve.
//
////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Disk shaped granulometries.
 */
void granulometries(MultiArray<2, int> const &inputImage, int radius, int count,
        std::vector<double> &granulometricCurve,
        std::vector<int> &axisAreas) {

    std::vector<StructuringElement2D> SEs;
    std::vector<double> areaRatios;

    createSEPyramid(radius, SEs, count);

    granulometricOpening(inputImage, SEs, areaRatios);

    granulometricCurve.push_back(areaRatios[0]);

    for (std::vector<StructuringElement2D>::iterator it = SEs.begin();
            it != SEs.end(); ++it) {
        axisAreas.push_back(getSEArea(*it));
    }

    for (unsigned i = 1; i < areaRatios.size(); ++i) {
        granulometricCurve.push_back(areaRatios[i] - areaRatios[i - 1]);
    }

}

/*
 * General granulometries.
 * The provided structuring element will be dilated to get bigger versions.
 */
void granulometries(MultiArray<2, int> const &inputImage, StructuringElement2D const &SE,
        int count, std::vector<double> &granulometricCurve,
        std::vector<int> &axisAreas) {
    std::vector<StructuringElement2D> SEs;
    std::vector<double> areaRatios;
    createSEPyramid(SE, SEs, count);
    granulometricOpening(inputImage, SEs, areaRatios);
    granulometricCurve.push_back(areaRatios[0]);

    for (std::vector<StructuringElement2D>::iterator it = SEs.begin();
            it != SEs.end(); ++it) {
        axisAreas.push_back(getSEArea(*it));
    }

    for (unsigned i = 1; i < areaRatios.size(); ++i) {
        granulometricCurve.push_back(areaRatios[i] - areaRatios[i - 1]);
    }

}


} // namespace morpho

} // namespace vigra

