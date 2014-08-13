
#include "multi_array.hxx"
#include "accumulator.hxx"
#include "mathutil.hxx"
#include "morpho_basic.hxx"

namespace vigra {
namespace morpho {

using namespace vigra::acc;


template<typename T>
int countNonZero(MultiArray<2, T> const &array) {

    vigra::MultiArray<2, double> data(array.width(), array.height());
    AccumulatorChainArray<CoupledArrays<2, double, int>,
            Select<LabelArg<2>, Count> > a;
    extractFeatures(data, array, a);

    int nonZero = array.width() * array.height() - get<Count>(a, 0);
    return nonZero;
}

void granulometricOpening(MultiArray<2, int> const &inputImage,
        std::vector<structuringElement2D> const &SEs,
        std::vector<double> &areas) {

    if (SEs.size() < 1) {
        return;
    }

    areas.clear();

    MultiArray<2, int> openImage(inputImage.width(), inputImage.height());

    int inputArea = countNonZero(inputImage);

    std::vector<structuringElement2D>::const_iterator it;
    for (it = SEs.begin(); it != SEs.end(); ++it) {
        morphoOpening(inputImage, openImage, *it);
        int openArea = countNonZero(openImage);
        double areaRatio = (double) (inputArea - openArea) / (double) inputArea;
        areas.push_back(areaRatio);

    }

}

/*
 * Generates image from SE
 */
void renderSE(structuringElement2D SE, MultiArray<2, int> &outputImage,
        Diff2D center) {

    std::vector<Diff2D>::const_iterator it;

    for (int y = 0; y < outputImage.height(); ++y) {
        for (int x = 0; x < outputImage.width(); ++x) {
            int posX = x - center.x;
            int posY = y - center.y;
            it = find(SE.begin(), SE.end(), Diff2D(posX, posY));
            if (it != SE.end()) {
                outputImage(x, y) = 1;
            }

        }
    }
}

/*
 * Generates SE from image
 */
void generateSE(MultiArray<2, int> seImage, structuringElement2D &se,
        Diff2D center) {

    std::vector<Diff2D> strElCoordinates;

    for (int y = 0; y < seImage.height(); ++y) {
        for (int x = 0; x < seImage.width(); ++x) {
            if (seImage(x, y) != 0) {
                int posX = x - center.x;
                int posY = y - center.y;
                strElCoordinates.push_back(Diff2D(posX, posY));
            }

        }
    }

    se = structuringElement2D(strElCoordinates);
}

/*
 * Dilates inputSE using SE as structuring element
 */
void dilateSE(structuringElement2D inputSE, structuringElement2D SE,
        structuringElement2D &outputSE) {

    int extensionY;
    int extensionX;

    extensionY = inputSE.maxOffset().y - inputSE.minOffset().y + 2;
    extensionX = inputSE.maxOffset().x - inputSE.minOffset().x + 2;

    MultiArray<2, int> SERender(extensionX * 2, extensionY * 2);
    MultiArray<2, int> SERenderDilated(extensionX * 2, extensionY * 2);
    Diff2D center(extensionX - 1, extensionY - 1);
    renderSE(inputSE, SERender, center);

    morphoDilation(SERender, SERenderDilated, SE);

    generateSE(SERenderDilated, outputSE, center);

}

void createSEPyramid(int radius, std::vector<structuringElement2D> &SEs,
        int count) {

    if (count < 1) {
        return;
    }

    SEs.clear();

    int currentRadius = radius;

    for (int i = 0; i < count; ++i) {
        structuringElement2D SE;
        generateDiscSE(currentRadius, SE);
        SEs.push_back(SE);
        currentRadius += radius;
    }

}

void createSEPyramid(structuringElement2D const &SE,
        std::vector<structuringElement2D> &SEs, int count) {

    if (count < 1) {
        return;
    }

    SEs.clear();

    SEs.push_back(SE);

    structuringElement2D dilatedSE = SE;

    for (int i = 1; i < count; ++i) {
        structuringElement2D tmp;
        dilateSE(dilatedSE, SE, tmp);
        SEs.push_back(tmp);
        dilatedSE = tmp;
    }

}

int getSEArea(structuringElement2D SE) {
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
        std::vector<double> &axisAreas) {

    std::vector<structuringElement2D> SEs;
    std::vector<double> areaRatios;

    createSEPyramid(radius, SEs, count);

    granulometricOpening(inputImage, SEs, areaRatios);

    granulometricCurve.push_back(areaRatios[0]);

    for (std::vector<structuringElement2D>::iterator it = SEs.begin();
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
void granulometries(MultiArray<2, int> const &inputImage, structuringElement2D const &SE,
        int count, std::vector<double> &granulometricCurve,
        std::vector<double> &axisAreas) {
    std::vector<structuringElement2D> SEs;
    std::vector<double> areaRatios;
    createSEPyramid(SE, SEs, count);
    granulometricOpening(inputImage, SEs, areaRatios);
    granulometricCurve.push_back(areaRatios[0]);

    for (std::vector<structuringElement2D>::iterator it = SEs.begin();
            it != SEs.end(); ++it) {
        axisAreas.push_back(getSEArea(*it));
    }

    for (unsigned i = 1; i < areaRatios.size(); ++i) {
        granulometricCurve.push_back(areaRatios[i] - areaRatios[i - 1]);
    }

}


} // namespace morpho

} // namespace vigra

