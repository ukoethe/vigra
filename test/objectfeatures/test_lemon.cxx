/************************************************************************/
/*                                                                      */
/*             Copyright 2011-2012 by Ullrich Koethe                    */
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

#include <iostream>
#include <sstream>
#include <map>
#include <set>

#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/accumulator.hxx>

namespace std {

template <unsigned int N, class T, class Stride>
ostream & operator<<(ostream & o, vigra::MultiArrayView<N, T, Stride> const & m)
{
    for(vigra::MultiArrayIndex k=0; k<m.size(); ++k)
        o << m[k] << " ";
    return o;
}

} // namespace std

using namespace vigra;

// mask cl.exe shortcomings
#if defined(_MSC_VER)
#pragma warning( disable : 4503 )
#endif

struct AccumulatorTest
{
    void testConvexHullFeatures2D()
    {
        using namespace vigra::acc;
        std::string prefix("testConvexHullFeatures2D(): ");

        int size = 6;
        MultiArray<2, int> mask(vigra::Shape2(size, size));

        mask(1, 1) = 1;
        mask(2, 1) = 1;
        mask(2, 2) = 1;
        mask(2, 3) = 1;
        mask(1, 3) = 1;
        mask(3, 1) = 1;
        mask(3, 3) = 1;
        mask(4, 1) = 1;
        mask(4, 3) = 1;

        AccumulatorChainArray<
                CoupledArrays<2, int>, 
                Select<LabelArg<1>, ConvexHullFeatures> > chf;
        chf.ignoreLabel(0);
        extractFeatures(mask, chf);

        getAccumulator<ConvexHullFeatures>(chf, 1).finalize();

        {
            TinyVector<double, 2> ref(2.5 - 1./18., 2.);
            shouldEqualSequenceTolerance(
                    get<ConvexHullFeatures>(chf, 1).inputCenter().begin(),
                    get<ConvexHullFeatures>(chf, 1).inputCenter().end(),
                    ref.begin(),
                    (std::numeric_limits<double>::epsilon() * 2));
        }
        {
            TinyVector<double, 2> ref(2.5, 2.);
            shouldEqualSequenceTolerance(
                    get<ConvexHullFeatures>(chf, 1).hullCenter().begin(),
                    get<ConvexHullFeatures>(chf, 1).hullCenter().end(),
                    ref.begin(),
                    (std::numeric_limits<double>::epsilon() * 2));
        }
        shouldEqual(get<ConvexHullFeatures>(chf, 1).inputVolume(), 9);
        shouldEqual(get<ConvexHullFeatures>(chf, 1).hullVolume(), 12);
        {
            TinyVector<double, 2> ref(8. / 3., 2.);
            shouldEqualSequenceTolerance(
                    get<ConvexHullFeatures>(chf, 1).defectCenter().begin(),
                    get<ConvexHullFeatures>(chf, 1).defectCenter().end(),
                    ref.begin(),
                    (std::numeric_limits<double>::epsilon() * 2));
        }
        shouldEqual(get<ConvexHullFeatures>(chf, 1).defectCount(), 2);
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeMean(),
                1.5,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeVariance(),
                0.5,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeSkewness(),
                0.0,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeKurtosis(),
                0.0,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectDisplacementMean(),
                ((2.5 - 1./18. - 1.) + 2*(3.5 - 2.5 + 1./18.))/3.,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).convexity(),
                9. / 12.,
                (std::numeric_limits<double>::epsilon() * 2));
    }

    void testConvexHullFeatures3D()
    {
        using namespace vigra::acc;
        std::string prefix("testConvexHullFeatures3D(): ");

        int size = 5;
        MultiArray<3, int> mask(vigra::Shape3(size, size, size), 0);
        for (int i = 0; i < 9; i++)
        {
            mask(i / 3 + 1, i % 3 + 1, 1) = 1;
            mask(i / 3 + 1, i % 3 + 1, 3) = 1;
        }
        mask(3, 1, 2) = 1; mask(3, 2, 2) = 1;
        mask(1, 3, 2) = 1; mask(2, 3, 2) = 1;

        // z = 0   | z = 1   | z = 2   | z = 3
        // --------+---------+---------+--------
        // 0 0 0 0 | 0 0 0 0 | 0 0 0 0 | 0 0 0 0
        // 0 0 0 0 | 0 x x x | 0 0 0 x | 0 x x x
        // 0 0 0 0 | 0 x x x | 0 0 0 x | 0 x x x
        // 0 0 0 0 | 0 x x x | 0 x x 0 | 0 x x x

        AccumulatorChainArray<
                CoupledArrays<3, int>,
                Select<LabelArg<1>, ConvexHullFeatures> > chf;
        chf.ignoreLabel(0);
        extractFeatures(mask, chf);

        getAccumulator<ConvexHullFeatures>(chf, 1).finalize();
        {
            // x and y coordinate: (7*1 + 7*2 + 8*3) / 22 = 45 / 22
            TinyVector<double, 3> ref(45. / 22., 45. / 22., 2.);
            shouldEqualSequenceTolerance(
                    get<ConvexHullFeatures>(chf, 1).inputCenter().begin(),
                    get<ConvexHullFeatures>(chf, 1).inputCenter().end(),
                    ref.begin(),
                    (std::numeric_limits<double>::epsilon() * 2));
        }
        {
            TinyVector<double, 3> ref(2., 2., 2.);
            shouldEqualSequenceTolerance(
                    get<ConvexHullFeatures>(chf, 1).hullCenter().begin(),
                    get<ConvexHullFeatures>(chf, 1).hullCenter().end(),
                    ref.begin(),
                    (std::numeric_limits<double>::epsilon() * 2));
        }
        shouldEqual(get<ConvexHullFeatures>(chf, 1).inputVolume(), 22);
        shouldEqual(get<ConvexHullFeatures>(chf, 1).hullVolume(), 27);
        {
            // x and y coordinate: (2*1 + 2*2 + 1*3) / 5 = 9 / 5
            TinyVector<double, 3> ref(9. / 5., 9. / 5., 2.);
            shouldEqualSequenceTolerance(
                    get<ConvexHullFeatures>(chf, 1).defectCenter().begin(),
                    get<ConvexHullFeatures>(chf, 1).defectCenter().end(),
                    ref.begin(),
                    (std::numeric_limits<double>::epsilon() * 2));
        }
        shouldEqual(get<ConvexHullFeatures>(chf, 1).defectCount(), 2);
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeMean(),
                2.5,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeVariance(),
                9. / 2.,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeSkewness(),
                0.0,
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectVolumeKurtosis(),
                0.0,
                (std::numeric_limits<double>::epsilon() * 2));
        // (sqrt(2) * 4 * (45 / 22 - 3 / 2) + sqrt(2) * 1 * (3 - 45 / 22)) / 5
        // = sqrt(2) / 5 * (90 / 11 - 3 - 45 / 22)
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).defectDisplacementMean(),
                sqrt(2.) / 5. * (90. / 11. - 3. - 45. / 22.),
                (std::numeric_limits<double>::epsilon() * 2));
        shouldEqualTolerance(
                get<ConvexHullFeatures>(chf, 1).convexity(),
                22. / 27.,
                (std::numeric_limits<double>::epsilon() * 2));
    }
};

struct FeaturesTestSuite : public vigra::test_suite
{
    FeaturesTestSuite()
        : vigra::test_suite("FeaturesTestSuite")
    {
        add(testCase(&AccumulatorTest::testConvexHullFeatures2D));
        add(testCase(&AccumulatorTest::testConvexHullFeatures3D));
    }
};

int main(int argc, char** argv)
{
    FeaturesTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
