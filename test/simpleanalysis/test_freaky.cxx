#include <sstream>
#include <iostream>
#include <functional>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/distancetransform.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/cornerdetection.hxx"
#include "vigra/symmetry.hxx"

using namespace vigra;

struct RegionGrowingTest
{
    typedef vigra::DImage Image;

    RegionGrowingTest()
    : img(7,7), seeds(7,7)
    {
        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        Image tmp(7,7);

        Image::ScanOrderIterator i = tmp.begin();
        Image::ScanOrderIterator end = tmp.end();
        Image::Accessor acc = tmp.accessor();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }

        distanceTransform(srcImageRange(tmp), destImage(img), 0.0, 2);

        seeds = 0;
        labelImageWithBackground(srcImageRange(tmp), destImage(seeds),
                                 false, 0.0);
    }

    struct DirectCostFunctor
    {
        typedef double argument_type;
        typedef double result_type;
        typedef double cost_type;

        void operator()(double const &) {}

        double const & cost(double const & v) const
        {
            return v;
        }
    };

    void voronoiTest()
    {
        Image res(img);

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing(srcImageRange(img), srcImage(seeds),
                            destImage(res), cost);

        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;

        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, vigra::Diff2D(x,y));
                double dist1 = VIGRA_CSTD::sqrt((2.0 - x)*(2.0 - x) +
                                                (2.0 - y)*(2.0 - y));
                double dist2 = VIGRA_CSTD::sqrt((5.0 - x)*(5.0 - x) +
                                                (5.0 - y)*(5.0 - y));

                if(!(dist1 == dist2))
                {
                    std::stringstream ss;
                    ss << x << "/" << y
                       << ": (dist1=" << dist1
                       << " == dist2=" << dist2 << ") == "
                       << ((dist1 == dist2) ? "tr" : "fal")
                       << (!(dist1 != dist2) ? "ue" : "se")
                       << "\n";

                    shouldMsg(dist == ((dist1 < dist2) ? 1.0 : 2.0),
                              ss.str().c_str());
                }
            }
        }
    }

    void voronoiWithBorderTest()
    {
        Image res(img);
        Image::value_type reference[] = {
            1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 0,
            1, 1, 1, 1, 1, 0, 2,
            1, 1, 1, 1, 0, 2, 2,
            1, 1, 1, 0, 2, 2, 2,
            1, 1, 0, 2, 2, 2, 2,
            1, 0, 2, 2, 2, 2, 2
        };

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing(srcImageRange(img), srcImage(seeds),
                            destImage(res), cost, KeepContours);

        shouldEqualSequence(res.begin(), res.end(), reference);
    }

    Image img, seeds;
};

struct SimpleAnalysisTestSuite
: public vigra::test_suite
{
    SimpleAnalysisTestSuite()
    : vigra::test_suite("SimpleAnalysisTestSuite")
    {
        add( testCase( &RegionGrowingTest::voronoiTest));
    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
