/************************************************************************/
/*                                                                      */
/*               Copyright 2011-2012 by Ullrich Koethe                  */
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

#include <cstdio>
#include <iostream>
#include <vector>
#include <unittest.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/k-means.hxx>


//added includes
#include "vigra/matrix.hxx"

using namespace vigra;

class UnsupervisedClusteringTest
{

private:

    /* \brief Run a k-means clustering unit test.

       \a clusterCenters must have been initialized with the starting centers;
        after this function returns it contains the final clusterCenters.

        \a expectedClusterCenters must contain the centroids expected at
        the completion of the k-means algorithm. If \a
        expectedClusterCenters does not match \a clusterCenters (within
        tolerance \a epsilon), the test fails.
    */

    template<class T1, class Stride1, class T2, class Stride2>
    void
    kMeansRunTest(MultiArrayView<2, T1, Stride1> const & data,
                  MultiArrayView<2, T2, Stride2> & clusterCenters,
                  MultiArrayView<2, T2, Stride2> const & expectedClusterCenters,
                  double epsilon = 1.0e-10,
                  KMeansOptions const & options = KMeansOptions())
    {
        int clusterCount = columnCount(clusterCenters);
        int expClusterCount = columnCount(expectedClusterCenters);
        vigra_precondition(clusterCount == expClusterCount,
        "kMeansRunTest(): number of initial and expected clusters must be equal");

        int numSamples = columnCount(data);
        MultiArray<1, unsigned int> assignments(numSamples);
        kMeans(data, assignments, clusterCenters, options);

        for (int i=0; i<clusterCount; i++) 
        {
            shouldEqualSequenceTolerance(columnVector(expectedClusterCenters,i).begin(),
            columnVector(expectedClusterCenters,i).end(),
            columnVector(clusterCenters, i).begin(),
            epsilon);
        }
    }


public:

    UnsupervisedClusteringTest()
    {}


    /*
        A simple sanity check: 2 samples in 3 dimensions, with 2
        clusters initialized one to each sample. Performing the
        clustering should not change the centroids.
    */

    void testKMeansSanity()
    {
        int numFeatures = 3;
        int numSamples = 2, clusterCount = 2;

        static const double dataVals[] =
        {
            0., 0., 0.,
            1., 1., 1.,
        };
        MultiArray<2, double> data(Shape2(numFeatures, numSamples), dataVals);
        MultiArray<2, double> clusterCenters(Shape2(numFeatures, clusterCount), dataVals);
        MultiArray<2, double> expectedClusterCenters(Shape2(numFeatures, clusterCount), dataVals);
        kMeansRunTest(data, clusterCenters, expectedClusterCenters);
    }

    /*
        A slightly more complicated example. 24 samples in 3
        dimensions. Each sample is close to one of the eight vertices of
        a unit cube.

        Each of the eight clusters initializes outside one of the vertices;
        clustering should quickly converge.
    */

    void testKMeansCube()
    {
        int numFeatures = 3;
        int numSamples = 24, clusterCount = 8;

        static const double dataVals[] =
        {
            0.5, 0.5, 0.5,
            0.51, 0.51, 0.51,
            0.49, 0.49, 0.49,
            0.5, 0.5, -0.5,
            0.51, 0.51, -0.51,
            0.49, 0.49, -0.49,
            0.5, -0.5, 0.5,
            0.51, -0.51, 0.51,
            0.49, -0.49, 0.49,
            0.5, -0.5, -0.5,
            0.51, -0.51, -0.51,
            0.49, -0.49, -0.49,
            -0.5, 0.5, 0.5,
            -0.51, 0.51, 0.51,
            -0.49, 0.49, 0.49,
            -0.5, 0.5, -0.5,
            -0.51, 0.51, -0.51,
            -0.49, 0.49, -0.49,
            -0.5, -0.5, 0.5,
            -0.51, -0.51, 0.51,
            -0.49, -0.49, 0.49,
            -0.5, -0.5, -0.5,
            -0.51, -0.51, -0.51,
            -0.49, -0.49, -0.49,
        };
        MultiArray<2, double> data(Shape2(numFeatures, numSamples), dataVals);

        static const double clusterVals[] =
        {
            0.6, 0.6, 0.6,
            0.6, 0.6, -0.6,
            0.6, -0.6, 0.6,
            0.6, -0.6, -0.6,
            -0.6, 0.6, 0.6,
            -0.6, 0.6, -0.6,
            -0.6, -0.6, 0.6,
            -0.6, -0.6, -0.6,
        };
        MultiArray<2, double> clusterCenters(Shape2(numFeatures, clusterCount), clusterVals);

        static const double expClusterVals[] =
        {
            0.5, 0.5, 0.5,
            0.5, 0.5, -0.5,
            0.5, -0.5, 0.5,
            0.5, -0.5, -0.5,
            -0.5, 0.5, 0.5,
            -0.5, 0.5, -0.5,
            -0.5, -0.5, 0.5,
            -0.5, -0.5, -0.5,
        };
        MultiArray<2, double> expectedClusterCenters(Shape2(numFeatures, clusterCount), expClusterVals);

        kMeansRunTest(data, clusterCenters, expectedClusterCenters);
        kMeansRunTest(data, clusterCenters, expectedClusterCenters, 1e-2,
                      KMeansOptions().maxIterations(1000).sampleSize(16));
    }


    /*
        Test for correct behaviour when there are more clusters than data
        points.

        The extra clusters should go to zero, but this is a design
        detail. The real point is to ensure that the algorithm
        implementation does not naively divide by 0 when a cluster is
        empty.
    */

    void testKMeansExtraClusters()
    {
        int numFeatures = 3;
        int numSamples = 2, clusterCount = 4;

        static const double dataVals[] =
        {
            0., 0., 0.,
            1., 1., 1.,
        };
        MultiArray<2, double> data(Shape2(numFeatures, numSamples), dataVals);

        static const double clusterVals[] =
        {
            0., 0., 0.,
            1., 1., 1.,
            1, 2., 1.,
            2., 4., 1.,
        };
        MultiArray<2, double> clusterCenters(Shape2(numFeatures, clusterCount), clusterVals);

        static const double expClusterVals[] =
        {
            0., 0., 0.,
            1., 1., 1.,
            0., 0., 0.,
            0., 0., 0.
        };

        MultiArray<2, double> expectedClusterCenters(Shape2(numFeatures, clusterCount), expClusterVals);
        kMeansRunTest(data, clusterCenters, expectedClusterCenters);
    }


    /*
        Check that the correct errors are thrown.
    */

    void testKMeansErrors()
    {
        int numFeatures = 3;
        int numSamples = 2, clusterCount = 2;

        static const double dataVals[] =
        {
            0., 0., 0.,
            1., 1., 1.,
        };
        MultiArray<2, double> data(Shape2(numFeatures, numSamples), dataVals);
        MultiArray<2, double> clusterCenters(Shape2(numFeatures, clusterCount), dataVals);
        MultiArray<2, double> expectedClusterCenters(Shape2(numFeatures, clusterCount), dataVals);
        MultiArray<1, unsigned int> assignments(numSamples);

        //not enough clusters
        MultiArray<2, double> clusterCenters2(Shape2(numFeatures, clusterCount - 1), dataVals);
        try 
        {
            kMeans(data, assignments, clusterCenters2);
            failTest("kMeans() failed to throw exception.");
        } 
        catch(PreconditionViolation & e)
        {
            std::string expected("\nPrecondition violation!\nkMeans(): At least 2 clusters are required."),
            actual(e.what());
            shouldEqual(actual.substr(0, expected.size()), expected);
        }
        
        //number of features mismatch
        MultiArray<2, double> data2(Shape2(numFeatures - 1, numSamples), dataVals);
        try 
        {
            kMeans(data2, assignments, clusterCenters);
            failTest("kMeans() failed to throw exception.");
        } 
        catch(PreconditionViolation & e)
        {
            std::string expected("\nPrecondition violation!\nkMeans(): Number of features mismatch between data and clusterCenter arrays."),
            actual(e.what());
            shouldEqual(actual.substr(0, expected.size()), expected);
        }

        //number of samples mismatch
        MultiArray<2, double> data3(Shape2(numFeatures, numSamples - 1), dataVals);
        try 
        {
            kMeans(data3, assignments, clusterCenters);
            failTest("kMeans() failed to throw exception.");
        } 
        catch(PreconditionViolation & e)
        {
            std::string expected("\nPrecondition violation!\nkMeans(): Number of samples mismatch between data and assignments arrays."),
            actual(e.what());
            shouldEqual(actual.substr(0, expected.size()), expected);
        }
    }
};

struct UnsupervisedLearningTestSuite : public vigra::test_suite
{
    UnsupervisedLearningTestSuite()
        : vigra::test_suite("UnsupervisedLearningTestSuite")
    {
        add(testCase(&UnsupervisedClusteringTest::testKMeansSanity));
        add(testCase(&UnsupervisedClusteringTest::testKMeansCube));
        add(testCase(&UnsupervisedClusteringTest::testKMeansExtraClusters));
        add(testCase(&UnsupervisedClusteringTest::testKMeansErrors));
    }
};


int main (int argc, char ** argv)
{
    UnsupervisedLearningTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
