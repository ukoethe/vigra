/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Ullrich Koethe and Philip Schill       */
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
#include <vigra/unittest.hxx>
#include <vigra/random_forest_3.hxx>
#include <vigra/random.hxx>
#ifdef HasHDF5
    #include <vigra/random_forest_3_hdf5_impex.hxx>
#endif

using namespace vigra;
using namespace vigra::rf3;

struct RandomForestTests
{
    void test_base_class()
    {
        typedef BinaryForest Graph;
        typedef Graph::Node Node;
        typedef LessEqualSplitTest<double> SplitTest;
        typedef ArgMaxAcc Acc;
        typedef RandomForest<MultiArray<2, double>, MultiArray<1, int>, SplitTest, Acc> RF;

        // Build a forest from scratch.
        Graph gr;
        RF::NodeMap<SplitTest>::type split_tests;
        RF::NodeMap<size_t>::type leaf_responses;
        {
            Node n0 = gr.addNode();
            Node n1 = gr.addNode();
            Node n2 = gr.addNode();
            Node n3 = gr.addNode();
            Node n4 = gr.addNode();
            Node n5 = gr.addNode();
            Node n6 = gr.addNode();
            gr.addArc(n0, n1);
            gr.addArc(n0, n2);
            gr.addArc(n1, n3);
            gr.addArc(n1, n4);
            gr.addArc(n2, n5);
            gr.addArc(n2, n6);

            split_tests.insert(n0, SplitTest(0, 0.6));
            split_tests.insert(n1, SplitTest(1, 0.25));
            split_tests.insert(n2, SplitTest(1, 0.75));
            leaf_responses.insert(n3, 0);
            leaf_responses.insert(n3, 0);
            leaf_responses.insert(n4, 1);
            leaf_responses.insert(n5, 2);
            leaf_responses.insert(n6, 3);
        }
        std::vector<int> distinct_labels;
        distinct_labels.push_back(0);
        distinct_labels.push_back(1);
        distinct_labels.push_back(-7);
        distinct_labels.push_back(3);
        auto const pspec = ProblemSpec<int>().num_features(2).distinct_classes(distinct_labels);
        RF rf = RF(gr, split_tests, leaf_responses, pspec);

        // Check if the given points are predicted correctly.
        double test_x_values[] = {
            0.2, 0.4, 0.2, 0.4, 0.7, 0.8, 0.7, 0.8,
            0.2, 0.2, 0.7, 0.7, 0.2, 0.2, 0.8, 0.8
        };
        MultiArray<2, double> test_x(Shape2(8, 2), test_x_values);
        int test_y_values[] = {
            0, 0, 1, 1, -7, -7, 3, 3
        };
        MultiArray<1, int> test_y(Shape1(8), test_y_values);
        MultiArray<1, int> pred_y(Shape1(8));
        rf.predict(test_x, pred_y, 1);
        shouldEqualSequence(pred_y.begin(), pred_y.end(), test_y.begin());
    }

    void test_default_rf()
    {
        typedef MultiArray<2, double> Features;
        typedef MultiArray<1, int> Labels;

        double train_x_values[] = {
            0.2, 0.4, 0.2, 0.4, 0.7, 0.8, 0.7, 0.8,
            0.2, 0.2, 0.7, 0.7, 0.2, 0.2, 0.8, 0.8
        };
        Features train_x(Shape2(8, 2), train_x_values);
        int train_y_values[] = {
            0, 0, 1, 1, -7, -7, 3, 3
        };
        Labels train_y(Shape1(8), train_y_values);
        Features test_x(train_x);
        Labels test_y(train_y);

        std::vector<RandomForestOptionTags> splits;
        splits.push_back(RF_GINI);
        splits.push_back(RF_ENTROPY);
        splits.push_back(RF_KSD);
        for (auto split : splits)
        {
            RandomForestOptions const options = RandomForestOptions()
                                                       .tree_count(1)
                                                       .bootstrap_sampling(false)
                                                       .split(split)
                                                       .n_threads(1);
            auto rf = random_forest(train_x, train_y, options);
            Labels pred_y(test_y.shape());
            rf.predict(test_x, pred_y, 1);
            shouldEqualSequence(pred_y.begin(), pred_y.end(), test_y.begin());
        }
    }

    void test_oob_visitor()
    {
        // Create a (noisy) grid with datapoints and assign classes as in a 4x4 chessboard.
        size_t const nx = 100;
        size_t const ny = 100;

        RandomNumberGenerator<MersenneTwister> rand;
        MultiArray<2, double> train_x(Shape2(nx*ny, 2));
        MultiArray<1, int> train_y(Shape1(nx*ny));
        for (size_t y = 0; y < ny; ++y)
        {
            for (size_t x = 0; x < nx; ++x)
            {
                train_x(y*nx+x, 0) = x + 2*rand.uniform()-1;
                train_x(y*nx+x, 1) = y + 2*rand.uniform()-1;
                if ((x/25+y/25) % 2 == 0)
                    train_y(y*nx+x) = 0;
                else
                    train_y(y*nx+x) = 1;
            }
        }

        RandomForestOptions const options = RandomForestOptions()
                                                   .tree_count(10)
                                                   .bootstrap_sampling(true)
                                                   .n_threads(1);
        OOBError oob;
        auto rf = random_forest(train_x, train_y, options, create_visitor(oob));
        should(oob.oob_err_ > 0.02 && oob.oob_err_ < 0.04); // FIXME: Use a statistical approach here.
    }
    
    void test_var_importance_visitor()
    {
        // Create a (noisy) grid with datapoints and split the classes according to an oblique line.
        size_t const nx = 20;
        size_t const ny = 20;

        RandomNumberGenerator<MersenneTwister> rand;
        MultiArray<2, double> train_x(Shape2(nx*ny, 2));
        MultiArray<1, int> train_y(Shape1(nx*ny));
        for (size_t y = 0; y < ny; ++y)
        {
            for (size_t x = 0; x < nx; ++x)
            {
                train_x(y*nx+x, 0) = x + 2*rand.uniform()-1;
                train_x(y*nx+x, 1) = y + 2*rand.uniform()-1;
                if (x - nx/2.0 + 4*y - 4*ny/2.0  <  0)
                    train_y(y*nx+x) = 0;
                else
                    train_y(y*nx+x) = 1;
            }
        }

        RandomForestOptions const options = RandomForestOptions()
                                                   .tree_count(10)
                                                   .bootstrap_sampling(true)
                                                   .n_threads(1);
        VariableImportance var_imp;
        auto rf = random_forest(train_x, train_y, options, create_visitor(var_imp));

        // The permutation importances of feature 1 should be about 
        // 10 times as big as the importances of feature 0.
        for (size_t i = 0; i < 4; ++i)
        {
            should(var_imp.variable_importance_(1, i) > 5 * var_imp.variable_importance_(0, i));
        }
    }

#ifdef HasHDF5
    void test_import()
    {
        typedef float FeatureType;
        typedef UInt32 LabelType;
        typedef MultiArray<2, FeatureType> Features;
        typedef MultiArray<1, LabelType> Labels;

        // Load the dummy random forest. This RF was trained on 500
        // 2-dimensional points in [0, 1]^2. All points with x<0.5
        // and y<0.5 have class 0, points with x<0.5 and y>=0.5 have
        // class 1, points with x>0.5 and y>=0.5 have class 2, and
        // points with x>=0.5 and y<0.5 have class 3.
        HDF5File hfile("data/rf.h5", HDF5File::ReadOnly);
        auto rf = random_forest_import_HDF5<Features, Labels>(hfile);
        
        // Create some test data.
        FeatureType test_x_data[] = {
            0.2f, 0.4f, 0.6f, 0.8f,  0.2f, 0.4f, 0.6f, 0.8f,  0.2f, 0.4f, 0.6f, 0.8f,  0.2f, 0.4f, 0.6f, 0.8f,
            0.2f, 0.2f, 0.2f, 0.2f,  0.4f, 0.4f, 0.4f, 0.4f,  0.6f, 0.6f, 0.6f, 0.6f,  0.8f, 0.8f, 0.8f, 0.8f
        };
        Features test_x(Shape2(16, 2), test_x_data);
        LabelType test_y_data[] = {
            0, 0, 3, 3,  0, 0, 3, 3,  1, 1, 2, 2,  1, 1, 2, 2
        };
        Labels test_y(Shape1(16), test_y_data);
        Labels pred_y(Shape1(16));

        // Use the RF to predict the data.
        rf.predict(test_x, pred_y);
        for (size_t i = 0; i < (size_t)test_y.size(); ++i)
            should(test_y(i) == pred_y(i));
    }

    void test_export()
    {
        typedef float FeatureType;
        typedef UInt32 LabelType;
        typedef MultiArray<2, FeatureType> Features;
        typedef MultiArray<1, LabelType> Labels;

        // Load the dummy random forest.
        HDF5File infile("data/rf.h5", HDF5File::ReadOnly);
        auto rf = random_forest_import_HDF5<Features, Labels>(infile);

        // Save the dummy random forest.
        HDF5File outfile("data/rf_out.h5", HDF5File::New);
        random_forest_export_HDF5(rf, outfile);
    }
#endif
};

struct RandomForestTestSuite : public test_suite
{
    RandomForestTestSuite()
        :
        test_suite("RandomForest test")
    {
        add(testCase(&RandomForestTests::test_base_class));
        add(testCase(&RandomForestTests::test_default_rf));
        add(testCase(&RandomForestTests::test_oob_visitor));
        add(testCase(&RandomForestTests::test_var_importance_visitor));
#ifdef HasHDF5
        add(testCase(&RandomForestTests::test_import));
        add(testCase(&RandomForestTests::test_export));
#endif
    }
};

int main(int argc, char** argv)
{
    RandomForestTestSuite forest_test;
    int failed = forest_test.run(testsToBeExecuted(argc, argv));
    std::cout << forest_test.report() << std::endl;
    return (failed != 0);
}
