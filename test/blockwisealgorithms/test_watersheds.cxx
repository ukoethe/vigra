/************************************************************************/
/*                                                                      */
/*     Copyright 2013-2014 by Martin Bidlingmaier and Ullrich Koethe    */
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

#define VIGRA_CHECK_BOUNDS

#include <vigra/blockwise_watersheds.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/unittest.hxx>
#include <vigra/multi_watersheds.hxx>

#include <iostream>
#include <sstream>

#include "utils.hxx"

using namespace std;
using namespace vigra;
using namespace blockwise_watersheds_detail;

struct BlockwiseWatershedTest
{
    void oneDimensionalTest()
    {
        typedef MultiArray<1, unsigned int> Array;
        typedef Array::difference_type Shape;

        vector<Array> data_sets;
        data_sets.push_back(Array(Shape(1)));
        data_sets.push_back(Array(Shape(5)));
        data_sets.push_back(Array(Shape(997)));
        data_sets.push_back(Array(Shape(10000)));
        for(int i = 0; i != data_sets.size(); ++i)
        {
            fillRandom(data_sets[i].begin(), data_sets[i].end(), 3);
        }
        
        for(int i = 0; i != data_sets.size(); ++i)
        {
            const Array& data = data_sets[i];
            
            NeighborhoodType neighborhood = DirectNeighborhood;
            Shape block_shape(1);


            typedef MultiArray<1, size_t> LabelArray;
            
            LabelArray tested_labels(data.shape());
            size_t tested_label_number = unionFindWatershedsBlockwise(data, tested_labels, neighborhood, block_shape);
            
            LabelArray correct_labels(data.shape());
            size_t correct_label_number = watershedsMultiArray(data, correct_labels, neighborhood,
                                                               WatershedOptions().unionFind());
            if(tested_label_number != correct_label_number ||
               !equivalentLabels(tested_labels.begin(), tested_labels.end(),
                                 correct_labels.begin(), correct_labels.end()))
            {
                ostringstream oss;
                oss << "labeling not equivalent" << endl;
                oss << "shape: " << data.shape() << endl;
                failTest(oss.str().c_str());
            }
        }
    }
    void fourDimensionalRandomTest()
    {
        typedef MultiArray<4, unsigned int> Array;
        typedef MultiArray<4, size_t> LabelArray;
        typedef Array::difference_type Shape;
        
        vector<Array> data_sets;
        data_sets.push_back(Array(Shape(1)));
        data_sets.push_back(Array(Shape(2)));
        data_sets.push_back(Array(Shape(4, 2, 1, 5)));
        data_sets.push_back(Array(Shape(4, 5, 7, 3)));
        data_sets.push_back(Array(Shape(6)));
        data_sets.push_back(Array(Shape(1, 10, 100, 1)));

        for(int i = 0; i != data_sets.size(); ++i)
        {
            fillRandom(data_sets[i].begin(), data_sets[i].end(), 3);
        }

        vector<Shape> block_shapes;
        block_shapes.push_back(Shape(1));
        block_shapes.push_back(Shape(2));
        block_shapes.push_back(Shape(1,10,10,2));
        block_shapes.push_back(Shape(1000000));
        block_shapes.push_back(Shape(4,3,10,1000));
        
        vector<NeighborhoodType> neighborhoods;
        neighborhoods.push_back(DirectNeighborhood);
        neighborhoods.push_back(IndirectNeighborhood);

        for(int i = 0; i != data_sets.size(); ++i)
        {
            const Array& data = data_sets[i];
            for(int j = 0; j != block_shapes.size(); ++j)
            {
                const Shape& block_shape = block_shapes[j];
                for(int k = 0; k != neighborhoods.size(); ++k)
                {
                    NeighborhoodType neighborhood = neighborhoods[k];
                    
                    LabelArray tested_labels(data.shape());
                    size_t tested_label_number = unionFindWatershedsBlockwise(data, tested_labels, neighborhood, block_shape);

                    LabelArray correct_labels(data.shape());
                    size_t correct_label_number = watershedsMultiArray(data, correct_labels, neighborhood,
                                                                       WatershedOptions().unionFind());
                    if(tested_label_number != correct_label_number ||
                       !equivalentLabels(tested_labels.begin(), tested_labels.end(),
                                         correct_labels.begin(), correct_labels.end()))
                    {
                        ostringstream oss;
                        oss << "labeling not equivalent" << endl;
                        oss << "array shape: " << data.shape() << endl;
                        oss << "block shape: " << block_shape << endl;
                        oss << "neighborhood: " << neighborhood << endl;
                        oss << "data:" << endl;
                        for(int i = 0; i != data.size(); ++i)
                            oss << data[i] << " ";
                        oss << endl;
                        oss << "expected labels:" << endl;
                        for(int i = 0; i != correct_labels.size(); ++i)
                            oss << correct_labels[i] << " ";
                        oss << endl;
                        oss << "got" << endl;
                        for(int i = 0; i != tested_labels.size(); ++i)
                            oss << tested_labels[i] << " ";
                        oss << endl;
                        failTest(oss.str().c_str());
                    }
                }
            }
        }
    }
    void chunkedTest()
    {
        typedef MultiArray<3, int> OldschoolArray;
        typedef MultiArray<3, size_t> OldschoolLabelArray;

        typedef ChunkedArrayLazy<3, int> Array;
        typedef ChunkedArrayLazy<3, size_t> LabelArray;
        
        typedef OldschoolArray::difference_type Shape;
        
        Shape shape(20, 30, 10);
        Shape chunk_shape(20, 30, 10);
        NeighborhoodType neighborhood = IndirectNeighborhood;

        OldschoolArray oldschool_data(shape);

        fillRandom(oldschool_data.begin(), oldschool_data.end(), 3);
        OldschoolLabelArray correct_labels(shape);
        size_t correct_label_number = watershedsMultiArray(oldschool_data, correct_labels, neighborhood,
                                                           WatershedOptions().unionFind());
        
        Array data(shape);
        data.commitSubarray(Shape(0), oldschool_data);
        LabelArray tested_labels(shape);
        size_t tested_label_number = unionFindWatershedsBlockwise(data, tested_labels, neighborhood);
        shouldEqual(correct_label_number, tested_label_number);
        shouldEqual(equivalentLabels(tested_labels.begin(), tested_labels.end(),
                                     correct_labels.begin(), correct_labels.end()),
                    true);
    }
};

struct BlockwiseWatershedTestSuite
  : public test_suite
{
    BlockwiseWatershedTestSuite()
      : test_suite("blockwise watershed test")
    {
        add(testCase(&BlockwiseWatershedTest::fourDimensionalRandomTest));
        add(testCase(&BlockwiseWatershedTest::oneDimensionalTest));
        add(testCase(&BlockwiseWatershedTest::chunkedTest));
    }
};

int main(int argc, char** argv)
{
    BlockwiseWatershedTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    cout << test.report() << endl;

    return failed != 0;
}
