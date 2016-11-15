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

#include <vigra/blockwise_labeling.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/multi_array_chunked.hxx>
#include <vigra/multi_labeling.hxx>
#include <vigra/unittest.hxx>

#include <vector>
#include <utility>
#include <functional>

#include "utils.hxx"

using namespace vigra;
using namespace std;

template <unsigned int N, class DatasIterator, class ShapesIterator>
void testOnData(DatasIterator datas_begin, DatasIterator datas_end,
                  ShapesIterator shapes_begin, ShapesIterator shapes_end)
{
    for(DatasIterator datas_it = datas_begin ; datas_it != datas_end; ++datas_it)
    {
        typedef typename DatasIterator::reference DataRef;
        DataRef data = *datas_it;
        for(ShapesIterator shapes_it = shapes_begin; shapes_it != shapes_end; ++shapes_it)
        {
            typedef typename ShapesIterator::reference ShapeRef;
            ShapeRef shape = *shapes_it;

            vector<NeighborhoodType> neighborhoods;
            neighborhoods.push_back(DirectNeighborhood);
            neighborhoods.push_back(IndirectNeighborhood);
            typedef vector<NeighborhoodType>::iterator NeighborhoodIterator;

            for(NeighborhoodIterator neighborhood_it = neighborhoods.begin();
                neighborhood_it != neighborhoods.end();
                ++neighborhood_it)
            {
                NeighborhoodType neighborhood = *neighborhood_it;
                typedef typename DatasIterator::value_type Data;

                std::vector<bool> with_backgrounds;
                with_backgrounds.push_back(true);
                with_backgrounds.push_back(false);
                for(std::vector<bool>::iterator backgrounds_it = with_backgrounds.begin();
                        backgrounds_it != with_backgrounds.end();
                        ++backgrounds_it)
                {
                    bool with_background = *backgrounds_it;
                    Data correct_labels(data.shape());
                    Data tested_labels(data.shape());

                    BlockwiseLabelOptions<N> options;
                    options.neighborhood(neighborhood);
                    options.blockShape(shape);

                    int correct_label_number;
                    int tested_label_number;
                    if(with_background)
                    {
                        correct_label_number = labelMultiArrayWithBackground(data, correct_labels, neighborhood, 1u);
                        options.ignoreBackgroundValue(1u);
                    }
                    else
                    {
                        correct_label_number = labelMultiArray(data, correct_labels, neighborhood);
                    }

                    tested_label_number = labelMultiArrayBlockwise(data, tested_labels, options);

                    if(!equivalentLabels(correct_labels.begin(), correct_labels.end(),
                                         tested_labels.begin(), tested_labels.end()) ||
                       correct_label_number != tested_label_number)
                    {
                        std::ostringstream oss;
                        oss << "labeling not equivalent" << endl;
                        oss << "with background: " << boolalpha << with_background << endl;
                        oss << "array shape: " << data.shape() << endl;
                        oss << "block shape: " << shape << endl;
                        oss << "neighborhood: " << neighborhood << endl;
                        oss << "data: " << endl;
                        for(int i = 0; i != data.size(); ++i)
                            oss << data[i] << " ";
                        oss << endl;
                        oss << "expected_labels: " << endl;
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
}

struct BlockwiseLabelingTest
{
    typedef MultiArray<5, unsigned int> Array5;
    typedef MultiArray<2, unsigned int> Array2;
    typedef MultiArray<1, unsigned int> Array1;

    typedef Array5::difference_type Shape5;
    typedef Array2::difference_type Shape2;
    typedef Array1::difference_type Shape1;

    vector<Array5> array_fives;
    vector<Array2> array_twos;
    vector<Array1> array_ones;

    vector<Shape5> shape_fives;
    vector<Shape2> shape_twos;
    vector<Shape1> shape_ones;

    BlockwiseLabelingTest()
    {
        array_fives.push_back(Array5(Shape5(1)));
        array_fives.push_back(Array5(Shape5(2,2,3,4,3)));
        array_fives.push_back(Array5(Shape5(5,6,2,2,3)));
        for(decltype(array_fives.size()) i = 0; i != array_fives.size(); ++i)
        {
            fillRandom(array_fives[i].begin(), array_fives[i].end(), 3);
        }

        array_twos.push_back(Array2(Shape2(1)));
        array_twos.push_back(Array2(Shape2(1,2)));
        array_twos.push_back(Array2(Shape2(2,2)));
        array_twos.push_back(Array2(Shape2(4,4)));
        array_twos.push_back(Array2(Shape2(6,10)));
        array_twos.push_back(Array2(Shape2(19,25)));
        for(decltype(array_twos.size()) i = 0; i != array_twos.size(); ++i)
        {
            fillRandom(array_twos[i].begin(), array_twos[i].end(), 3);
        }

        array_ones.push_back(Array1(Shape1(1)));
        array_ones.push_back(Array1(Shape1(2)));
        array_ones.push_back(Array1(Shape1(47)));
        array_ones.push_back(Array1(Shape1(81)));
        array_ones.push_back(Array1(Shape1(997)));
        for(decltype(array_ones.size()) i = 0; i != array_ones.size(); ++i)
        {
            fillRandom(array_ones[i].begin(), array_ones[i].end(), 3);
        }

        shape_fives.push_back(Shape5(1));
        shape_fives.push_back(Shape5(2));
        shape_fives.push_back(Shape5(5,4,3,2,1));
        shape_fives.push_back(Shape5(100000));

        shape_twos.push_back(Shape2(1));
        shape_twos.push_back(Shape2(2));
        shape_twos.push_back(Shape2(5,4));
        shape_twos.push_back(Shape2(1000000));

        shape_ones.push_back(Shape1(1));
        shape_ones.push_back(Shape1(2));
        shape_ones.push_back(Shape1(3));
        shape_ones.push_back(Shape1(5));
        shape_ones.push_back(Shape1(213));
    }

    void debugTest()
    {
        typedef MultiArray<2, int> Array;
        typedef Array::difference_type Shape;

        Shape shape = Shape(2);
        Array data(shape);

        data(0,0) = 1;
        data(1,0) = 1;
        data(0,1) = 1;
        data(1,1) = 1;

        MultiArray<2, size_t> blockwise_labels(shape);
        MultiArray<2, size_t> labels(shape);

        Shape block_shape(1, 1);
        //TinyVector<MultiArrayIndex, 3> block_shape(1, 1, 1);

        NeighborhoodType neighborhood = IndirectNeighborhood;
        using namespace vigra::blockwise;

        size_t count = labelMultiArrayWithBackground(data, labels, neighborhood, 1);
        size_t blockwise_count = labelMultiArrayBlockwise(data, blockwise_labels,
                                                          BlockwiseLabelOptions<2>().neighborhood(neighborhood)
                                                                                 .ignoreBackgroundValue(1)
                                                                                 .blockShape(block_shape));
        shouldEqual(count, blockwise_count);
        shouldEqual(equivalentLabels(labels.begin(), labels.end(),
                                     blockwise_labels.begin(), blockwise_labels.end()),
                    true);
    }
    void chunkedArrayTest()
    {
        typedef ChunkedArrayLazy<3, int> DataArray;
        typedef ChunkedArrayLazy<3, size_t> LabelArray;
        typedef DataArray::shape_type Shape;

        Shape shape = Shape(100, 200, 300);
        Shape chunk_shape = Shape(32);

        DataArray data(shape, chunk_shape);
        fillRandom(data.begin(), data.end(), 3);
        LabelArray labels(shape, chunk_shape);

        using namespace vigra::blockwise;

        BlockwiseLabelOptions<3> options;
        options.neighborhood(IndirectNeighborhood).ignoreBackgroundValue(1);
        options.blockShape(Shape(40));

        size_t tested_label_number = labelMultiArrayBlockwise(data, labels, options);
        MultiArray<3, size_t> checked_out_labels(shape);
        labels.checkoutSubarray(Shape(0), checked_out_labels);

        MultiArray<3, int> oldschool_data_array(shape);
        MultiArray<3, size_t> oldschool_label_array(shape);

        data.checkoutSubarray(Shape(0), oldschool_data_array);
        size_t actual_label_number = labelMultiArrayWithBackground(oldschool_data_array, oldschool_label_array,
                                                                   IndirectNeighborhood, 1);

        shouldEqual(tested_label_number, actual_label_number);
        shouldEqual(equivalentLabels(checked_out_labels.begin(), checked_out_labels.end(),
                                     oldschool_label_array.begin(), oldschool_label_array.end()), true);
    }

    void fiveDimensionalRandomTest()
    {
        testOnData<5>(array_fives.begin(), array_fives.end(),
                     shape_fives.begin(), shape_fives.end());
    }
    void twoDimensionalRandomTest()
    {
        testOnData<2>(array_twos.begin(), array_twos.end(),
                     shape_twos.begin(), shape_twos.end());
    }
    void oneDimensionalRandomTest()
    {
        testOnData<1>(array_ones.begin(), array_ones.end(),
                     shape_ones.begin(), shape_ones.end());
    }
};

struct BlockwiseLabelingTestSuite
  : public test_suite
{
    BlockwiseLabelingTestSuite()
      : test_suite("blockwise labeling test")
    {
        add(testCase(&BlockwiseLabelingTest::oneDimensionalRandomTest));
        add(testCase(&BlockwiseLabelingTest::twoDimensionalRandomTest));
        add(testCase(&BlockwiseLabelingTest::fiveDimensionalRandomTest));
        add(testCase(&BlockwiseLabelingTest::debugTest));
        add(testCase(&BlockwiseLabelingTest::chunkedArrayTest));
    }
};

int main(int argc, char** argv)
{
    BlockwiseLabelingTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    cout << test.report() << endl;

    return failed != 0;
}

