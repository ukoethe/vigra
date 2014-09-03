#include <vigra/blockwise_labeling.hxx>


#include <vigra/multi_array.hxx>
#include <vigra/multi_labeling.hxx>
#include <vigra/unittest.hxx>

#include <vector>
#include <utility>
#include <functional>

using namespace vigra;
using namespace std;

template <class Iterator1,class Iterator2>
bool equivalentLabels(Iterator1 begin1, Iterator1 end1,
                      Iterator2 begin2, Iterator2 end2)
{
    if(end1 - begin1 != end2 - begin2)
        return false;
    
    typedef vector<int> LabelMap;
    LabelMap left_to_right; // from range 1 to range 2
    LabelMap right_to_left;
    for( ; begin1 != end1; ++begin1, ++begin2)
    {
        if(left_to_right.size() <= *begin1)
            left_to_right.resize(*begin1 + 1, -1); // "-1" means unmapped
        if(right_to_left.size() <= *begin2)
            right_to_left.resize(*begin2 + 1, -1);

        if(left_to_right[*begin1] == -1) // unmapped -> map it
            left_to_right[*begin1] = *begin2;
        else if(left_to_right[*begin1] != *begin2) // already mapped to different value -> not equivalent labels
            return false;

        if(right_to_left[*begin2] == -1)
            right_to_left[*begin2] = *begin1;
        else if(right_to_left[*begin2] != *begin1)
            return false;
    }
    
    return true;
}

template <class Iterator>
void fill_random(Iterator begin, Iterator end, int maximum)
{
    for( ; begin != end; ++begin)
        *begin = rand() % maximum;
}

template <class DatasIterator, class ShapesIterator>
void test_on_data(DatasIterator datas_begin, DatasIterator datas_end,
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
            //neighborhoods.push_back(IndirectNeighborhood);
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
                    
                    int correct_label_number;
                    int tested_label_number;
                    if(with_background)
                    {
                        correct_label_number = labelMultiArrayWithBackground(data, correct_labels, neighborhood, 1u);
                        tested_label_number = labelMultiArrayWithBackgroundBlockwise(data, tested_labels, neighborhood, 1u, shape);
                    }
                    else
                    {
                        correct_label_number = labelMultiArray(data, correct_labels, neighborhood);
                        tested_label_number = labelMultiArrayBlockwise(data, tested_labels, neighborhood, shape);
                    }

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
    typedef MultiArray<4, unsigned int> Array4;
    typedef MultiArray<1, unsigned int> Array1;

    typedef Array5::difference_type Shape5;
    typedef Array4::difference_type Shape4;
    typedef Array1::difference_type Shape1;

    vector<Array5> array_fives;
    vector<Array4> array_fours;
    vector<Array1> array_ones;

    vector<Shape5> shape_fives;
    vector<Shape4> shape_fours;
    vector<Shape1> shape_ones;

    BlockwiseLabelingTest()
    {
        array_fives.push_back(Array5(Shape5(1)));
        array_fives.push_back(Array5(Shape5(2,2,3,4,3)));
        array_fives.push_back(Array5(Shape5(5,6,3,2,3)));
        for(int i = 0; i != array_fives.size(); ++i)
        {
            fill_random(array_fives[i].begin(), array_fives[i].end(), 3);
        }

        array_fours.push_back(Array4(Shape4(1)));
        array_fours.push_back(Array4(Shape4(1,2,3,4)));
        array_fours.push_back(Array4(Shape4(6,10,8,7)));
        for(int i = 0; i != array_fours.size(); ++i)
        {
            fill_random(array_fours[i].begin(), array_fours[i].end(), 3);
        }
        
        array_ones.push_back(Array1(Shape1(1)));
        array_ones.push_back(Array1(Shape1(2)));
        array_ones.push_back(Array1(Shape1(47)));
        array_ones.push_back(Array1(Shape1(81)));
        array_ones.push_back(Array1(Shape1(997)));
        
        for(int i = 0; i != array_ones.size(); ++i)
        {
            fill_random(array_ones[i].begin(), array_ones[i].end(), 3);
        }

        shape_fives.push_back(Shape5(1));
        shape_fives.push_back(Shape5(2));
        shape_fives.push_back(Shape5(5,4,3,2,1));
        shape_fives.push_back(Shape5(100000));

        shape_fours.push_back(Shape4(1));
        shape_fours.push_back(Shape4(2));
        shape_fours.push_back(Shape4(5,4,2,1));
        shape_fours.push_back(Shape4(1000000));

        shape_ones.push_back(Shape1(1));
        shape_ones.push_back(Shape1(2));
        shape_ones.push_back(Shape1(3));
        shape_ones.push_back(Shape1(5));
        shape_ones.push_back(Shape1(213));
    }
    
    void debug_test()
    {
        typedef MultiArray<2, int> Array;
        typedef Array::difference_type Shape;

        Shape shape = Shape(4);
        Array data(shape);

        data(0,0) = 1;
        data(1,1) = 2;
        data(2,1) = 2;
        data(2,2) = 2;
        data(3,3) = 3;
        data(3,2) = 3;
        
        Array blockwise_labels(shape);
        Array labels(shape);
    
        Shape block_shape(2, 2);

        NeighborhoodType neighborhood = IndirectNeighborhood;
    
        size_t count = labelMultiArrayWithBackground(data, labels, neighborhood);
        size_t blockwise_count = labelMultiArrayWithBackgroundBlockwise(data, blockwise_labels, neighborhood, 0, block_shape);
        shouldEqual(count, blockwise_count);
        shouldEqual(equivalentLabels(labels.begin(), labels.end(),
                                     blockwise_labels.begin(), blockwise_labels.end()),
                    true);
    }

    void five_dimensional_test()
    {
        test_on_data(array_fives.begin(), array_fives.end(),
                     shape_fives.begin(), shape_fives.end());
    }
    void four_dimensional_test()
    {
        test_on_data(array_fours.begin(), array_fours.end(),
                     shape_fours.begin(), shape_fours.end());
    }
    void one_dimensional_test()
    {
        test_on_data(array_ones.begin(), array_ones.end(),
                     shape_ones.begin(), shape_ones.end());
    }
};

struct BlockwiseLabelingTestSuite
  : public test_suite
{
    BlockwiseLabelingTestSuite()
      : test_suite("blockwise labeling test")
    {
        add(testCase(&BlockwiseLabelingTest::five_dimensional_test));
        add(testCase(&BlockwiseLabelingTest::four_dimensional_test));
        add(testCase(&BlockwiseLabelingTest::one_dimensional_test));
        add(testCase(&BlockwiseLabelingTest::debug_test));
    }
};

int main(int argc, char** argv)
{
    BlockwiseLabelingTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    cout << test.report() << endl;

    return failed != 0;
}

