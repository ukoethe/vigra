#include <vigra/label_multi_range_blockwise.hxx>

#include <vigra/range/multi_range.hxx>
#include <vigra/range/multi_blockify.hxx>

#include <vigra/range/multi_range.hxx>
#include <vigra/range/multi_range_traits.hxx>

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
    LabelMap mapping; //from range 1 to range 2
    for( ; begin1 != end1; ++begin1, ++begin2)
    {
        if(mapping.size() <= *begin1)
            mapping.resize(*begin1 + 1, -1);

        if(mapping[*begin1] == -1)
            mapping[*begin1] = *begin2;
        else if(mapping[*begin1] != *begin2)
            return false;
    }
    
    return true;
}

template <unsigned int N>
typename MultiArrayShape<N>::type random_coordinates(const typename MultiArrayShape<N>::type& within)
{
    typename MultiArrayShape<N>::type result;
    for(int i = 0; i != N; ++i)
        result[i] = rand() % within[i];
    return result;
}

template <unsigned int N>
typename MultiArrayShape<N>::type random_neighbor_offset(NeighborhoodType neighborhood)
{
    typename MultiArrayShape<N>::type result = 0;
    if(neighborhood == DirectNeighborhood)
    {
        unsigned int dim = rand() % N;
        unsigned int direction = rand() % 2;
        if(direction == 0)
            result[dim] = - 1;
        else
            result[dim] = 1;
    }
    else
    {
        for(int i = 0; i != N; ++i)
            result[i] = (rand() % 3) - 1;
    }
    return result;
}
template <unsigned int N>
void set_within(typename MultiArrayShape<N>::type& coordinates, typename MultiArrayShape<N>::type bounds)
{
    for(int i = 0; i != N; ++i)
    {
        coordinates[i] = max(coordinates[i], 0);
        coordinates[i] = min(coordinates[i], bounds[i] - 1);
    }
}

template <unsigned int N, class S>
void create_random_path(MultiArrayView<N, int, S> view, NeighborhoodType neighborhood, unsigned int length, int value)
{
    typedef typename MultiArrayShape<N>::type Shape;
    
    Shape path = random_coordinates(view.shape());
    for(unsigned int i = 0; i != length; ++i)
    {
        view[path] = value;
        path += random_neighbor_offset<N>(neighborhood);
        set_within(path);
    }
}

template <class Iterator>
void fill_random(Iterator begin, Iterator end, int maximum)
{
    for( ; begin != end; ++begin)
        *begin = rand() % maximum;
}

template <unsigned int N, class T, class S, class LabelFunctor1, class LabelFunctor2>
bool equivalentFunctors(MultiArrayView<N, T, S> data, LabelFunctor1 functor1, LabelFunctor2 functor2)
{
    MultiArray<N, T> first_labels(data.shape());
    functor1(data, first_labels);

    MultiArray<N, T> second_labels(data.shape());
    functor2(data, second_labels);
    
    return equivalent_labels(first_labels.begin(), first_labels.end(),
                             second_labels.begin(), second_labels.end());
}

template <class DatasIterator, class ShapesIterator>
void test_on_data(DatasIterator datas_begin, DatasIterator datas_end,
               ShapesIterator shapes_begin, ShapesIterator shapes_end)
{
    typedef typename DatasIterator::reference DataRef;
    std::for_each(datas_begin, datas_end, [&](DataRef data)
    {
        typedef typename ShapesIterator::reference ShapeRef;
        std::for_each(shapes_begin, shapes_end, [&](ShapeRef shape)
        {
            vector<NeighborhoodType> neighborhoods;
            neighborhoods.push_back(DirectNeighborhood);
            neighborhoods.push_back(IndirectNeighborhood);
            std::for_each(neighborhoods.begin(), neighborhoods.end(), [&](NeighborhoodType neighborhood)
            {
                typedef typename DatasIterator::value_type Data;
                
                Data correct_labels(data.shape());
                int correct_label_number = labelMultiArray(data, correct_labels, neighborhood);
                
                Data tested_labels(data.shape());
                auto data_blocks = multi_blockify(multi_range(data), shape);
                auto tested_labels_blocks = multi_blockify(multi_range(tested_labels), shape);
                int tested_label_number = label_multi_range_blockwise(data_blocks, tested_labels_blocks, neighborhood, std::equal_to<int>());
                if(!equivalentLabels(correct_labels.begin(), correct_labels.end(),
                                     tested_labels.begin(), tested_labels.end()) ||
                   correct_label_number != tested_label_number)
                {
                    std::ostringstream oss;
                    oss << "labeling not equivalent" << endl;
                    oss << "array shape: " << data.shape() << endl;
                    oss << "block shape: " << shape << endl;
                    failTest(oss.str().c_str());
                }
            });
        });
    });
}

struct BlockwiseLabelingTest
{   
    typedef MultiArray<5, int> Array5;
    typedef MultiArray<4, int> Array4;
    typedef MultiArray<1, int> Array1;

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
        //array_fives.push_back(Array5(Shape5(5,6,3,2,3)));
        std::for_each(array_fives.begin(), array_fives.end(), [](Array5& arr)
        {
            fill_random(arr.begin(), arr.end(), 3);
        });

        array_fours.push_back(Array4(Shape4(1)));
        array_fours.push_back(Array4(Shape4(1,2,3,4)));
        array_fours.push_back(Array4(Shape4(6,10,8,7)));
        std::for_each(array_fours.begin(), array_fours.end(), [](Array4& arr)
        {
            fill_random(arr.begin(), arr.end(), 3);
        });
        
        array_ones.push_back(Array1(Shape1(1)));
        array_ones.push_back(Array1(Shape1(2)));
        array_ones.push_back(Array1(Shape1(47)));
        array_ones.push_back(Array1(Shape1(81)));
        
        std::for_each(array_ones.begin(), array_ones.end(), [](Array1& arr)
        {
            fill_random(arr.begin(), arr.end(), 3);
        });

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
        typedef MultiArray<2, size_t> Array;
        typedef Array::difference_type Shape;

        Shape shape = Shape(4);
        Array data(shape);

        data(0,0) = 1;
        data(1,1) = 2;
        data(2,1) = 2;
        data(2,2) = 2;
        data(3,3) = 3;
        data(3,2) = 3;

        auto data_range = multi_range(data);
        Array blockwise_labels(shape);
        auto label_range = multi_range(blockwise_labels);
    
        Shape block_shape(2, 2);
        auto data_blocks_range = multi_blockify(data_range, block_shape);
        auto label_blocks_range = multi_blockify(label_range, block_shape);

        NeighborhoodType neighborhood = IndirectNeighborhood;
        auto equal = [](size_t a, size_t b)
        {
            return a == b;
        };
    
        Array labels(shape);
        size_t count = labelMultiArray(data, labels, neighborhood, equal);
        size_t blockwise_count = label_multi_range_blockwise(data_blocks_range, label_blocks_range, neighborhood, equal);
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

    typedef MultiArray<2, size_t> Array;
    typedef Array::difference_type Shape;

    Shape shape = Shape(4);
    Array data(shape);
    
    data(0,0) = 1;
    data(1,1) = 2;
    data(2,1) = 2;
    data(2,2) = 2;
    data(3,3) = 3;
    data(3,2) = 3;

    auto data_range = multi_range(data);
    Array blockwise_labels(shape);
    auto label_range = multi_range(blockwise_labels);
    
    Shape block_shape(2, 1);
    auto data_blocks_range = multi_blockify(data_range, block_shape);
    auto label_blocks_range = multi_blockify(label_range, block_shape);

    NeighborhoodType neighborhood = IndirectNeighborhood;
    auto equal = [](size_t a, size_t b)
    {
        return a == b;
    };
    
    Array labels(shape);
    size_t count = labelMultiArray(data, labels, neighborhood, equal);
    size_t blockwise_count = label_multi_range_blockwise(data_blocks_range, label_blocks_range, neighborhood, equal);
    //vigra_assert(count == blockwise_count, "");
    //vigra_assert(equivalent_labels(multi_range(labels), multi_range(blockwise_labels)), "");
}
