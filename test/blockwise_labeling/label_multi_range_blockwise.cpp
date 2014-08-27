#include <vigra/multi_array.hxx>
#include <vigra/range/multi_range.hxx>

#include <vigra/label_multi_range_blockwise.hxx>

#include <vigra/range/multi_blockify.hxx>
#include <vigra/range/multi_range.hxx>
#include <vigra/range/multi_range_traits.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/error.hxx>
#include <vigra/multi_labeling.hxx>

#include <iostream>
#include <vector>

using namespace vigra;
using namespace std;

template <class LabelRange1,class LabelRange2>
bool equivalent_labels(LabelRange1 range1, LabelRange2 range2)
{
    typedef typename multi_range_traits::value_type<LabelRange1>::type Label;
    typedef vector<typename std::decay<Label>::type> LabelMap;
    LabelMap left_to_right;
    LabelMap right_to_left;

    bool result = true;
    for_each_flattened(multi_zip(range1, range2), [&](std::pair<Label, Label> p)
    {
        Label left = p.first;
        Label right = p.second;
        left_to_right.resize(max<size_t>(left_to_right.size(), left + 1));
        right_to_left.resize(max<size_t>(right_to_left.size(), right + 1));

        if(left_to_right[left] == 0)
            left_to_right[left] = right;
        else if(left_to_right[left] != right)
            result = false;
        
        if(right_to_left[right] == 0)
            right_to_left[right] = left;
        else if(right_to_left[right] != left)
            result = false;
    });

    return result;
}

int main()
{
  try {
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
    vigra_assert(count == blockwise_count, "");
    vigra_assert(equivalent_labels(multi_range(labels), multi_range(blockwise_labels)), "labelings differ");
  }
  catch(std::exception & e)
  {
      std::cerr << "Program failed:\n" << e.what() << "\n";
  }
}
