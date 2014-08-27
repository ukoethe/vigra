#include <label_multi_range.hxx>

#include <vigra/range/multi_range.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/error.hxx>

#include <iostream>

using namespace vigra;
using namespace std;


int main()
{
    typedef MultiArray<2, int> Array;
    typedef Array::difference_type Shape;

    Array data(Shape(4));
    Array labels(Shape(4));
    
    data.subarray(Shape(0,0), Shape(2,2)) = 1;
    data.subarray(Shape(2,2), Shape(4,4)) = 2;
    int highest_label = labelMultiRange(multi_range(data), multi_range(labels),
                                        DirectNeighborhood, 
                                        [](int a, int b){ return a == b; });
    vigra_assert(highest_label == 4, "");
    vigra_assert(labels(1,1) != labels(2,2), "");
    vigra_assert(labels(1,1) != labels(2,1), "");
    vigra_assert(labels(1,1) != labels(1,2), "");
    vigra_assert(labels(0,0) == labels(1,1), "");
}
