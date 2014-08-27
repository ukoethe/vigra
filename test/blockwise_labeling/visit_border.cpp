#include <visit_border.hxx>

#include <vigra/range/multi_range.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/union_find.hxx>

#include <vector>

using namespace vigra;
using namespace std;

int main()
{
    typedef MultiArray<2, int> Array;
    typedef Array::difference_type Shape;
    Shape block_shape(2);
    Array l(block_shape);
    Array r(block_shape);
    
    // give each element a unique value
    int index = 0;
    for(int& i : l)
    {
        i = index;
        ++index;
    }
    for(int& i : r)
    {
        i = index;
        ++index;
    }
    typedef vector<pair<int, int> > vec;
    vec visited_values;
    auto visitor = [&visited_values](int left, int right)
    {
        visited_values.push_back({left, right});
    };

    Shape difference;
    NeighborhoodType neighborhood;
    auto visit = [&]()
    {
        visited_values.clear();
        visitBorder(multi_range(l), multi_range(r),
                    difference, neighborhood, visitor);
    };

    difference = Shape(1, 0);
    neighborhood = DirectNeighborhood;
    visit();
    vigra_assert((visited_values == vec{{l(1,0), r(0,0)}, {l(1,1), r(0,1)}}), "");
    
    difference = Shape(0, 1);
    visit();
    vigra_assert((visited_values == vec{{l(0,1), r(0,0)}, {l(1,1), r(1,0)}}), "");
    
    difference = Shape(0, 1);
    neighborhood = IndirectNeighborhood;
    visit();
    vigra_assert((visited_values == vec{{l(0,1), r(0,0)}, {l(0,1), r(1,0)},
                                        {l(1,1), r(0,0)}, {l(1,1), r(1,0)}}), "");
}

