#define VIGRA_CHECK_BOUNDS

#include <vigra/blockwise_watersheds.hxx>
#include <vigra/blockwise_watersheds.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>

#include <iostream>

using namespace std;
using namespace vigra;
using namespace blockwise_watersheds_detail;

int main()
{
    typedef MultiArray<2, int> Array;
    typedef MultiArray<2, int> DirectionsArray;
    typedef Array::difference_type Shape;

    Shape shape(4);
    Shape block_shape(2);
    
    Array data(shape);
    Array labels(shape);

    data(0,0) = 0;
    data(1,0) = 1;
    data(2,0) = 2;
    data(3,0) = 3;
    data(0,1) = 10;
    data(1,1) = 10;
    data(2,1) = 10;
    data(3,1) = 10;

    DirectionsArray directions(shape);
    MultiArray<2, DirectionsArray::view_type> directions_blocks = blockify(directions, block_shape);
    prepareBlockwiseWatersheds(OverlapsGenerator<2, int, StridedArrayTag>(data, block_shape), directions_blocks.begin(), DirectNeighborhood);
    unionFindWatershedBlockwise(data, labels, DirectNeighborhood, block_shape);

    for(auto it = directions.begin(); it != directions.end(); ++it)
        cout << *it << " ";
    cout << endl;
    for(auto it = labels.begin(); it != labels.end(); ++it)
        cout << *it << " ";
}
