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
    typedef ChunkedArrayLazy<2, int> Array;
    typedef ChunkedArrayLazy<2, unsigned short> DirectionsArray;
    typedef Array::shape_type Shape;

    Shape shape(4);
    Shape block_shape(2);
    
    Array data(shape, block_shape);
    Array labels(shape, block_shape);
    
    MultiArray<2, int> oldschool_data_array(shape);


    oldschool_data_array(0,0) = 0;
    oldschool_data_array(1,0) = 1;
    oldschool_data_array(2,0) = 2;
    oldschool_data_array(3,0) = 3;
    oldschool_data_array(0,1) = 10;
    oldschool_data_array(1,1) = 10;
    oldschool_data_array(2,1) = 10;
    oldschool_data_array(3,1) = 10;
    
    data.commitSubarray(Shape(0), oldschool_data_array);

    DirectionsArray directions(shape, block_shape);
    unionFindWatershedBlockwise(data, labels, DirectNeighborhood, directions);

    for(auto it = directions.begin(); it != directions.end(); ++it)
        cout << *it << " ";
    cout << endl;
    for(auto it = labels.begin(); it != labels.end(); ++it)
        cout << *it << " ";
}
