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
   

    typedef GridGraph<2, undirected_tag> Graph;
    typedef typename Graph::NodeIt GraphScanner;
    typedef typename Graph::OutArcIt NeighborIterator;

    Graph graph(Shape(2), IndirectNeighborhood);
    cout << "max degree " << graph.maxDegree() << endl;
    for(GraphScanner node(graph); node != lemon::INVALID; ++node)
    {
        cout << "at " << *node << endl;
        for(NeighborIterator arc(graph, *node); arc != lemon::INVALID; ++arc)
        {
            cout << arc.neighborIndex() << endl;
        }
    }

    Shape shape(4);
    Shape block_shape(2);
    
    Array data(shape);
    Array labels(shape);
    DirectionsArray directions(shape);
    MultiArray<2, DirectionsArray::view_type> directions_blocks = blockify(directions, block_shape);
    prepareBlockwiseWatersheds(OverlapsGenerator<2, int, StridedArrayTag>(data, block_shape), directions_blocks.begin());
}
