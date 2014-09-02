#include <vigra/blockwise_labeling.hxx>
#include <vigra/blockify.hxx>

using namespace vigra;
using namespace blockwise_labeling_detail;

class BlockwiseLabelingTest
{
public:
    void debugTest()
    {
        typedef MultiArray<2, int> Array;
        typedef MultiArray<2, Array::view_type> Blocks;
        typedef Array::difference_type Shape;
        
        Shape block_shape(2);

        Array data(Shape(4));
        Blocks data_blocks = blockify(data, block_shape);

        Array labels(Shape(4));
        Blocks label_blocks = blockify(labels, block_shape);
        
        //label_blocks.begin()[Shape(1,1)];
        /*
        Shape blocks_shape = data_blocks.shape();
        int l = blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                  label_blocks.begin(), label_blocks.end(),
                                  DirectNeighborhood, std::equal_to<int>(), blocks_shape);
        */
    }
};

int main()
{
}
