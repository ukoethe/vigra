#include <vigra/blockwise_labeling.hxx>
#include <vigra/blockify.hxx>

using namespace vigra;

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
        
        int l = labelMultiArrayBlockwise(data, labels);
    }
};

int main()
{
}
