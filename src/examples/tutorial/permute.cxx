#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/linear_algebra.hxx>

using namespace vigra;

int main (int argc, char ** argv) {

    // initialize array
    vigra::MultiArray<4, int> intarray(Shape4(1,2,3,4));
    
    // print the shape of the original array
    std::cout << "Shape of intarray:\n";
    for (int i = 0; i < intarray.shape().size(); i++) {
        std::cout << intarray.size(i);
    }
    std::cout << std::endl;

    // permute dimensions
    vigra::MultiArray<4, int> permuted = intarray.permuteDimensions(Shape4(1,2,3,0));

    // print the shape of the permuted array
    std::cout << "Shape of permuted intarray:\n";
    for (int i = 0; i < permuted.shape().size(); i++) {
        std::cout << permuted.size(i);
    }
    std::cout << std::endl;

    return 0;
}
