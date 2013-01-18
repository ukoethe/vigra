#include <vigra/multi_array.hxx>
#include <iostream>

using namespace vigra;
int main (int argc, char ** argv) {

    MultiArray<2, int> array(Shape2(4,4));
    array = 1;
    for (int i = 0; i< array.size(); i++) {
        std::cout << array[i];
    }
    std::cout << std::endl;
    MultiArrayView<2,int, StridedArrayTag> view;
    view = array.subarray(Shape2(0,0), Shape2(2,2));
    view = 2;
    for (int i = 0; i< array.size(); i++) {
        std::cout << array[i];
    }
    std::cout << std::endl;
    view = array.subarray(Shape2(2,2), Shape2(4,4));
    view = 4;
    for (int i = 0; i< array.size(); i++) {
        std::cout << array[i];
    }
    std::cout << std::endl;
    return 0;
}
