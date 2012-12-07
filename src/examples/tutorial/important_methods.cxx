#include <vigra/multi_array.hxx>
#include <iostream>

using namespace vigra;

int main (int argc, char ** argv) {
    // initialize 2x3-matrix and declare iterator
    vigra::MultiArray<2, int> intArray(Shape2(3,2));
    intArray.init(0);
    vigra::MultiArray<2, int>::iterator iter;

    // set 2nd row (equivalent to dimension 1 and index 1) to 5
    intArray.bind<1>(1) = 5;

    // print the array on console
    for (iter = intArray.begin(); iter != intArray.end(); iter++) {
		std::cout << *iter << "  ";
	}
    std::cout << std::endl;

    // initialize new array with 3rd column of intArray
    vigra::MultiArray<1, int> lowArray = intArray.bind<0>(2);

    // copy the upper-left subsquare of a 4x4-matrix
    vigra::MultiArray<2, int> _44Matrix(Shape2(4,4));
    _44Matrix.init(1);
    vigra::MultiArray<2, int> subsquare = _44Matrix.subarray(Shape2(0,0), Shape2(2,2));

    // change the elements of a subarray
    _44Matrix.subarray(Shape2(0,2),Shape2(4,4)) = 0;
    for (iter = _44Matrix.begin(); iter != _44Matrix.end(); iter++) {
		std::cout << *iter << "  ";
	}
    std::cout << std::endl;


    return 0;
}
