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

    for (iter = lowArray.begin(); iter != lowArray.end(); iter++) {
		std::cout << *iter << "  ";
	}
    std::cout << std::endl;


    return 0;
}
