#include <vigra/multi_array.hxx>
#include <iostream>

using namespace vigra;

int main (int argc, char ** argv) {
	
	// declare 2-dimensional MultiArray with 3 integer elements in first dimension and 2 
	// integer elements in second dimension
	vigra::MultiArray<2, int> intArray(Shape2(3,2));
	
	// set all elements on 3
	intArray.init(3);

	// print the first element on console
	std::cout << intArray[0] << std::endl;

	return 0;
}
