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

	// declare iterator
	vigra::MultiArray<2, int>::iterator iter;
	int count = 0;
	// iterate over intArray, set a new value and print the element
	for (iter = intArray.begin(); iter != intArray.end(); iter++) {
		*iter = ++count;
		std::cout << *iter << std::endl;
	}

	// section: add two 3x3-matrices
	vigra::MultiArray<2, int> matrix1(Shape2(3,3));
	vigra::MultiArray<2, int> matrix2(Shape2(3,3));
	matrix1.init(1);
	matrix2.init(3);
	
	int elements = matrix1.size();			// number of array elements
	for (int i=0; i < elements; i++) {
		matrix1[i] += matrix2[i];
	}

	// section: wrong math (add 2x3- and 3x2-matrices)
	vigra::MultiArray<2, int> _23Matrix(Shape2(2,3));
	vigra::MultiArray<2, int> _32Matrix(Shape2(3,2));
	_23Matrix.init(1);
	_32Matrix.init(3);
	
	int size = _23Matrix.size();			// number of array elements
	for (int i=0; i < size; i++) {
		_23Matrix[i] += _32Matrix[i];
	}

    // section: correct addition of 2 matrices
    matrix1 += matrix2;                     // works fine!
/*  _23Matrix += _32Matrix;                 // error: wrong matrix sizes! */

    // indexing via coordinates
    // print element in second column and third row
    std::cout << matrix1(1,2) << std::endl;

    // set all elements of second row to 13
    for (int i = 0; i < _23Matrix.size(1); i++) {
        _23Matrix(1,i) = 13;
    }

    // access via Shape2(x,y), set element in second column and third row
    matrix1[Shape2(1,2)] = 22;

	return 0;
}
