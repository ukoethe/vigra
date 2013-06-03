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
		std::cout << *iter << " ";
	}
    std::cout << std::endl;

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
	vigra::MultiArray<2, int> matrix32(Shape2(2,3));
	vigra::MultiArray<2, int> matrix23(Shape2(3,2));
	matrix32.init(1);
	matrix23.init(3);
	
	int size = matrix32.size();			// number of array elements
	for (int i=0; i < size; i++) {
		matrix32[i] += matrix23[i];
	}

    // section: correct addition of 2 matrices
    matrix1 += matrix2;                     // works fine!
/*  matrix32 += matrix23;                 // error: wrong matrix sizes! */

    // indexing via coordinates
    // print element in second column and third row
    std::cout << matrix1(1,2) << std::endl;

    // set all elements of second row to 13
    for (int i = 0; i < matrix32.size(1); i++) {
        matrix32(1,i) = 13;
    }

    // access via Shape2(x,y), set element in second column and third row
    matrix1[Shape2(1,2)] = 22;

    // iterating over a Shape-Object
    vigra::MultiArray<2, int> matrix84(Shape2(4,8));
    matrix84.init(5);
    // instantiate Shape-Object
    Shape2 p;
    // iterate over 3rd column :
    // set first dimension on 2 (equals 3rd column)
    // then iterate over second dimension (equals rows)
    p[0] = 2;                                  
    for(p[1]=0; p[1]<matrix84.size(1); p[1]++) {   
        matrix84[p] = 7;
    }

    // flatten an array in scan order
    // set first row of matrix23 to 1 2 3, second row to 4 5 6
    count = 0;
    for (iter = matrix23.begin(); iter != matrix23.end(); iter++) {
        *iter = count++;
    }
    // create 1D-array of appropriate size
    vigra::MultiArray<1, int> flatArray(Shape1(matrix23.size()));
    // copy 2D-array into 1D-array
    std::copy(matrix23.begin(), matrix23.end(), flatArray.begin());
    // print 1D-array on console; 
    for (iter = flatArray.begin(); iter != flatArray.end(); iter++) {
         std::cout << *iter << " ";
    }
    std::cout << std::endl;
	return 0;
}
