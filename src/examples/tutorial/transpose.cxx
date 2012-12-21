#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/linear_algebra.hxx>

using namespace vigra;

void print(vigra::MultiArrayView<2,int, StridedArrayTag> array) {
    vigra::MultiArrayView<2,int>::iterator iter;
    for (iter = array.begin(); iter != array.end(); iter++) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
}

void print(vigra::MultiArrayView<3,int, StridedArrayTag> array) {
    vigra::MultiArrayView<3,int>::iterator iter;
    for (iter = array.begin(); iter != array.end(); iter++) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
}

int main (int argc, char ** argv) {

    // initialize array
    vigra::MultiArray<2, int> intarray(Shape2(4,4));
    for (int i = 0; i < intarray.size(); i++) {
        intarray[i] = i % intarray.size(0);
    }

    std::cout << "intarray:\n";
    print(intarray);

    // create a transposed array and a transposed view
    vigra::MultiArray<2, int> transarray = intarray.transpose();
    vigra::MultiArrayView<2, int, StridedArrayTag> transarrayView = intarray.transpose();

    std::cout << "transarray:\n";
    print(transarray);
    std::cout << "transArrayView:\n";
    print(transarrayView);

    // set transarray to 5    
    transarray = 5;
    std::cout << "intarray after setting transarray to 5:\n";
    print(intarray);

    // set transarrayView to 5
    transarrayView = 5;
    std::cout << "intarray after setting transarrayView to 5:\n";
    print(intarray);

    // transposing a 5D array
    vigra::MultiArray<5, int> array5D(Shape5(1,2,3,4,5));

    std::cout << "Shape of Array5D:\n";
    for (int i = 0; i < array5D.shape().size(); i++) {
        std::cout << array5D.size(i);
    }
    std::cout << std::endl;

    array5D = array5D.transpose();

    std::cout << "Shape of transposed Array5D:\n";
    for (int i = 0; i < array5D.shape().size(); i++) {
        std::cout << array5D.size(i);
    }
    std::cout << std::endl;
   
    return 0;
}




