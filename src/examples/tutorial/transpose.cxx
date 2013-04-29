#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/linear_algebra.hxx>

using namespace vigra;

template <unsigned int N, class T, class Stride>
void print(MultiArrayView<N, T, Stride> array) 
{
    typedef typename MultiArrayView<N, T, Stride>::iterator iterator;
    
    for (iterator i = array.begin(); i != array.end(); ++i) {
        std::cout << *i << " ";
    }
    std::cout << std::endl;
}

int main (int argc, char ** argv) {

    // initialize array
    MultiArray<2, int> intarray(Shape2(4,4));
    for (int i = 0; i < intarray.size(); i++) {
        intarray[i] = i % intarray.shape(0);
    }

    std::cout << "intarray:\n";
    print(intarray);

    // create a transposed array and a transposed view
    MultiArray<2, int> transarray = intarray.transpose();
    MultiArrayView<2, int> transarrayView = intarray.transpose();

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
    // instantiate 5D array
    MultiArray<5, int> array5D(Shape5(1,2,3,4,5));

    // print the shape of the original array
    std::cout << "Shape of Array5D:\n";
    for (int i = 0; i < 5; i++) {
        std::cout << array5D.shape(i);
    }
    std::cout << std::endl;

    // transpose array
    MultiArrayView<5, int> arrayview5D = array5D.transpose();

    // print the shape of transposed array
    std::cout << "Shape of transposed Array5D:\n";
    for (int i = 0; i < 5; i++) {
        std::cout << arrayview5D.shape(i);
    }
    std::cout << std::endl;
    
    // transpose to an explicitly specified axis permutation
    MultiArrayView<5, int> arrayview5D_permuted = array5D.permuteDimensions(Shape5(2,1,3,4,0));

    // print the shape of transposed array
    std::cout << "Shape of transposed Array5D:\n";
    for (int i = 0; i < 5; i++) {
        std::cout << arrayview5D_permuted.shape(i);
    }
    std::cout << std::endl;
    

    return 0;
}




