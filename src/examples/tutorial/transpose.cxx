#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/linear_algebra.hxx>

using namespace vigra;

template <class T>
void print(MultiArrayView<2, T> array) 
{
    for(int y=0; y<array.shape(1); ++y)
    {
        for(int x=0; x<array.shape(0); ++x)
            std::cout << array(x, y) << "  ";
        std::cout << "\n";
	}
}

int main (int argc, char ** argv)
{
    // create array
    MultiArray<2, int> base_array(Shape2(4,4));
    
    // init array such that pixel values are equal to their x coordinate
    for (int i = 0; i < base_array.size(); i++) 
    {
        base_array[i] = i % base_array.shape(0);
    }

    std::cout << "base_array:\n";
    print(base_array);

    // create a transposed array and a transposed view
    MultiArray<2, int> transarray = base_array.transpose();
    MultiArrayView<2, int> transarrayView = base_array.transpose();

    std::cout << "transarray:\n";
    print(transarray);
    std::cout << "transArrayView:\n";
    print(transarrayView);

    // set transarray to 5    
    transarray = 5;
    std::cout << "base_array after setting transarray to 5\n(no change, since transarray is a copy):\n";
    print(base_array);

    // set transarrayView to 5
    transarrayView = 5;
    std::cout << "base_array after setting transarrayView to 5\n(base_array changes because transarrayView is a view):\n";
    print(base_array);

    // transposing a 5D array
    // instantiate 5D array
    MultiArray<5, int> array5D(Shape5(1,2,3,4,5));

    // print the shape of the original array
    std::cout << "Shape of array5D: " << array5D.shape() << "\n";

    // transpose array
    MultiArrayView<5, int> arrayview5D = array5D.transpose();

    // print the shape of transposed array
    std::cout << "Shape of array5D view after default transpose(): " << arrayview5D.shape() << "\n";
    
    // transpose to an explicitly specified axis permutation
    MultiArrayView<5, int> arrayview5D_permuted = array5D.transpose(Shape5(2,1,3,4,0));

    // print the shape of transposed array
    std::cout << "Shape of array5D view after user-defined transpose(): " << arrayview5D_permuted.shape() << "\n";
    std::cout << "    (applied permutation 2 => 0, 1 => 1, 3 => 2, 4 => 3, 0 => 4 to the axes)\n";

    return 0;
}




