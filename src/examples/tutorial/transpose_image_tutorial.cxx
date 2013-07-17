#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/basicgeometry.hxx>

using namespace vigra; 

int main(int argc, char ** argv)
{
    if(argc != 3) 
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(grayscale images only, supported formats: " << impexListFormats() << ")" << std::endl;
        return 1;
    }
    
    // choose diagonal for transpose
    std::cout << "Transpose about which diagonal?\n";
    std::cout << "1 - major\n";
    std::cout << "2 - minor\n";
    int mode;
    std::cin >> mode;

    try 
    {
        // read image given as first command line argument
        MultiArray<2, UInt8> imageArray;
        importImage(argv[1], imageArray);

        if(mode == 1)
        {
            // when transposing about the major diagonal, we can simply
            // write a transposed view to the file given as second argument
            exportImage(imageArray.transpose(), argv[2]);
        }
        else
        {
            // when transposing about the minor diagonal, we need a new
            // image with reversed shape to hold the transposed data
            MultiArray<2, UInt8> transposed(reverse(imageArray.shape()));
            transposeImage(imageArray, transposed, minor);
            exportImage(transposed, argv[2]);
        }
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
