#include <iostream>
#include "vigra/impex.hxx"

using namespace vigra; 


int main(int argc, char ** argv) {
    
    // read image given as first argument
    vigra::ImageImportInfo imageInfo(argv[1]);

    // instantiate array for image data
    vigra:MultiArray<2, unsigned char> imageArray;

    // copy image data into array
    importImage(info, destImage(imageArray));

    // to be continued...
}
