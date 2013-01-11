#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include "vigra/impex.hxx"

using namespace vigra; 


int main(int argc, char ** argv) {
    
    // instantiate array for image data of size 180x180 pixels
    vigra:MultiArray<2, unsigned char> imageArray(Shape2(180,180));

    // set all pixels to a certain gray value
    imageArray = 120;

    // write Image data to "testimage.gif"
    exportImage(srcImageRange(imageArray), vigra::ImageExportInfo("testimage.gif"));

    return 0;
}
