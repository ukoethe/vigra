#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include "vigra/impex.hxx"

using namespace vigra; 


int main(int argc, char ** argv) {
    
    // instantiate array for image data of size 160x160 pixels
    vigra:MultiArray<2, unsigned char> imageArray(Shape2(160,160));

    // set the pixels to black (0) or white (255)
        for (int i = 0; i < imageArray.size(0); i++) {
        for (int j = 0; j < imageArray.size(1); j++) {

        if ((i%20)/10 == (j%20)/10) imageArray[Shape2(i,j)] = 0;
        else                        imageArray[Shape2(i,j)] = 255; 
        }
    }

    // write Image data to "testimage.gif"
    exportImage(srcImageRange(imageArray), vigra::ImageExportInfo("testimage.gif"));

    return 0;
}
