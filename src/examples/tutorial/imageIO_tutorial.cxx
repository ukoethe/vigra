#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include "vigra/impex.hxx"

using namespace vigra; 


int main(int argc, char ** argv) {
    
    // read image given as first argument
    vigra::ImageImportInfo imageInfo(argv[1]);

    // instantiate array for image data
    vigra::MultiArray<2, unsigned char> imageArray(imageInfo.shape());

    // copy image data into array
    importImage(imageInfo, destImage(imageArray));

    // set every second horizontal line to black
    for (int i = 0; i<imageInfo.height(); i+=2) {
        imageArray.bind<1>(i)=0;
    }

    // write image data to the file given as second argument
    exportImage(srcImageRange(imageArray), vigra::ImageExportInfo(argv[2]));

    return 0;
}
