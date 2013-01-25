#include <vigra/multi_array.hxx>
#include "vigra/impex.hxx"
#include <vigra/multi_math.hxx>


int main (int argc, char ** argv) {

    // read image given as first argument
    vigra::ImageImportInfo info(argv[1]);

    // process grayscale image
    if (info.isGrayscale()) {
        // use vigra::multi_math functions
        using namespace vigra::multi_math;
        // instantiate array for image data
        vigra::MultiArray<2, unsigned char> imageArray(info.shape());
        // copy image data into array
        importImage(info, destImage(imageArray));
        // invert image
        imageArray = 255-imageArray;
        // write image data to the file given as second argument
        exportImage(srcImageRange(imageArray), vigra::ImageExportInfo(argv[2]));
    }

    // process color image
    else {
        // use vigra::multi_math functions
        using namespace vigra::multi_math;
        // instantiate array for image data
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > imageArray(info.shape());
        // copy image data into array
        importImage(info, destImage(imageArray));
        // invert image
        vigra::RGBValue<unsigned char> temp(255,255,255);
        imageArray = temp - imageArray;
        // write image data to the file given as second argument
        exportImage(srcImageRange(imageArray), vigra::ImageExportInfo(argv[2]));
    }
    return 0;
}
