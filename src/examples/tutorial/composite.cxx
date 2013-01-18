#include <vigra/multi_array.hxx>
#include "vigra/impex.hxx"
#include <vigra/multi_math.hxx>
#include <vigra/basicgeometry.hxx>

using namespace vigra;


int main (int argc, char ** argv) {

    // read image given as first argument
    vigra::ImageImportInfo info(argv[1]);

    // process grayscale image
    if (info.isGrayscale()) {
       
        // instantiate arrays for image data and for composite image of appropriate size
        vigra::MultiArray<2, unsigned char> imageArray(info.shape());
        vigra::MultiArray<2, unsigned char> exportArray(2*info.shape());

        // copy image data into array
        importImage(info, destImage(imageArray));

        // topleft corner of composite image contains the original image
        exportArray.subarray(Shape2(0,0), info.shape()) = imageArray;

        // topright corner contains the vertically reflected image
        vigra::MultiArrayView<2, unsigned char> topright = 
            exportArray.subarray(Shape2(info.width(),0), Shape2(2*info.width(), info.height()));
        reflectImage(srcImageRange(imageArray), destImage(topright), vigra::vertical);

        // reflect the upperhalf horizontally to the lowerhalf
        vigra::MultiArrayView<2, unsigned char> upperhalf = 
            exportArray.subarray(Shape2(0,0), Shape2(2*info.width(), info.height()));
        vigra::MultiArrayView<2, unsigned char> lowerhalf = 
            exportArray.subarray(Shape2(0,info.height()), 2*info.shape());
        reflectImage(srcImageRange(upperhalf), destImage(lowerhalf), vigra::horizontal);

        // write image data to the file given as second argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));
    }

    // process color image
    else {

        // instantiate arrays for image data and for composite image of appropriate size
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > imageArray(info.shape());
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > exportArray(2*info.shape());

        // copy image data into array
        importImage(info, destImage(imageArray));

        // topleft corner of composite image contains the original image
        exportArray.subarray(Shape2(0,0), info.shape()) = imageArray;

        // topright corner contains the vertically reflected image
        vigra::MultiArrayView<2, vigra::RGBValue<unsigned char> > topright = 
            exportArray.subarray(Shape2(info.width(),0), Shape2(2*info.width(), info.height()));
        reflectImage(srcImageRange(imageArray), destImage(topright), vigra::vertical);

        // reflect the upperhalf horizontally to the lowerhalf
        vigra::MultiArrayView<2, vigra::RGBValue<unsigned char> > upperhalf = 
            exportArray.subarray(Shape2(0,0), Shape2(2*info.width(), info.height()));
        vigra::MultiArrayView<2, vigra::RGBValue<unsigned char> > lowerhalf = 
            exportArray.subarray(Shape2(0,info.height()), 2*info.shape());
        reflectImage(srcImageRange(upperhalf), destImage(lowerhalf), vigra::horizontal);

        // write image data to the file given as second argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));
    }
    return 0;
}
