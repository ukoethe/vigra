#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_math.hxx>
#include <vigra/basicgeometry.hxx>

using namespace vigra;


int main (int argc, char ** argv) 
{
    // read image given as first argument
    vigra::ImageImportInfo info(argv[1]);
    
    int w = info.width(),
        h = info.height();

    // process grayscale image
    if (info.isGrayscale()) 
    {
        // instantiate arrays for image data and for composite image of appropriate size
        vigra::MultiArray<2, UInt8> imageArray(info.shape());
        vigra::MultiArray<2, UInt8> exportArray(2*info.shape());

        // copy image data into array
        importImage(info, imageArray);

        // topleft corner of composite image contains the original image
        exportArray.subarray(Shape2(0,0), info.shape()) = imageArray;

        // topright corner contains the vertically reflected image
        MultiArrayView<2, UInt8> topright = exportArray.subarray(Shape2(w,0), Shape2(2*w, h));
        reflectImage(imageArray, topright, vertical);

        // reflect the upperhalf horizontally to the lowerhalf
        MultiArrayView<2, UInt8> upperhalf = exportArray.subarray(Shape2(0,0), Shape2(2*w, h));
        MultiArrayView<2, UInt8> lowerhalf = exportArray.subarray(Shape2(0,h), 2*info.shape());
        reflectImage(upperhalf, lowerhalf, horizontal);

        // write image data to the file given as second argument
        exportImage(exportArray, ImageExportInfo(argv[2]));
    }
    // process color image
    else 
    {
        // instantiate arrays for image data and for composite image of appropriate size
        MultiArray<2, RGBValue<UInt8> > imageArray(info.shape());
        MultiArray<2, RGBValue<UInt8> > exportArray(2*info.shape());

        // copy image data into array
        importImage(info, imageArray);

        // topleft corner of composite image contains the original image
        exportArray.subarray(Shape2(0,0), info.shape()) = imageArray;

        // topright corner contains the vertically reflected image
        MultiArrayView<2, RGBValue<UInt8> > topright = exportArray.subarray(Shape2(w,0), Shape2(2*w, h));
        reflectImage(imageArray, topright, vertical);

        // reflect the upperhalf horizontally to the lowerhalf
        MultiArrayView<2, RGBValue<UInt8> > upperhalf = exportArray.subarray(Shape2(0,0), Shape2(2*w, h));
        MultiArrayView<2, RGBValue<UInt8> > lowerhalf = exportArray.subarray(Shape2(0,h), 2*info.shape());
        reflectImage(upperhalf, lowerhalf, horizontal);

        // write image data to the file given as second argument
        exportImage(exportArray, ImageExportInfo(argv[2]));
    }
    return 0;
}
